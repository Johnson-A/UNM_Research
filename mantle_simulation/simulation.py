"""
This version of the code runs a swarm of simulations of various viscosities
and temperatures per viscosity. Since 5-29, it also includes an adiabatic
temperature variance at the LAB.
"""

from __future__ import print_function

import argparse
import errno
import itertools
import math
import os
import shutil
import sys
from collections import defaultdict
from functools import partial
from os import makedirs
from time import (clock, strftime)

import git
from dolfin import (Constant, DirichletBC, ERROR, Expression, Function,
                    FunctionSpace, MPI, Point, SubDomain,
                    TestFunctions, XDMFFile, assign, div, dot, dx, exp,
                    grad, inner, interpolate, mpi_comm_world, near, project, set_log_level,
                    solve, split, sym, parameters, FiniteElement, VectorElement, MixedElement,
                    BoxMesh)

import LAB
from constants import mesh_width, mesh_height, nx, ny, nz

parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -march=native'
parameters['form_compiler']['representation'] = 'quadrature'
# To show parameter information; print(info(parameters, True))

set_log_level(ERROR)

COMM = mpi_comm_world()
RANK = MPI.rank(COMM)
IS_MAIN_PROC = RANK == 0

rho_0 = 3300.0
rho_melt = 2900.0
darcy = 1e-13  # k over (mu*phi) in darcy
alpha = 2.5e-5
g = 9.81
# not actually used.. (in cian's notes he is using the non dimensional
# kappa in the equations, which here i have defined as 1 (as
# kappa/kappa_0) so i think i can get away with this)
kappa = 1.0E-6

b = 12.7
cc = math.log(128)
# Ep  = 0.014057
theta = 0.5
h = 1000000.0
kappa_0 = 1.0E-6

output_every = 1


def main_proc(work):
    def no_op(*_, **__): pass

    def work_with_delayed_execution(*args):
        reduced = (f() if callable(f) else f for f in args)
        return work(*reduced)

    return work_with_delayed_execution if IS_MAIN_PROC else no_op


log = main_proc(print)


@main_proc
def time_left(steps_finished, total_steps, start_time):
    completed = steps_finished / total_steps
    rate = completed / (clock() - start_time)
    seconds_left = (1.0 - completed) / rate if steps_finished != 0 else 0.0
    print('{:.4f}% Completed, Time Left {:.2f} minutes'.format(completed, seconds_left / 60.0))


def linear_interpolant(x0, y0, x1, y1, x_val):
    return y0 + (y1 - y0) / (x1 - x0) * (x_val - x0)


class TemperatureProfile(Expression):
    dz = mesh_height / nz

    def __init__(self, temperatures, **kwargs):  # TODO: kwargs
        self.delta = max(temperatures) - min(temperatures)
        temperatures = [x / self.delta for x in temperatures]  # Make temperature non-dimensional

        self.surface = temperatures[0]
        self.lithosphere_lab = temperatures[1]
        self.asthenosphere_lab = temperatures[2]
        self.bottom = temperatures[3]

    def temperature(self, x):
        lab_height_at = LAB.height_at(x)

        if x[2] >= lab_height_at:
            return linear_interpolant(lab_height_at, self.lithosphere_lab, mesh_height, self.surface, x[2])
        else:
            return linear_interpolant(0.0, self.bottom, lab_height_at, self.asthenosphere_lab, x[2])

    def eval(self, value, x):
        offsets = [-2.0, -1.0, 0.0, 1.0, 2.0]  # TODO: reduce window size
        samples = [self.temperature(x + (0.0, 0.0, delta_z * self.dz)) for delta_z in offsets]
        value[0] = sum(samples) / len(samples)


class DefaultDictByKey(defaultdict):
    def __missing__(self, key):
        self[key] = value = self.default_factory(key)
        return value


def create_xdmf(path, file_name):
    f = XDMFFile(COMM, path + '/' + file_name + '.xdmf')
    # Write out data at every step, at a small performance cost
    f.parameters['flush_output'] = True
    f.parameters['rewrite_function_mesh'] = False
    return f


def top(x, on_boundary):
    return on_boundary and near(x[2], mesh_height)


def bottom(x, on_boundary):
    return on_boundary and near(x[2], 0.0)


def left(x, on_boundary):
    return on_boundary and near(x[0], 0.0)


def right(x, on_boundary):
    return on_boundary and near(x[0], mesh_width)


def back(x, on_boundary):
    return on_boundary and near(x[1], mesh_width)


def front(x, on_boundary):
    return on_boundary and near(x[1], 0.0)


class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return left(x, on_boundary)

    def map(self, x, y):
        y[0] = x[0] - mesh_width
        y[1] = x[1]
        y[2] = x[2]

# TODO: Go over 3d conversion
def run_with_params(Tb, mu_value, k_s, path):
    run_time_init = clock()

    mesh = BoxMesh(Point(0.0, 0.0, 0.0), Point(mesh_width, mesh_width, mesh_height), nx, ny, nz)

    pbc = PeriodicBoundary()

    WE = VectorElement('CG', mesh.ufl_cell(), 2)
    SE = FiniteElement('CG', mesh.ufl_cell(), 1)
    WSSS = FunctionSpace(mesh, MixedElement(WE, SE, SE, SE), constrained_domain=pbc)
    # W = FunctionSpace(mesh, WE, constrained_domain=pbc)
    # S = FunctionSpace(mesh, SE, constrained_domain=pbc)
    W = WSSS.sub(0).collapse()
    S = WSSS.sub(1).collapse()

    temperature_vals = [27.0 + 273, Tb + 273, 1300.0 + 273, 1305.0 + 273]
    temp_prof = TemperatureProfile(temperature_vals, element=S.ufl_element())

    mu_a = mu_value  # this was taken from the Blankenbach paper, can change

    Ep = b / temp_prof.delta

    mu_bot = exp(-Ep * (temp_prof.bottom * temp_prof.delta - 1573.0) + cc) * mu_a

    # TODO: verify exponentiation
    Ra = rho_0 * alpha * g * temp_prof.delta * h ** 3 / (kappa_0 * mu_a)
    w0 = rho_0 * alpha * g * temp_prof.delta * h ** 2 / mu_a
    tau = h / w0
    p0 = mu_a * w0 / h

    log(mu_a, mu_bot, Ra, w0, p0)

    slip_vx = 1.6E-09 / w0  # Non-dimensional
    slip_velocity = Constant((slip_vx, 0.0, 0.0))
    zero_slip = Constant((0.0, 0.0, 0.0))

    time_step = 3.0E11 / tau * 2

    dt = Constant(time_step)
    t_end = 3.0E15 / tau / 5.0  # Non-dimensional times

    u = Function(WSSS)

    # Instead of TrialFunctions, we use split(u) for our non-linear problem
    v, p, T, Tf = split(u)
    v_t, p_t, T_t, Tf_t = TestFunctions(WSSS)

    T0 = interpolate(temp_prof, S)

    mu_exp = Expression('exp(-Ep * (T_val * dTemp - 1573.0) + cc * x[2] / mesh_height)',
                       Ep=Ep, dTemp=temp_prof.delta, cc=cc, mesh_height=mesh_height, T_val=T0,
                       element=S.ufl_element())

    Tf0 = interpolate(temp_prof, S)

    mu = Function(S)
    v0 = Function(W)

    v_theta = (1.0 - theta) * v0 + theta * v

    T_theta = (1.0 - theta) * T0 + theta * T

    Tf_theta = (1.0 - theta) * Tf0 + theta * Tf

    # TODO: Verify forms

    r_v = (inner(sym(grad(v_t)), 2.0 * mu * sym(grad(v)))
           - div(v_t) * p
           - T * v_t[2]) * dx

    r_p = p_t * div(v) * dx

    heat_transfer = Constant(k_s) * (Tf_theta - T_theta) * dt

    r_T = (T_t * ((T - T0) + dt * inner(v_theta, grad(T_theta)))  # TODO: Inner vs dot
           + (dt / Ra) * inner(grad(T_t), grad(T_theta))
           - T_t * heat_transfer) * dx

    v_melt = Function(W)
    z_hat = Constant((0.0, 0.0, 1.0))

    # TODO: inner -> dot, take out Tf_t
    r_Tf = (Tf_t * ((Tf - Tf0) + dt * inner(v_melt, grad(Tf_theta)))
            + Tf_t * heat_transfer) * dx

    r = r_v + r_p + r_T + r_Tf

    bcv0 = DirichletBC(WSSS.sub(0), zero_slip, top)
    bcv1 = DirichletBC(WSSS.sub(0), slip_velocity, bottom)
    bcv2 = DirichletBC(WSSS.sub(0).sub(1), Constant(0.0), back)
    bcv3 = DirichletBC(WSSS.sub(0).sub(1), Constant(0.0), front)

    bcp0 = DirichletBC(WSSS.sub(1), Constant(0.0), bottom)
    bct0 = DirichletBC(WSSS.sub(2), Constant(temp_prof.surface), top)
    bct1 = DirichletBC(WSSS.sub(2), Constant(temp_prof.bottom), bottom)
    bctf1 = DirichletBC(WSSS.sub(3), Constant(temp_prof.bottom), bottom)

    bcs = [bcv0, bcv1, bcv2, bcv3, bcp0, bct0, bct1, bctf1]

    t = 0
    count = 0
    files = DefaultDictByKey(partial(create_xdmf, path))

    while t < t_end:
        mu.interpolate(mu_exp)
        rhosolid = rho_0 * (1.0 - alpha * (T0 * temp_prof.delta - 1573.0))
        deltarho = rhosolid - rho_melt
        # TODO: project (accuracy) vs interpolate
        assign(v_melt, project(v0 - darcy * (grad(p) * p0 / h - deltarho * z_hat * g) / w0, W))
        # TODO: Written out one step later?
        # v_melt.assign(v0 - darcy * (grad(p) * p0 / h - deltarho * yvec * g) / w0)
        # TODO: use nP after to avoid projection?

        solve(r == 0, u, bcs)
        nV, nP, nT, nTf = u.split()  # TODO: write with Tf, ... etc

        if count % output_every == 0:
            time_left(count, t_end / time_step, run_time_init)  # TODO: timestep vs dt

            # TODO: Make sure all writes are to the same function for each time step
            files['T_fluid'].write(nTf, t)
            files['p'].write(nP, t)
            files['v_solid'].write(nV, t)
            files['T_solid'].write(nT, t)
            files['mu'].write(mu, t)
            files['v_melt'].write(v_melt, t)
            files['gradp'].write(project(grad(nP), W), t)
            files['rho'].write(project(rhosolid, S), t)
            files['Tf_grad'].write(project(grad(Tf), W), t)
            files['advect'].write(project(dt * dot(v_melt, grad(nTf))), t)
            files['ht'].write(project(heat_transfer, S), t)

        assign(T0, nT)
        assign(v0, nV)
        assign(Tf0, nTf)

        t += time_step
        count += 1

    log('Case mu={}, Tb={}, k={} complete. Run time = {:.2f} minutes'.format(mu_a, Tb, k_s, (clock() - run_time_init) / 60.0))


def main():
    # TODO: Argument parsing
    parser = argparse.ArgumentParser(description='Run mantle simulation')
    parser.add_argument('base')
    parser.add_argument('mu_vals', metavar='MU', type=float, nargs='+')
    base = sys.argv[1] if len(sys.argv) > 1 else 'run'

    setup_base_directory(base)

    T_vals = [1300.0]
    mu_vals = [5e21]
    k_s = [2e-2, 1e-2, 1e-3, 0.0]

    for (mu, T, k) in itertools.product(mu_vals, T_vals, k_s):
        sub_dir = base + '/mu={}/Tb={}/k={}'.format(mu, T, k)

        log('Creating ' + sub_dir)
        main_proc(makedirs)(sub_dir)

        COMM.barrier()
        run_with_params(T, mu, k, sub_dir)


@main_proc
def setup_base_directory(base):
    try:
        makedirs(base)
    except OSError as err:
        if err.errno == errno.EEXIST:
            print('Could not create base directory because it already exists')
            current_time = strftime('%Y-%m-%d_%H:%M:%S')
            backup_location = 'backup_runs/' + base + '_copy_' + current_time
            shutil.move(base, backup_location)
            print('Moved {}/ to backup location at {}/'.format(base, backup_location))
        else:
            print('E: ABORT Could not setup base environment')
            raise

    copy_dir = base + '/code_copy'
    makedirs(copy_dir)
    repo = git.Repo(search_parent_directories=True)

    with open(copy_dir + '/repo_sha', mode='w') as repo_info:
        repo_info.write(repo.head.object.hexsha + '\n')
        for d in repo.head.commit.diff(None, create_patch=True, paths='mantle_simulation'):
            repo_info.write(d.a_path + '\n')

            with open(copy_dir + '/diff_' + os.path.basename(d.a_path), mode='w') as diff_file:
                diff_file.write(d.diff)


if __name__ == '__main__':
    main()
