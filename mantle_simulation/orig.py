"""
This version of the code runs a swarm of simulations of various viscosities
and temperatures per viscosity. Since 5-29, it also includes an adiabatic
temperature variance at the LAB.
"""

from __future__ import print_function

import errno
import itertools
import math
import shutil
import sys
from collections import defaultdict
from functools import partial
from os import makedirs
from time import (clock, strftime)

from dolfin import (Constant, DirichletBC, ERROR, Expression, Function,
                    FunctionSpace, MPI, MixedFunctionSpace, Point, RectangleMesh, SubDomain,
                    TestFunctions, VectorFunctionSpace, XDMFFile, assign, div, dot, dx, exp,
                    grad, inner, interpolate, mpi_comm_world, near, project, set_log_level,
                    solve, split, sym, tanh)

set_log_level(ERROR)

COMM = mpi_comm_world()
RANK = MPI.rank(COMM)
IS_MAIN_PROC = RANK == 0

rho_0 = 3300.0
rhomelt = 2900.0
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

output_every = 10
# Number of cells in each dimension
nx = 30
ny = nx

# Non-dimensional mesh size
mesh_width = 1.0
mesh_height = 0.4 * mesh_width

def main_proc(work):
    def no_op(*_, **__): pass

    return work if IS_MAIN_PROC else no_op


log = main_proc(print)


@main_proc
def time_left(steps_finished, total_steps, start_time):
    completed = steps_finished / total_steps
    rate = completed / (clock() - start_time)
    seconds_left = (1.0 - completed) / rate if steps_finished != 0 else 0.0
    print('{0:.4f} {1:.2f}'.format(completed, seconds_left / 60.0))


class LithosExp(Expression):
    height = 0.05
    width = 0.2
    scale = 0.05
    LAB_height = 0.75 * mesh_height

    def ridge(self, x, offset):
        return self.height * (1.0 - tanh((x[0] - (0.5 + offset) * mesh_width) / self.scale))

    def eval(self, values, x):
        hump = self.ridge(x, self.width) - self.ridge(x, -self.width)
        values[0] = self.LAB_height - hump


def linear_interpolant(x0, y0, x1, y1, x_val):
    return y0 + (y1 - y0) / (x1 - x0) * (x_val - x0)


class TemperatureProfile(Expression):
    LAB = LithosExp()
    dy = mesh_height / ny

    def __init__(self, temperatures):
        self.delta = max(temperatures) - min(temperatures)
        temperatures = [x / self.delta for x in temperatures]  # Make temperature non-dimensional

        self.surface = temperatures[0]
        self.lithosphere_lab = temperatures[1]
        self.asthenosphere_lab = temperatures[2]
        self.bottom = temperatures[3]

    def temperature(self, x):
        lab_height = self.LAB(x)

        if x[1] >= lab_height:
            return linear_interpolant(lab_height, self.lithosphere_lab, mesh_height, self.surface, x[1])
        else:
            return linear_interpolant(0.0, self.bottom, lab_height, self.asthenosphere_lab, x[1])

    def eval(self, value, x):
        offsets = [-2.0, -1.0, 0.0, 1.0, 2.0]
        samples = [self.temperature(x + (0.0, delta_y * self.dy)) for delta_y in offsets]
        value[0] = sum(samples) / len(samples)


class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return left(x, on_boundary)

    def map(self, x, y):
        y[0] = x[0] - mesh_width
        y[1] = x[1]


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
    return on_boundary and near(x[1], mesh_height)


def bottom(x, on_boundary):
    return on_boundary and near(x[1], 0.0)


def left(x, on_boundary):
    return on_boundary and near(x[0], 0.0)


def right(x, on_boundary):
    return on_boundary and near(x[0], mesh_width)


def run_with_params(Tb, mu_value, k_s, path):
    run_time_init = clock()

    temperature_vals = [27.0 + 273, Tb + 273, 1300.0 + 273, 1305.0 + 273]
    temp_prof = TemperatureProfile(temperature_vals)

    mu_a = mu_value  # this was taken from the Blankenbach paper, can change

    Ep = b / temp_prof.delta

    mu_bot = exp(-Ep * (temp_prof.bottom * temp_prof.delta - 1573.0) + cc) * mu_a

    Ra = rho_0 * alpha * g * temp_prof.delta * h ** 3 / (kappa_0 * mu_a)
    w0 = rho_0 * alpha * g * temp_prof.delta * h ** 2 / mu_a
    tau = h / w0
    p0 = mu_a * w0 / h

    log(mu_a, mu_bot, Ra, w0, p0)

    vslipx = 1.6e-09 / w0
    vslip = Constant((vslipx, 0.0))  # Non-dimensional
    noslip = Constant((0.0, 0.0))

    time_step = 3.0E11 / tau / 10.0

    dt = Constant(time_step)
    tEnd = 3.0E15 / tau / 5.0  # Non-dimensionalising times

    mesh = RectangleMesh(Point(0.0, 0.0), Point(mesh_width, mesh_height), nx, ny)

    pbc = PeriodicBoundary()
    W = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=pbc)
    S = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)
    WSSS = MixedFunctionSpace([W, S, S, S])  # WSSS -> W

    u = Function(WSSS)

    # Instead of TrialFunctions, we use split(u) for our non-linear problem
    v, p, T, Tf = split(u)
    v_t, p_t, T_t, Tf_t = TestFunctions(WSSS)

    T0 = interpolate(temp_prof, S)

    FluidTemp = Expression('max(T0, 1.031)', T0=T0)

    muExp = Expression('exp(-Ep * (T_val * dTemp - 1573.0) + cc * x[1] / mesh_height)',
                       Ep=Ep, dTemp=temp_prof.delta, cc=cc, mesh_height=mesh_height, T_val=T0)

    Tf0 = interpolate(temp_prof, S)

    mu = Function(S)
    v0 = Function(W)

    v_theta = (1.0 - theta) * v0 + theta * v

    T_theta = (1.0 - theta) * T0 + theta * T

    Tf_theta = (1.0 - theta) * Tf0 + theta * Tf

    r_v = (inner(sym(grad(v_t)), 2.0 * mu * sym(grad(v)))
           - div(v_t) * p
           - T * v_t[1]) * dx

    r_p = p_t * div(v) * dx

    heat_transfer = k_s * (Tf_theta - T_theta) * dt

    r_T = (T_t * ((T - T0) + dt * inner(v_theta, grad(T_theta)))
           + (dt / Ra) * inner(grad(T_t), grad(T_theta))
           - T_t * heat_transfer) * dx

    # yvec = Constant((0.0, 1.0))
    # rhosolid = rho_0 * (1.0 - alpha * (T_theta * temp_prof.delta - 1573.0))
    # deltarho = rhosolid - rhomelt
    # v_f = v_theta - darcy * (grad(p) * p0 / h - deltarho * yvec * g) / w0

    v_melt = Function(W)
    yvec = Constant((0.0, 1.0))

    # TODO: inner -> dot, take out Tf_t
    r_Tf = (Tf_t * ((Tf - Tf0) + dt * dot(v_melt, grad(Tf_theta)))
            + Tf_t * heat_transfer) * dx

    r = r_v + r_p + r_T + r_Tf

    bcv0 = DirichletBC(WSSS.sub(0), noslip, top)
    bcv1 = DirichletBC(WSSS.sub(0), vslip, bottom)
    bcp0 = DirichletBC(WSSS.sub(1), Constant(0.0), bottom)
    bct0 = DirichletBC(WSSS.sub(2), Constant(temp_prof.surface), top)
    bct1 = DirichletBC(WSSS.sub(2), Constant(temp_prof.bottom), bottom)
    bctf1 = DirichletBC(WSSS.sub(3), Constant(temp_prof.bottom), bottom)

    bcs = [bcv0, bcv1, bcp0, bct0, bct1, bctf1]

    t = 0
    count = 0
    files = DefaultDictByKey(partial(create_xdmf, path))

    while t < tEnd:
        mu.interpolate(muExp)
        rhosolid = rho_0 * (1.0 - alpha * (T0 * temp_prof.delta - 1573.0))
        deltarho = rhosolid - rhomelt
        assign(v_melt, project(v0 - darcy * (grad(p) * p0 / h - deltarho * yvec * g) / w0, W))
        # use nP after to avoid projection?
        # pdb.set_trace()

        solve(r == 0, u, bcs)
        nV, nP, nT, nTf = u.split()

        if count % output_every == 0:
            time_left(count, tEnd / time_step, run_time_init)

            # TODO: Make sure all writes are to the same function for each time step
            files['T_fluid'].write(nTf)
            files['p'].write(nP)
            files['v_solid'].write(nV)
            files['T_solid'].write(nT)
            files['mu'].write(mu)
            files['v_melt'].write(v_melt)
            files['gradp'].write(project(grad(nP), W))
            files['rho'].write(project(rhosolid, S))
            files['Tf_grad'].write(project(grad(Tf), W))
            files['advect'].write(project(dt * dot(v_melt, grad(nTf))))
            files['ht'].write(project(heat_transfer, S))

        assign(T0, nT)
        assign(v0, nV)
        assign(Tf0, nTf)

        t += time_step
        count += 1

    log('Case mu=%g, Tb=%g complete. Run time = %g s' % (mu_a, Tb, clock() - run_time_init))


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else 'run'

    setup_base_directory(base)

    T_vals = [1300.0]
    mu_vals = [5e21]
    # k_s     = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0.0]
    k_s = [2e-2]

    for (mu, T, k) in itertools.product(mu_vals, T_vals, k_s):
        sub_dir = base + '/mu={0}/T={1}/k={2}'.format(mu, T, k)

        log('Creating ' + sub_dir)
        main_proc(makedirs)(sub_dir)

        COMM.barrier()
        run_with_params(T, mu, k, sub_dir)


@main_proc
def setup_base_directory(base):
    try:
        makedirs(base)
        shutil.copyfile(__file__, base + '/code_copy.py')
    except OSError as err:
        if err.errno == errno.EEXIST:
            print('Could not create base directory because it already exists')
            current_time = strftime('%Y-%m-%d_%H:%M:%S')
            backup_location = 'backup_runs/' + base + '_copy_' + current_time
            shutil.move(base, backup_location)
            print('Moved {0}/ to backup location at {1}/'.format(base, backup_location))
        else:
            print('E: ABORT Could not setup base environment')
            raise


if __name__ == '__main__':
    main()
