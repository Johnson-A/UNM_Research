'''
This version of the code runs a swarm of simulations of various viscosities
and temperatures per viscosity. Since 5-29, it also includes an adiabatic
temperature variance at the LAB.
'''

import errno
import itertools
import math
from os import makedirs
from shutil import copyfile
from time import clock

from dolfin import *

set_log_level(ERROR)

rho_0 = 3300.
rhomelt = 2900.
darcy = 1e-13  # k over (mu*phi) in darcy
alpha = 2.5e-5
g = 9.81
# not actually used.. (in cian's notes he is using the non dimensional
# kappa in the equations, which here i have defined as 1 (as
# kappa/kappa_0) so i think i can get away with this)
kappa = 1.E-6

b = 12.7
cc = math.log(128)
# Ep  = 0.014057
theta = 0.5
h = 1000000.0
kappa_0 = 1.E-6

output_every = 100
nx = 30
ny = nx
# non-dimensional mesh size

mesh_width = 1.0
mesh_height = 0.4 * mesh_width

class LithosExp(Expression):
    height = 0.05
    width = 0.2
    scale = 0.05
    LAB_height = 0.75 * mesh_height

    def ridge(self, x, offset):
        return self.height * (1 - tanh((x[0] - (0.5 + offset) * mesh_width) / self.scale))

    def eval(self, values, x):
        hump = self.ridge(x, self.width) - self.ridge(x, -self.width)
        values[0] = self.LAB_height - hump

LAB = LithosExp()

def top(x, on_boundary):
    return on_boundary and near(x[1], mesh_height)

def bottom(x, on_boundary):
    return on_boundary and near(x[1], 0.0)

def left(x, on_boundary):
    return on_boundary and near(x[0], 0.0)

def right(x, on_boundary):
    return on_boundary and near(x[0], mesh_width)

def RunJob(Tb, mu_value, k_s, path):
    runtimeInit = clock()

    def file_in(f):
        return File(path + '/' + f + '.pvd')

    T_solid_file = file_in('T_solid')
    T_fluid_file = file_in('T_fluid')
    mu_file = file_in('mu')
    v_solid_file = file_in('v_solid')
    gradp_file = file_in('gradp')
    p_file = file_in('pstar')
    v_melt_file = file_in('v_melt')
    rho_file = file_in('rho_solid')

    temp_values = [27. + 273, Tb + 273, 1300. + 273, 1305. + 273]
    dTemp = temp_values[3] - temp_values[0]
    temp_values = [x / dTemp for x in temp_values]  # non dimensionalising temp

    mu_a = mu_value  # this was taken from the blankenbach paper, can change

    Ep = b / dTemp

    mu_bot = exp(-Ep * (temp_values[3] * dTemp - 1573) + cc) * mu_a

    Ra = rho_0 * alpha * g * dTemp * h**3 / (kappa_0 * mu_a)
    w0 = rho_0 * alpha * g * dTemp * h**2 / mu_a
    tau = h / w0
    p0 = mu_a * w0 / h

    print(mu_a, mu_bot, Ra, w0, p0)

    vslipx = 1.6e-09 / w0
    vslip = Constant((vslipx, 0.0))  # nondimensional
    noslip = Constant((0.0, 0.0))

    dt = 3.E11 / tau
    tEnd = 3.E15 / tau  # non-dimensionalising times

    class PeriodicBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return left(x, on_boundary)

        def map(self, x, y):
            y[0] = x[0] - mesh_width
            y[1] = x[1]

    pbc = PeriodicBoundary()

    class TempExp(Expression):
        def eval(self, value, x):
            if x[1] >= LAB(x):
                value[0] = temp_values[0] + \
                    (temp_values[1] - temp_values[0]) * \
                    (mesh_height - x[1]) / (mesh_height - LAB(x))
            else:
                value[0] = temp_values[3] - \
                    (temp_values[3] - temp_values[2]) * (x[1]) / (LAB(x))

    mesh = RectangleMesh(Point(0.0, 0.0), Point(mesh_width, mesh_height), nx, ny)

    Svel = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=pbc)
    Spre = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)
    Stemp = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)
    Smu = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)
    Sgradp = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=pbc)
    Srho = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)
    S0 = MixedFunctionSpace([Svel, Spre, Stemp])

    u = Function(S0)
    v, p, T = split(u)

    v_t, p_t, T_t = TestFunctions(S0)

    T0 = interpolate(TempExp(), Stemp)

    FluidTemp = Expression('max(T0, 1.031)', T0=T0)

    muExp = Expression('exp(-Ep * (T_val * dTemp - 1573) + cc * x[1] / mesh_height)',
                       Ep=Ep, dTemp=dTemp, cc=cc, mesh_height=mesh_height, T_val=T0)

    mu = interpolate(muExp, Smu)

    rhosolid = Function(Srho)
    deltarho = Function(Srho)

    v0 = Function(Svel)
    vmelt = Function(Svel)

    Tf = Function(Stemp)

    v_theta = (1. - theta) * v0 + theta * v

    T_theta = (1. - theta) * T + theta * T0

    r_v = (inner(sym(grad(v_t)), 2. * mu * sym(grad(v)))
           - div(v_t) * p
           - T * v_t[1]) * dx

    r_p = p_t * div(v) * dx

    r_T = (T_t * ((T - T0) + dt * inner(v_theta, grad(T_theta)))
           + (dt / Ra) * inner(grad(T_t), grad(T_theta))
           + T_t * k_s * (Tf - T_theta) * dt) * dx

    r = r_v + r_p + r_T

    bcv0 = DirichletBC(S0.sub(0), noslip, top)
    bcv1 = DirichletBC(S0.sub(0), vslip, bottom)
    bcp0 = DirichletBC(S0.sub(1), Constant(0.0), bottom)
    bct0 = DirichletBC(S0.sub(2), Constant(temp_values[0]), top)
    bct1 = DirichletBC(S0.sub(2), Constant(temp_values[3]), bottom)

    bcs = [bcv0, bcv1, bcp0, bct0, bct1]

    t = 0
    count = 0
    while (t < tEnd):
        Tf.interpolate(FluidTemp)
        # pdb.set_trace()
        # print(Tf.vector().array())

        solve(r == 0, u, bcs)
        nV, nP, nT = u.split()

        gp = grad(nP)
        rhosolid = rho_0 * (1 - alpha * (nT * dTemp - 1573))
        deltarho = rhosolid - rhomelt
        yvec = Constant((0.0, 1.0))
        vmelt = nV * w0 - darcy * (gp * p0 / h - deltarho * yvec * g)
        # TODO vrel

        if count % output_every == 0:
            T_fluid_file << Tf
            p_file << nP
            v_solid_file << nV
            T_solid_file << nT
            mu_file << mu
            # TODO mu_file << project(mu * mu_a, Smu)
            gradp_file << project(grad(nP), Sgradp)
            rho_file << project(rhosolid, Srho)
            v_melt_file << project(vmelt, Svel)

        assign(T0, nT)
        assign(v0, nV)
        mu.interpolate(muExp)

        t += dt
        count += 1

    print('Case mu=%g, Tb=%g complete.' % (mu_a, Tb), ' Run time =', clock() - runtimeInit, 's')

if __name__ == '__main__':
    base = 'run'

    comm = mpi_comm_world()
    rank = MPI.rank(comm)

    if rank == 0:
        try:
            makedirs(base)
            copyfile(__file__, base + '/code_copy.py')
        except OSError as err:
            if err.errno == errno.EEXIST:
                print('Could not create base directory')
                # dolfin_error('Init', 'Create base dir', 'Error: Base directory already exists')

            print('Could not setup base environment')
            raise

    T_vals = [1300]
    mu_vals = [5e21]
    k_s = [1e-7, 1e-6, 1e-5, 1e-4]

    for (mu, T, k) in itertools.product(mu_vals, T_vals, k_s):
        sub_dir = base + '/mu=' + str(mu) + '/T=' + str(T) + '/k=' + str(k)

        if rank == 0:
            print('Creating ' + sub_dir)
            makedirs(sub_dir)

        comm.barrier()
        RunJob(T, mu, k, sub_dir)
