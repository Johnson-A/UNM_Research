""" This version of the code runs a swarm of simulations of various viscosities and
temperatures per viscosity. Since 5-29, it also includes an adiabatic temperature variance at the LAB.
All quantities are scaled using option B in Cian_RaNotes.pdf until they are written out."""
import math
import os
import errno
from shutil import copyfile
from time import clock
from dolfin import *

# Specify optimization parameters
set_log_level(PROGRESS)
parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['representation'] = 'quadrature'

# Constant Definitions
rho_0 = 3300.0
alpha = 2.5e-5
g = 9.81
# kappa = 1.0E-6 not actually used.. (in cian's notes he is using the non dimensional kappa,
# which here i have defined as 1 (as kappa/kappa_0) so i think i can get away with this)
b = 12.7
cc = math.log(128)
theta = 0.5
h = 1000000.0
kappa_0 = 1.0E-6
outputInterval = 1

# Number of separations in each dimension
nx = 10
ny = 10
nz = 10

# Relative dimensions of the domain
meshWidth = 1.0 # defined to be 1.0 -> 1000 km roughly
meshHeight = 0.4
LABHeight = 0.75 * meshHeight # 100 km from surface


class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[0], 0) and on_boundary)

    def map(self, x, y):
        y[0] = x[0] - meshWidth
        y[1] = x[1]
        y[2] = x[2]


class LithosExp(Expression):
    def eval(self, values, x):
        height = 0.25 * meshHeight
        width = 0.2 * meshWidth # radius of protrusion
        scale = width / 4 # slope of edge of protrusion
        tanhStep = lambda radius: tanh((radius - width) / scale) + tanh((-radius - width) / scale)
        normalize = height / abs(tanhStep(0))
        r = sqrt((x[0] - meshWidth / 2) ** 2 + (x[1] - meshWidth / 2) ** 2)
        values[0] = LABHeight + normalize * tanhStep(r)


def runJob(T_b, mu_value, path):
    # Create data files -- write all data to one file?
    comm = mpi_comm_world()
    rank = MPI.rank(comm)

    def createXDMF(filePath):
        f = XDMFFile(comm, filePath)
        f.parameters['flush_output'] = True  # Write out data at every step, at a small performance cost
        f.parameters['rewrite_function_mesh'] = False  # Avoid rewriting the same function mesh
        # f.parameters['multi_file'] = 10 # issue 278 https://bitbucket.org/fenics-project/dolfin/issue/278/hdf5-file-integrity
        return f

    fileNames = ['t6t', 'mu', 'velocity', 'gradp', 'pstar']
    [tFile, muFile, uFile, gradpFile, pFile] = [createXDMF(path + name + '.xdmf') for name in fileNames]

    logFile = open(path + 'LogFile' + str(rank) + '.txt', 'w', 1) # bufSize = 1 -> line buffered
    def logText(s):
        print(s)
        logFile.write(s + '\n')

    # Calculate all values dependent on Tb and mu_value
    temps = [27.0 + 273, T_b + 273, 1300.0 + 273, 1500.0 + 273]
    dTemp = temps[3] - temps[0]
    temps = [x / dTemp for x in temps]  # non-dimensionalising temperatures

    Ep = b / dTemp
    mu_a = mu_value  # Assumption taken from the blankenbach paper

    mu_bot = exp(-Ep * (temps[3] * dTemp - 1573) + cc) * mu_a
    Ra = rho_0 * alpha * g * dTemp * h ** 3 / (kappa_0 * mu_a)
    w0 = rho_0 * alpha * g * dTemp * h ** 2 / mu_a
    tau = h / w0
    p0 = mu_a * w0 / h

    vslipx = 1.6e-09 / w0  # non-dimensional
    vslip = Constant((vslipx, 0.0, 0.0))
    noslip = Constant((0.0, 0.0, 0.0))
    dt = 3.0E11 / tau * 10 # TODO:
    tEnd = 3.0E15 / tau  # non-dimensionalising times

    # Specify Mesh and Functions
    mesh = BoxMesh(0, 0, 0, meshWidth, meshWidth, meshHeight, nx, ny, nz)

    pbc = PeriodicBoundary()
    Svel = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain = pbc)
    Spre = FunctionSpace(mesh, 'CG', 1, constrained_domain = pbc)
    Stemp = FunctionSpace(mesh, 'CG', 1, constrained_domain = pbc)
    Smu = FunctionSpace(mesh, 'CG', 1, constrained_domain = pbc)
    Sgradp = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain = pbc)

    # Create the corresponding functions, S0 is used for convenience
    S0 = MixedFunctionSpace([Svel, Spre, Stemp])
    u = Function(S0)
    du = TrialFunction(S0)
    v, p, T = split(u)
    v_t, p_t, T_t = TestFunctions(S0)

    LAB = LithosExp(element=Stemp.ufl_element())

    class TempExp(Expression):
        def eval(self, value, x):
            if x[2] >= LAB(x):
                value[0] = temps[0] + (temps[1] - temps[0]) * (meshHeight - x[2]) / (meshHeight - LAB(x))
            else:
                value[0] = temps[3] - (temps[3] - temps[2]) * x[2] / LAB(x)

    T0 = Function(Stemp, name='Temperature')
    T0.interpolate(TempExp(element=Stemp.ufl_element()))

    v0 = Function(Svel, name='Velocity')

    mu = Function(Smu, name='mu')
    muExp = Expression('exp(-Ep * (T_val * dTemp - 1573) + cc * x[2] / meshHeight)', Smu.ufl_element(),
                        Ep=Ep, dTemp=dTemp, cc=cc, meshHeight=meshHeight, T_val=T0)
    mu.interpolate(muExp)

    # Specify boundaries and boundary conditions for velocity, temperature, and pressure
    def top(x, on_boundary): return near(x[2], meshHeight)
    def bottom(x, on_boundary): return near(x[2], 0)

    def left(x, on_boundary): return (x[2] <= LAB(x) and near(x[0], 0))
    def right(x, on_boundary): return (x[2] <= LAB(x) and near(x[0], meshWidth))

    def back(x, on_boundary): return (x[2] <= LAB(x) and near(x[1], meshWidth))
    def front(x, on_boundary): return (x[2] <= LAB(x) and near(x[1], 0))

    bcv0 = DirichletBC(S0.sub(0), noslip, top)
    bcv1 = DirichletBC(S0.sub(0), vslip, bottom)
    bcv2 = DirichletBC(S0.sub(0).sub(1), Constant(0.0), back)
    bcv3 = DirichletBC(S0.sub(0).sub(1), Constant(0.0), front)

    bcp0 = DirichletBC(S0.sub(1), Constant(0.0), top)
    bct0 = DirichletBC(S0.sub(2), Constant(temps[0]), top)
    bct1 = DirichletBC(S0.sub(2), Constant(temps[3]), bottom)

    bcs = [bcv0, bcv1, bcv2, bcv3, bcp0, bct0, bct1]

    # Form definitions
    v_theta = (1.0 - theta) * v0 + theta * v

    T_theta = (1.0 - theta) * T0 + theta * T

    r_v = (inner(sym(grad(v_t)), 2.0 * mu * sym(grad(v))) - div(v_t) * p - T * v_t[2]) * dx

    r_p = p_t * div(v) * dx

    r_T = (T_t * ((T - T0) + dt * inner(v_theta, grad(T_theta))) \
        + (dt / Ra) * inner(grad(T_t), grad(T_theta))) * dx

    r = r_v + r_p + r_T

    J = derivative(r, u, du)

    problem = NonlinearVariationalProblem(r, u, bcs, J)
    solver = NonlinearVariationalSolver(problem)
    prm = solver.parameters
    prm['newton_solver']['linear_solver'] = 'mumps'
    print(info(prm, True))
    # prm['nonlinear_solver'] = 'snes'
    # prm['snes_solver']['line_search'] = 'basic'
    # prm['snes_solver']['linear_solver']= 'lu'
    # prm['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
    # prm['newton_solver']['krylov_solver']['monitor_convergence'] = True

    # Begin the simulation loop
    runtimeInit = clock()
    t = 0
    count = 0
    while t < tEnd:
        solver.solve()

        nV, nP, nT = u.split()

        assign(T0, nT)
        assign(v0, nV)
        mu.interpolate(muExp)

        t += dt
        count += 1

        if count % outputInterval == 0:
            # file << (data, timeStamp)
            pFile << (nP, t)
            uFile << (nV, t)
            tFile << (nT, t)
            muFile << (mu, t)
            gradpFile << (project(grad(nP), Sgradp), t)
            # TODO: Output to same function so we can run simulation in paraview

            #TODO: write out melt-solid relative motion vectors

            timeElapsed = clock() - runtimeInit
            rate = timeElapsed / count
            perComp = t / tEnd * 100
            tEst = (tEnd - t) / t * (clock() - runtimeInit) / 3600
            logText('Step %g' % count + ': rate = %g' % rate +
                    ' sec/step --- %g' % perComp + '%%, %g' % tEst + ' hrs')
            # TODO: Write with the scaled versions
            # Learn how to script paraview to apply common operations

    logText('Case mu=%g, Tb=%g complete.' % (mu_a, T_b) + ' Run time = %g' % (clock() - runtimeInit) + 's')
    logFile.close()


if __name__ == '__main__':
    base = 'run/'
    try:
        copyfile(__file__, base + 'code_copy.py')
    except:
        pass
    # Mus = [1e19, 1e20, 1e21]
    # Tbs = [800, 1000, 1300]
    Mus = [1e20];
    Tbs = [1300];

    for mu in Mus:
        for temp in Tbs:
            work_path = base + 'mu=' + str(mu) + '/Tb=' + str(temp) + '/'
            try:
                os.makedirs(work_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

            runJob(temp, mu, work_path)
