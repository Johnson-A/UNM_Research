from dolfin import *

class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[0], 0) and on_boundary)

    def map(self, x, y):
        y[0] = x[0] - 1.0
        y[1] = x[1]


mesh = UnitSquareMesh(10, 10)
pbc = PeriodicBoundary()

V = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=pbc)
S = FunctionSpace(mesh, 'CG', 1, constrained_domain=pbc)

fluid_vel = Function(V)

v_t = TestFunction(V)
s_t = TestFunction(S)

r_v = inner(grad(v_t), grad(fluid_vel)) * dx

def top(x, on_boundary): return near(x[1], 1.0)
def bottom(x, on_boundary): return near(x[1], 0.0)

bcv1 = DirichletBC(V, (1.0, 1.0), top)
bcv0 = DirichletBC(V, (0.0, 1.0), bottom)

solve(r_v == 0.0, fluid_vel, [bcv0, bcv1])

File('output.pvd') << fluid_vel

print(type(fluid_vel))
