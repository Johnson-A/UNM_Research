from dolfin import *

mesh = UnitSquareMesh(8,8)

V = VectorFunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)

# Define boundary conditions
tol = 1E-14
def lower_boundary(x, on_boundary):
    return on_boundary and abs(x[1]) < tol

def upper_boundary(x, on_boundary):
    return on_boundary and abs(x[1]-1) < tol

Gamma_0 = DirichletBC(V, (0.0,0.0), lower_boundary)
Gamma_1 = DirichletBC(V, (1.0,0.0), upper_boundary)
bcs = [Gamma_0, Gamma_1]

a = div(v) * div(u) * dx

L = 0.0

# Compute solution
t = Function(V)

solve(a==L,t,bcs)

plot(u, titel='Velocity', interactive=True)
