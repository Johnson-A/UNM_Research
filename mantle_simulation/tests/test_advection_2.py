from dolfin import (Constant, DirichletBC, Expression, Function, FunctionSpace,
    MixedFunctionSpace, Point, RectangleMesh, TestFunction, TestFunctions,
    assign, dot, dx, grad, interpolate, near, plot, solve, split, TrialFunctions)


v1 = Constant((1.0, 0.0))
v2 = Constant((-1.0, 0.0))

def left(x, on_boundary):
    return on_boundary and near(x[0], 0.0)

def right(x, on_boundary):
    return on_boundary and near(x[0], 1.0)

mesh = RectangleMesh(Point(0.0, 0.0), Point(1.0, 1.0), 50, 50)
dt = 0.01

S = FunctionSpace(mesh, 'CG', 1)
W = MixedFunctionSpace([S, S])
u = Function(W)

T, Tf = split(u)
# T, Tf = TrialFunctions(W)
T_t, Tf_t = TestFunctions(W)

T0 = interpolate(Expression('1.0 / (10.0 * x[0] + 1.0)'), S)
Tf0 = interpolate(Expression('1.0 / (10.0 * (1.0 - x[0]) + 1.0)'), S)

heat_transfer = 10.0 * (Tf - T) * dt

r_1 = T_t * (
    (T - T0)
    + dt * dot(v1, grad(T))
    - heat_transfer) * dx

r_2 = Tf_t * (
    (Tf - Tf0)
    + dt * dot(v2, grad(Tf))
    + heat_transfer) * dx

r = r_1 + r_2

bc1 = DirichletBC(W.sub(0), 1.0, left)
bc2 = DirichletBC(W.sub(1), 1.0, right)

t = 0.0
while t <= 10.0:
    solve(r == 0, u, [bc1, bc2])

    nT, nTf = u.split()

    assign(T0, nT)
    assign(Tf0, nTf)

    plot(T0, mesh)
    plot(Tf0, mesh)

    t += dt
