from dolfin import (Constant, DirichletBC, Function, FunctionSpace, Point,
    RectangleMesh, TestFunction, dx, grad, inner, interactive, near, plot,
    solve, assign, interpolate, Expression, dot)

v = Constant((1.0, 0.0))

def left(x, on_boundary):
    return on_boundary and near(x[0], 0.0)

mesh = RectangleMesh(Point(0.0, 0.0), Point(1.0, 1.0), 50, 50)
dt = 0.001

S = FunctionSpace(mesh, 'CG', 1)
T = Function(S)
T0 = interpolate(Expression("1.0 / (10.0 * x[0] + 1.0)"), S)
T_t = TestFunction(S)

r = (T_t * ((T - T0) + dt * dot(v, grad(T)))) * dx
bc1 = DirichletBC(S, 1.0, left)

t = 0.0
while t <= 1.0:
    solve(r == 0, T, [bc1])

    assign(T0, T)
    plot(T, mesh)

    t += dt
