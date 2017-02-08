from dolfin import tanh, sqrt

from constants import mesh_width, mesh_height

keel_width = 0.2 * mesh_width
scale = keel_width / 4

LAB_height = 0.75 * mesh_height
keel_height = mesh_height / 8


def ridge(r, offset):
    return 1.0 - tanh((r - offset) / scale)

def hump(r):
    return ridge(r, keel_width) - ridge(-r, -keel_width)

def height_at(x):
    r = sqrt((x[0] - mesh_width / 2) ** 2 + (x[1] - mesh_width / 2) ** 2)
    return LAB_height - keel_height * hump(r) / hump(0)
