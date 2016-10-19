from dolfin import tanh

from constants import mesh_width, mesh_height

keel_width = 0.2 * mesh_width
scale = keel_width / 4

LAB_height = 0.75 * mesh_height
keel_height = mesh_height / 8


def ridge(x, offset):
    return keel_height * (1.0 - tanh((x[0] - (mesh_width / 2 + offset)) / scale))


def height_at(x):
    hump = ridge(x, keel_width) - ridge(x, -keel_width)
    return LAB_height - hump
