from _typeshed import Incomplete

class Edge:
    first_point: Incomplete
    second_point: Incomplete
    eps_in: Incomplete
    eps_out: Incomplete
    z: Incomplete
    depth: Incomplete
    normal: Incomplete
    def __init__(self, first_point, second_point, eps_in, eps_out, z, depth) -> None: ...
    def derivative(self, gradient_fields, n_points): ...
    def derivative_3D(self, gradient_fields, n_points): ...
    def derivative_2D(self, gradient_fields, n_points): ...
