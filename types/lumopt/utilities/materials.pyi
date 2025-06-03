from _typeshed import Incomplete
from lumopt.utilities.wavelengths import Wavelengths as Wavelengths

class Material:
    object_dielectric: Incomplete
    base_epsilon: Incomplete
    name: Incomplete
    mesh_order: Incomplete
    def __init__(self, base_epsilon: float = 1.0, name=..., mesh_order: Incomplete | None = None) -> None: ...
    wavelengths: Incomplete
    permittivity: Incomplete
    def set_script(self, sim, poly_name) -> None: ...
    def get_eps(self, wavelengths): ...
    @staticmethod
    def get_wavelengths(sim): ...
