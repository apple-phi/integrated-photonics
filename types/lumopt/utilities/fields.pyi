from _typeshed import Incomplete
from lumopt.utilities.scipy_wrappers import wrapped_GridInterpolator as wrapped_GridInterpolator

class Fields:
    x: Incomplete
    y: Incomplete
    z: Incomplete
    E: Incomplete
    D: Incomplete
    H: Incomplete
    wl: Incomplete
    eps: Incomplete
    pointing_vect: Incomplete
    normalized: bool
    getfield: Incomplete
    geteps: Incomplete
    getDfield: Incomplete
    getHfield: Incomplete
    evals: int
    def __init__(self, x, y, z, wl, E, D, eps, H) -> None: ...
    def scale(self, dimension, factors) -> None: ...
    def make_field_interpolation_object(self, F): ...
    def plot(self, ax, title, cmap) -> None: ...
    def plot_full(self, D: bool = False, E: bool = True, eps: bool = False, H: bool = False, wl: float = 1.55e-06, original_grid: bool = True) -> None: ...
    def plot_field(self, field_func: Incomplete | None = None, original_grid: bool = True, wl: float = 1.55e-06, name: str = 'field') -> None: ...

class FieldsNoInterp(Fields):
    x: Incomplete
    y: Incomplete
    z: Incomplete
    deltas: Incomplete
    E: Incomplete
    D: Incomplete
    H: Incomplete
    wl: Incomplete
    eps: Incomplete
    pointing_vect: Incomplete
    normalized: bool
    getfield: Incomplete
    geteps: Incomplete
    getDfield: Incomplete
    getHfield: Incomplete
    evals: int
    def __init__(self, x, y, z, wl, deltas, E, D, eps, H) -> None: ...
    def make_field_interpolation_object_nointerp(self, F): ...
    def plot(self, ax, title, cmap) -> None: ...
    def scale(self, dimension, factors) -> None: ...
