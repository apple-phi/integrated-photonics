from _typeshed import Incomplete

class Simulation:
    fdtd: Incomplete
    workingDir: Incomplete
    def __init__(self, workingDir, use_var_fdtd, hide_fdtd_cad) -> None: ...
    def save(self, name): ...
    def load(self, name) -> None: ...
    def save_index_to_vtk(self, filename) -> None: ...
    def save_fields_to_vtk(self, filename) -> None: ...
    def remove_data_and_save(self) -> None: ...
    def __del__(self) -> None: ...
