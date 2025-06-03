from _typeshed import Incomplete
from lumopt.optimizers.minimizer import Minimizer as Minimizer

class ScipyOptimizers(Minimizer):
    method: Incomplete
    pgtol: Incomplete
    ftol: Incomplete
    def __init__(self, max_iter, method: str = 'L-BFGS-B', scaling_factor: Incomplete | None = None, pgtol: float = 1e-05, ftol: float = 1e-12, scale_initial_gradient_to: int = 0, penalty_fun: Incomplete | None = None, penalty_jac: Incomplete | None = None) -> None: ...
    def run(self): ...
    def concurrent_adjoint_solves(self): ...
