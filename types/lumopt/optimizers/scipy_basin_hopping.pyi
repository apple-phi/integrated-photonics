from _typeshed import Incomplete
from lumopt.optimizers.minimizer import Minimizer as Minimizer

class ScipyBasinHopping(Minimizer):
    T: Incomplete
    stepsize: Incomplete
    minimizer_kwargs: Incomplete
    take_step: Incomplete
    accept_test: Incomplete
    interval: Incomplete
    disp: Incomplete
    niter_success: Incomplete
    seed: Incomplete
    def __init__(self, niter: int = 100, T: float = 1.0, stepsize: float = 0.5, minimizer_kwargs=..., take_step: Incomplete | None = None, accept_test: Incomplete | None = None, interval: int = 50, disp: bool = True, niter_success: int = 10, seed: int = 1234567890, scaling_factor: float = 1.0, scale_initial_gradient_to: float = 0.0, penalty_fun: Incomplete | None = None, penalty_jac: Incomplete | None = None) -> None: ...
    def run(self): ...
    def report_writing(self) -> None: ...
