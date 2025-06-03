from _typeshed import Incomplete
from lumopt.optimizers.maximizer import Maximizer as Maximizer

class AdaptiveGradientDescent(Maximizer):
    max_dx: Incomplete
    all_params_equal: Incomplete
    predictedchange_hist: Incomplete
    min_dx: Incomplete
    dx_regrowth_factor: Incomplete
    dx: Incomplete
    def __init__(self, max_dx, min_dx, max_iter, dx_regrowth_factor, all_params_equal, scaling_factor) -> None: ...
    current_params: Incomplete
    current_fom: Incomplete
    def run(self): ...
    def calculate_change(self, gradients, dx): ...
    def reduce_step_size(self) -> None: ...
    def enforce_bounds(self, params): ...
