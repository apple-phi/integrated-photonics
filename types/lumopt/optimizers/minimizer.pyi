from _typeshed import Incomplete
from lumopt.optimizers.optimizer import Optimizer as Optimizer

class Minimizer(Optimizer):
    current_fom: Incomplete
    current_params: Incomplete
    current_gradients: Incomplete
    def define_callables(self, callable_fom, callable_jac): ...
