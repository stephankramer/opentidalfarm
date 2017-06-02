from .base_turbine import BaseTurbine
from .controls import Controls

class SmearedTurbine(BaseTurbine):
    def __init__(self, thrust_coefficient=None,
                 thrust_curve=None,
                 diameter=20.): 


        # Initialize the base class.
        super(SmearedTurbine, self).__init__(
                thrust_coefficient=thrust_coefficient,
                thrust_curve=thrust_curve,
                diameter=diameter)
