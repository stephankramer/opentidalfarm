import numpy
import dolfin
import dolfin_adjoint
import abc
from ..optimisation_helpers import MinimumDistanceConstraints
from ..optimisation_helpers import MinimumDistanceConstraintsLargeArrays

class BaseFarm(object):
    """A base Farm class from which other Farm classes should be derived."""
    def __init__(self, domain, turbine=None, site_ids=None, 
                 function_space=None, n_time_steps=None):
        """Create an empty Farm."""
        # Create a chaching object for the interpolated turbine density fields
        # (as their computation is very expensive)
        # n_time_steps is only used with dynamic thrust.

        if function_space is None:
            function_space = FunctionSpace(self.domain.mesh, "CG", 2)
        self._turbine_function_space = function_space
        self._turbine_density = Function(self._turbine_function_space)

        self.domain = domain
        self._set_turbine_specification(turbine)

        # The measure of the farm site
        self.site_dx = self.domain.dx(site_ids)


    @property
    def density_function(self):
        return self._turbine_density


    def _get_turbine_specification(self):
        if self._turbine_specification is None:
            raise ValueError("The turbine specification has not yet been set.")
        return self._turbine_specification


    def _set_turbine_specification(self, turbine_specification):
        self._turbine_specification = turbine_specification


    turbine_specification = property(_get_turbine_specification,
                                     _set_turbine_specification,
                                     "The turbine specification.")

    @property
    @abc.abstractmethod
    def control_array_global(self):
        """A serialized representation of the farm based on the controls.

        :returns: A serialized representation of the farm based on the controls.
        :rtype: numpy.ndarray
        """
        pass

    def force(self, u):
        """Return the thrust force exerted by the farm for given velocity or speed u"""
        return self._turbine_density * self.turbine_specification.force(u)

    def power(self, u):
        """Return the amount of power produced by the farm for given speed u"""
        return self._turbine_density * self.turbine_specification.power(u)

    def power_integral(self, u):
        """Return the UFL expression for the total amount of power produced by the farm for given speed u"""
        return self.power(u, turbine_density=self.turbine_density) * self.site_dx
