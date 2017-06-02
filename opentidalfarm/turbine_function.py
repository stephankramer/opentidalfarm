import copy
import numpy
from dolfin import *
from dolfin_adjoint import *

__all__ = ["TurbineFunction"]

class TurbineFunction(object):

    def __init__(self, cache, V, turbine_specification):

        self._parameters = copy.deepcopy(cache._parameters)
        self._turbine_specification = turbine_specification
        self._cache = cache

        # Precompute some turbine parameters for efficiency.
        self.x = interpolate(Expression("x[0]", degree=1), V).vector().array()
        self.y = interpolate(Expression("x[1]", degree=1), V).vector().array()
        self.V = V


    def __call__(self, name="", derivative_index=None, derivative_var=None,
                 timestep=None):
        """If the derivative selector is i >= 0, the Expression will compute the
        derivative of the turbine with index i with respect to either the x or y
        coorinate."""

        params = self._parameters

        if derivative_index is None:
            position = params["position"]
        else:
            position = [params["position"][derivative_index]]

        ff = numpy.zeros(len(self.x))
        # Ignore division by zero.
        numpy.seterr(divide="ignore")
        eps = 1e-12


        for x_pos, y_pos in position:
            radius = self._turbine_specification.radius
            x_unit = numpy.minimum(
                numpy.maximum((self.x-x_pos)/radius, -1+eps), 1-eps)
            y_unit = numpy.minimum(
                numpy.maximum((self.y-y_pos)/radius, -1+eps), 1-eps)

            # Apply chain rule to get the derivative with respect to the turbine
            # friction.
            exp = numpy.exp(-1./(1-x_unit**2)-1./(1-y_unit**2)+2)

            if derivative_index is None:
                ff += exp

            if derivative_var == "turbine_pos_x":
                ff += exp*(-2*x_unit/((1.0-x_unit**2)**2))*(-1.0/radius)

            elif derivative_var == "turbine_pos_y":
                ff += exp*(-2*y_unit/((1.0-y_unit**2)**2))*(-1.0/radius)

        # Reset numpy to warn for zero division errors.
        numpy.seterr(divide="warn")

        f = Function(self.V, name=name, annotate=False)
        f.vector().set_local(ff)
        f.vector().apply("insert")
        return f
