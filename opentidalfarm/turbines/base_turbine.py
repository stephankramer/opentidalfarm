import dolfin
import math

class BaseTurbine(object):
    """A base turbine class from which others are derived."""
    def __init__(self, thrust_coefficient=None,
                 thrust_curve=None,
                 diameter=20., 
                 controls=None):
        # Possible turbine parameters.
        self._diameter = diameter
        self._controls = controls
        self._swept_area = math.pi * self._diameter**2 / 4.0

        if thrust_curve is None:
            if thrust_coefficient is None:
                self._thrust_coefficient = 0.8
            else:
                self._thrust_coefficient = thrust_coefficient
            self._thrust_curve = None
        elif thrust_coefficient is None:
            self._thrust_curve = thrust_curve
            self._thrust_coefficient = None
        else:
            raise ValueError("Either thrust_coefficient or thrust_curve should be provided, not both!")


    @property
    def diameter(self):
        """The diameter of a turbine.
        :returns: The diameter of a turbine.
        :rtype: float
        """
        if self._diameter is None:
            raise ValueError("Diameter has not been set!")
        return self._diameter


    @property
    def radius(self):
        """The radius of a turbine.
        :returns: The radius of a turbine.
        :rtype: float
        """
        return self.diameter*0.5

    def get_thrust_coefficient(self, u):
        """Get thrust coefficient C_t for given speed `u`""" 
        if self._thrust_curve is None:
            return self._thrust_coefficient
        else:
            return tabulated_expression(u, self._thrust_curve)

    def get_power_coefficient(self, u):
        """Get power coefficient C_P for given speed `u`""" 
        C_t = self.get_thrust_coefficient(u)
        return 0.5 * C_t * (1 + sqrt(1-C_t))

    def _set_controls(self, controls):
        self._controls = controls

    def _get_controls(self):
        if self._controls is not None:
            return self._controls
        else:
            raise ValueError("The controls have not been set.")

    controls = property(_get_controls, _set_controls, "The turbine controls.")

    def force(self, u):
        """Return the thrust force exerted by the turbines for given velocity or speed u
        
        :param u: velocity vector or speed
        :type u: dolfin.Function or float"""
        unorm = dolfin.dot(u,u)**0.5
        return 0.5 * self.get_thrust_coefficient(u_norm) * self._swept_area * u_norm * u


    def power(self, u):
        """Return the amount of power produced by the turbines for given speed u

        :param u: speed (scalar)
        :type u: dolfin.Function or float"""
        return 0.5 * self.get_power_coefficient(u) * self._swept_area * u**3
