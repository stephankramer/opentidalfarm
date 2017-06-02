from .base_turbine import BaseTurbine
from .controls import Controls
from dolfin import pi, dot
from ..helpers import tabulated_expression
import numpy
import math

class BumpTurbine(BaseTurbine):
    """Turbine modelled as a bump function of increased bottom friction.

    The amount of friction is worked out from the thrust_coefficient 
    and diameter in such a way that the force exerted is given by:

    .. math:: F = \\tfrac 12 C_t A_t u^2,
    
    where :math:`C_t` is the `thrust_coefficient`, :math:`A_t = \pi D^2/4` 
    the wetted cross section for `diameter` :math:`D`, :math:`u` the
    speed. Instead of a thrust_coefficient a thrust_curve can be
    provided as a list of speed, thrust pairs `(u, C_t)`.
    
    If depth is provided it is used to apply the correction from
    [1] S.C. Kramer, M.D. Piggott, Renewable Energy, Vol. 92, http://doi.org/10.1016/j.renene.2016.02.022,
    to estimate the free stream velocity from the local 
    depth-averaged velocity in the force and power calculations."""

    def __init__(self, thrust_coefficient=None,
                 thrust_curve=None,
                 diameter=20., minimum_distance=None,
                 depth=None,
                 controls=Controls(position=True)):
        # Check for a given minimum distance.
        if minimum_distance is None:
            minimum_distance = diameter*1.5

        self.depth = depth

        # Initialize the base class.
        super(BumpTurbine, self).__init__(
                thrust_coefficient=thrust_coefficient,
                thrust_curve=thrust_curve,
                diameter=diameter,
                controls=controls)
        if depth is not None:
            self._correct_thrust_curve(depth)

    def _correct_thrust_curve(self, depth):
        """Set thrust curve given as a list of pairs `(u, Ct)`."""
        if self._thrust_curve is None:
            return
        if not hasattr(self._original_thrust_curve):
            self._original_thrust_curve = self._thrust_curve
        # translate thrust curve from (u_inf, C_t) to (u_bar, C_t)
        # where u_inf and u_bar are free-stream and depth-averaged speeds resp.
        self._thrust_curve = []
        for u_inf, Ct in self._original_thrust_curve:
            # eqn (13) from [1]
            self.area_ratio = pi * self.diameter / (4*self.depth)
            u_bar = (1+(1-self.area_ratio*Ct)**0.5)/2 * u_inf
            self._thrust_curve.append((u_bar, Ct))

    def force(self, u):
        """Return the thrust force exerted by a turbine for given velocity u

        :param u: velocity vector or speed
        :type u: dolfin.Function or float."""
        if self.depth is None:
            correction = 1
        else:
            # correction from eqn (15) in [1]
            u_norm = dot(u,u) ** 0.5
            C_t = self.get_thrust_coefficient(u_norm)
            # ratio = A_t/\hat{A_t} = pi (D/2)^2 / D*H = pi * D / (4 * H)
            # note that here we always use linear depth to avoid adding unnec. non-linearities
            area_ratio = pi * self.diameter / (4*self.depth)
            correction = 4/(1+(1-area_ratio*C_t)**0.5)**2
        return correction * super(BumpTurbine, self).force(self, u)


    def power(self, u):
        """Return the amount of power produced by a turbine for given speed u

        :param u: speed (scalar)
        :type u: dolfin.Function or float."""
        if self.depth is None:
            correction = 1
        else:
            # eqn (C.2) from [1]
            # [1] S.C. Kramer, M.D. Piggott http://doi.org/10.1016/j.renene.2016.02.022
            # where turbine_field = C_t A_t / (2*Integral(bump))
            C_t = self.get_thrust_coefficient(u)
            # ratio = A_t/\hat{A_t} = pi (D/2)^2 / D*H = pi * D / (4 * H)
            area_ratio = pi * self.diameter / (4*self.depth)
            correction = 4 * (1+(1-C_t)**0.5)/(1+(1-area_ratio*C_t)**0.5)**3

        return correction * super(BumpTurbine, self).power(self, u)


def thrust_from_power_coefficient(Cp):
    # we want to solve Cp = Ct (1+sqrt(1-Ct))/2 for Ct
    # by substitution we have: Cp^2-Ct*Cp = -Ct^3/4
    # i.o.w.: Ct^3 - 4*Cp*Ct + 4 Cp^2 = 0
    r = numpy.roots([1., 0., -4*Cp, 4*Cp**2])

    # check which roots are actually solutions to the original equation
    Ct = []
    for root in r:
        if root>=0. and root<=1. and abs(root*(1+math.sqrt(1-root))/2.-Cp)<1e-12:
            Ct.append(root)

    # for Cp near the Betz limit, we get 2 solutions, so we take the smallest
    return min(Ct)


def standard_thrust_curve(cut_in_speed=1.0, rated_speed=3.0, cut_out_speed=6.0,
        thrust_coefficient_below_rated_speed=0.8, power_coefficient_below_rated_speed=None,
        delta_u=0.1):
    # for u<u0 and u>un we assume Ct = 0 - so we automatically get a cut_in and out
    # to control the interval over which the thrust kicks in around the cut_in_speed
    # we explcitly set the high speed at which Ct is still 0
    curve = [[cut_in_speed-delta_u, 0.0],
             [cut_in_speed, thrust_coefficient_below_rated_speed],
             [rated_speed, thrust_coefficient_below_rated_speed]]

    if power_coefficient_below_rated_speed is None:
        Ct = thrust_coefficient_below_rated_speed
        power_coefficient_below_rated_speed = Ct * (1.+math.sqrt(1.-Ct))/2.

    for u in numpy.arange(rated_speed+delta_u, cut_out_speed+1e-6, delta_u):
        Cp = power_coefficient_below_rated_speed * rated_speed**3/u**3
        Ct = thrust_from_power_coefficient(Cp)
        curve.append([u, Ct])

    curve.append([cut_out_speed+delta_u, 0])
    return curve
