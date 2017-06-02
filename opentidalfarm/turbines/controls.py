class Controls(object):
    """Specifies the controls for optimisation.

    This class holds the controls for optimisation, such as the position and the
    density of the turbines. The user initializes this class with their desired
    control parameters.
    """
    def __init__(self, position=False, density=False, dynamic_thrust=False):
        """Initialize with the desired controls.

        :param bool position: Whether or not turbine position is a control.
        :param bool density: Whether or not turbine density is a control.
        :param bool dynamic_thrust: Whether or not dynamic thrust is a control.
        """

        self._controls = {"position": False,
                          "density": False,
                          "dynamic_thrust": False}

        def _process(key, value):
            """Check value is of type bool. Raise ValueError if it is not."""
            try:
                assert isinstance(value, bool)
                # Change the control value in the dictionary.
                self._controls[key] = value
            # Raise an error if a boolean was not given.
            except AssertionError:
                raise ValueError("%s must be a boolean (%s)." %
                                 (key.capitalize(), str(type(value))))

        # Process the given values
        _process("position", position)
        _process("density", density)
        _process("dynamic thrust", dynamic_thrust)


    def __str__(self):
        """Returns a string representation of the enabled control parameters."""
        string = "Control parameters:"
        # Get enabled controls.
        enabled = [key for key in self._controls if self._controls[key]]
        # Add the keys to the string to be returned.
        if len(enabled) > 0:
            for control in enabled:
                string += "\n - %s" % control.capitalize()
        else:
            string += " no control parameters have been enabled!"
        return string


    @property
    def position(self):
        """Whether position is enabled as a control parameter.

        :returns: True if position is enabled as a control parameter.
        :rtype: bool
        """
        return self._controls["position"]


    @property
    def density(self):
        """Whether density is enabled as a control parameter.

        :returns: True if density is enabled as a control parameter.
        :rtype: bool
        """
        return self._controls["density"]


    @property
    def dynamic_thrust(self):
        """Whether dynamic thrust is enabled as a control parameter.

        :returns: True if dynamic thrust is enabled as a control parameter.
        :rtype: bool
        """
        return self._controls["dynamic thrust"]
