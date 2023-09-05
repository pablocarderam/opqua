
"""Contains class Intervention."""

class Intervention(object):
    """Class defines a new intervention to be done at a specified time.

    Attributes:
        time (number): time at which intervention will take place.
        method_name (String): intervention to be carried out, must correspond to the
            name of a method of the Model object.
        args (array-like): contains arguments for function in positinal order.
        model (Model object): Model object this intervention is associated to.
    """

    def __init__(self, time, method_name, args, model):
        """Create a new Intervention.

        Arguments:
            time (number): time at which intervention will take place.
            method_name (String): intervention to be carried out, must correspond to the
                name of a method of the Model object.
            args (array-like): contains arguments for function in positinal order.
            model (Model object): Model object this intervention is associated to.
        """
        super(Intervention, self).__init__()
        self.time = time
        self.intervention = method_name
        self.args = args
        self.model = model

    def doIntervention(self):
        """Execute intervention function with specified arguments."""

        getattr(self.model, self.intervention)(*self.args)
