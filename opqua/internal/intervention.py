
"""Contains class Intervention."""

class Intervention(object):
    """Class defines a new intervention to be done at a specified time.

    Methods:
    doIntervention -- executes intervention function with specified arguments
    """

    def __init__(self, time, method_name, args, model):
        """Create a new Intervention.

        Arguments:
        time -- time at which intervention will take place (number)
        method_name -- intervention to be carried out, must correspond to the
            name of a method of the Model object (String)
        args -- contains arguments for function in positinal order (array-like)
        model -- Model object this intervention is associated to (Model)
        """
        super(Intervention, self).__init__()
        self.time = time
        self.intervention = method_name
        self.args = args
        self.model = model

    def doIntervention(self):
        """Execute intervention function with specified arguments."""

        getattr(self.model, self.intervention)(*self.args)
