
"""Contains class Intervention."""

class Intervention(object):
    """Class defines a new intervention to be done at a specified time.

    Methods:
    doIntervention -- executes intervention function with specified arguments
    """

    def __init__(self, time, function, args):
        """Create a new Intervention.

        Arguments:
        time -- time at which intervention will take place (number)
        function -- intervention to be carried out (method of class Model)
        args -- contains arguments for function in positinal order (array-like)
        """
        super(Intervention, self).__init__()
        self.time = time
        self.intervention = function
        self.args = args

    def doIntervention(self):
        """Execute intervention function with specified arguments."""

        self.intervention(*self.args)
