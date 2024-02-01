
"""Contains class Setup."""

# import importlib.util
# import importlib_resources
import sys
import io
import pandas as pd

class Setup(object):
    """Class defines a setup with population parameters.

    Methods:
    setParameters -- sets values for all parameters in the simulation
    save -- saves Setup parameters to given file location as a CSV file
    load -- loads Setup parameters from CSV file at given location
    """

    def __init__(self):
        """Create a new Setup."""

        super(Setup, self).__init__()

        self.parameter_names = [
            'id',
            'num_loci','possible_alleles',
            'fitnessHost','contactHost','receiveContactHost','mortalityHost',
            'natalityHost','recoveryHost','migrationHost',
            'populationContactHost','receivePopulationContactHost',
            'mutationHost','recombinationHost',
            'fitnessVector','contactVector','receiveContactVector',
            'mortalityVector','natalityVector','recoveryVector',
            'migrationVector',
            'populationContactVector','receivePopulationContactVector',
            'mutationVector','recombinationVector',
            'contact_rate_host_vector',
            'transmission_efficiency_host_vector',
            'transmission_efficiency_vector_host',
            'contact_rate_host_host',
            'transmission_efficiency_host_host',
            'mean_inoculum_host','mean_inoculum_vector',
            'variance_inoculum_host','variance_inoculum_vector',
            'recovery_rate_host','recovery_rate_vector',
            'mortality_rate_host','mortality_rate_vector',
            'recombine_in_host','recombine_in_vector',
            'num_crossover_host','num_crossover_vector',
            'mutate_in_host','mutate_in_vector',
            'death_rate_host','death_rate_vector',
            'birth_rate_host','birth_rate_vector',
            'vertical_transmission_host','vertical_transmission_vector',
            'inherit_protection_host','inherit_protection_vector',
            'protection_upon_recovery_host','protection_upon_recovery_vector'
            ]

    def setParameters(self,**kwargs):
        """Initializes a new Setup.

        Arguments:
        loaded_from_file -- whether parameters already imported (Boolean)

        Keyword arguments:
        **kwargs -- setup parameters and values
        """
        for parameter, value in kwargs.items():
            setattr( self, parameter, value )

        self.num_loci = int(self.num_loci)

        if isinstance(self.possible_alleles, list):
            self.possible_alleles = self.possible_alleles
        else:
            self.possible_alleles = [self.possible_alleles] * self.num_loci
                # possible_alleles must be a list with all available alleles for
                # each position

    def save(self,save_to_file):
        """
        Saves Setup parameters to given file location as a CSV file.

        Functions (e.g. fitness functions) cannot be saved in this format.

        Arguments:
        save_to_file -- file path and name to save  parameters under (String)
        """
        out = 'Parameter,Value\n'
        function_counter = 0
        for parameter in self.parameter_names:
            if hasattr( getattr(self,parameter), '__call__'):
                    # checks if parameter is function
                out = out + parameter + ',#FUNCTION:'+parameter+'Function\n'
                function_counter += 1
            else:
                out = out + parameter + ',' + str( getattr(self,parameter) )+'\n'

        file = open(save_to_file,'w')
        file.write(out)
        file.close()

        print('Parameter file saved.')

    def load(self,file,preset=None,**kwargs):
        """
        Loads Setup parameters from CSV file at given location.

        Functions (e.g. fitness functions) are loaded as function names preceded
        by the keyword: "#FUNCTION:" . The functions themselves are defined in a
        separate file.

        Arguments:
        file -- file path to CSV file with parameters (String)

        Keyword arguments:
        preset -- if using preset parameters, 'host-host' or 'vector-borne'
            (String, default None)
        **kwargs -- setup parameters and values
        """
        if preset is None:
            df = pd.read_csv(file)
        else:
            file_name = preset+'.csv'
            parameter_dir = importlib_resources.files('opqua') / 'parameters'
            parameter_bytes = (parameter_dir / file_name).read_bytes()
            parameter_str = str(parameter_bytes,'utf-8')
            parameter_data = io.StringIO(parameter_str)
            df = pd.read_csv(parameter_data)

        # if function_file_path is not None:
        #     spec = importlib.util.spec_from_file_location(
        #         'function_params', function_file_path
        #         )
        #     function_params = importlib.util.module_from_spec(spec)
        #     sys.modules['function_params'] = function_params
        #     spec.loader.exec_module(function_params)

        for i,row in df.iterrows():
            if '#FUNCTION:' in str(row['Value']): # checks if parameter is function
                function_name = row['Value'][len('#FUNCTION'):].strip()
                if ( #function_file_path is None or
                        not hasattr(function_params, function_name) ):
                    setattr( self, row['Parameter'], lambda g:1 )
                else:
                    setattr( self, row['Parameter'], getattr(
                            function_params, function_name
                            ) )
            elif not isinstance(row['Value'], str) and pd.isna(row['Value']):
                setattr( self, row['Parameter'], None )
            elif isinstance(row['Value'], str) and (
                    row['Value'].replace('.','',1).isdigit()
                        # checks if number, https://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-represents-a-number-float-or-int
                    or 'e-' in row['Value'] or 'e+' in row['Value']
                        # or if scientific notation
                    ):
                setattr( self, row['Parameter'], float(row['Value']) )
            else:
                setattr( self, row['Parameter'], row['Value'] )

        self.setParameters(**kwargs)

        print('Parameter file loaded.')
