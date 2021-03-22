import numpy as np

def model_DMO(pars):
    '''generate model vector given cosmologicl parameters from cocoa
        Args:
            pars: input parameters to run the the theoritical DMO model vector
            e.g. pars = {'Omega_m': 0.3, 
                         'n_s': 0.97, 
                         'A_s': 2.19e-9}
        Returns:
            modelv: 1D array
                (scale-cutted) model vector
    '''

    # Fill cocoa interaction code bit here.
    
    return modelv

def model_bary(pars, scenario):
    '''generate baryon contaminated model vector given the input parameter and the baryonic scenario
        Args:
            pars: parameter dict.
            scenario: string
                available baryon scenario keywords are: 
                'TNG100', 'mb2', 'eagle', 'illustris', 'HzAGN', 'cowls_AGN_T80', 'cowls_AGN_T85', 'cowls_AGN_T87', 'BAHAMAS_T78', 'BAHAMAS_T76', 'BAHAMAS_T80'
    '''
    # Fill cocoa interaction code bit here.
    return modelv  

class PCA():
    '''generating PCs from given baryonic scearnios'''
    def __init__(self, pars, simKeys, COV):
        '''
            Args:
                par: parameter dictionary
                    the base parameters to generate the baryonic model vectors to constructs PCs.
                    e.g. pars = {'Omega_m': 0.3, 
                                'n_s': 0.97, 
                                'A_s': 2.19e-9}

                simKeys: 
                    a list of baryoinc scenarios
                    e.g. simKeys = ['TNG100', 'mb2', 'eagle', 'illustris', 'HzAGN',
                                    'cowls_AGN_T80', 'cowls_AGN_T85', 'cowls_AGN_T87',
                                    'BAHAMAS_T78', 'BAHAMAS_T76', 'BAHAMAS_T80']
                
                COV: ndarray, NxN
                    covariance subjected to the (scale-cutted) data vector
        '''

        self.pars = pars
        self.simKeys = simKeys
        self.Nscenario = len(simKeys)
        self.Ndata = len(COV)

        self._init_COV(COV)
        

    def _init_COV(self, COV):
        '''init quantities related to covariance matrix'''
        self.COV = COV
        self.invCOV = np.linalg.inv(COV)   # inverse COV

        # Cholesky decomposition on COV
        self.L = np.linalg.cholesky(COV)
        self.invL = np.linalg.inv(self.L)

    def build_Ratio(self, pars=None):
        '''Build the Ratio matrix (at the given pars)
            Returns:
                Ratio: ndarray (self.Ndata x self.Nscenario)
        '''

        if pars is None:
            pars = self.pars
        
        Ratio = np.zeros((self.Ndata, self.Nscenario))
        modelv_dmo = model_DMO(pars)

        #for j in range()


