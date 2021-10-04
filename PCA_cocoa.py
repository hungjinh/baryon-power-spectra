import numpy as np

def model_dmo(pars):
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
            
            Properties:
                self.PCs: the dictionary that stores principal components
                    self.PCs[0]: the 1st PC
                    self.PCs[1]: the 2nd PC, ...
                    self.PCs[self.Nscenario-1]: the last PC
        '''

        self.pars = pars
        self.simKeys = simKeys
        self.Nscenario = len(simKeys)
        self.Ndata = len(COV)

        self._init_COV(COV)
        
        self.PCs = self.gen_PC_dict()

    def _init_COV(self, COV):
        '''init quantities related to covariance matrix'''
        self.COV = COV
        self.invCOV = np.linalg.inv(COV)   # inverse COV

        # Cholesky decomposition on COV
        self.L = np.linalg.cholesky(COV)
        self.invL = np.linalg.inv(self.L)

    def build_Ratio(self):
        '''Build the Ratio matrix
            self.Ratio: ndarray (self.Ndata x self.Nscenario)
        '''
        
        self.Ratio = np.zeros((self.Ndata, self.Nscenario))
        self.modelv_dmo = model_dmo(self.pars)

        for j, scenario in enumerate(self.simKeys):
            modelv_bary = model_bary(self.pars, scenario)
            self.Ratio.T[j] = modelv_bary/self.modelv_dmo
        
        self.Ratio -= 1.
    
    def build_Delta(self):
        '''Build the difference matrix
            self.Delta: ndarray (self.Ndata x self.Nscenario)
        '''
        DeltaT = self.Ratio.T*self.modelv_dmo-self.modelv_dmo
        self.Delta = DeltaT.T
    
    def build_wDelta(self):
        '''Build the covariance weighted difference matrix'''
        self.wDelta = np.dot(self.invL, self.Delta)

    def SVD(self, Matrix=None):
        '''compute PCA given the input matrix
            default Matrix to build PCs: self.wDelta

            self.U matrix stroes the PCs.
            PC1 = self.U.T[0]
            PC2 = self.U.T[1]
        '''

        if Matrix is None:
            Matrix = self.wDelta

        self.U, self.Sdig, VT = np.linalg.svd(Matrix, full_matrices=True)
    
    def gen_PC_dict(self):

        self.build_Ratio()
        self.build_Delta()
        self.build_wDelta()
        self.SVD(Matrix=self.wDelta)

        PCs = {}
        for j in range(self.Nscenario):
            PCs[j] = self.U.T[j]
        
        return PCs

if __name__ == '__main__':
    pars = {}
    pars['Omega_m'] = 0.3
    pars['A_s'] = 2.19e-9
    pars['n_s'] = 0.97

    #cov = np.load()
    pca = PCA(pars, simKeys=['mb2', 'eagle', 'TNG100', 'illustris', 'cowls_AGN_T87'], cov=cov)
    PC1 = pca.PCs[0]
    PC2 = pca.PCs[1]

    # etc...