import numpy as np
import sys
import pathlib
import pandas as pd

dirPkR = str(pathlib.Path(__file__).parent.absolute())+'/logPkRatio/'

def print_to_c_1d(array, arrName=None, max_line_width=80):
    s = np.array2string(array, formatter={'float_kind': lambda x: "%.6e" % x}, separator=',', max_line_width=max_line_width)
    N = len(array)
    if arrName is not None:
        print('static double '+arrName+f'[{N}] = ', end='')
        print('{'+s[1:-1]+'};', end='\n\n')
    else:
        print('{'+s[1:-1]+'}', end='')


def print_to_c_2d(array, arrName=None, max_line_width=80):
    Nrows, Ncols = array.shape
    if arrName is not None:
        print('static double '+arrName+f'[{Nrows}][{Ncols}] = ', end='')
    print('{', end='')
    for i in range(Nrows):
        print_to_c_1d(array[i, :], max_line_width=80)
    print('};', end='\n\n')


if __name__ == '__main__':
    ''' Terminal calling example:
        $ python to_c_array.py TNG100 > baryons.h
        $ python to_c_array.py mb2 >> baryons.h

        All available sim_key:
            'TNG100', 'mb2', 'eagle', 'illustris', 'HzAGN'
            'cowls_AGN_T80', 'cowls_AGN_T85', 'cowls_AGN_T87'
            'BAHAMAS_T78', 'BAHAMAS_T76', 'BAHAMAS_T80'
            'owls_AGN', 'owls_DBLIMFV1618', 'owls_NOSN_NOZCOOL', 'owls_NOSN', 'owls_NOZCOOL', 'owls_REF', 'owls_WDENS', 'owls_WML1V848', 'owls_WML4'
    '''
    sim_key = str(sys.argv[1])
    filename = f'logPkRatio_{sim_key}.dat'

    TbArr = pd.read_csv(dirPkR+filename, sep='\s+').to_numpy()
    logk = TbArr[:, 0]      # the 0th col is logk
    logPkR = TbArr[:, 1:]   # the 1th~end-th cols are body of logPkRatio

    print_to_c_1d(logk, arrName=f'logkBins_{sim_key}')
    print_to_c_2d(logPkR, arrName=f'logPkR_{sim_key}')
