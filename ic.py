# % ########################################################
# % #                  University of Kassel                #
# % # Department for RF-Techniques / Communication Systems #
# % #         Wilhelmshoeher Allee 73, 34121 Kassel        #
# % #              S. Semmelrodt / August 2003             #
# % # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# % #  Comments and bug reports to: semmelr@uni-kassel.de  #
# % ########################################################
# %
# % Function:   Interference cancellation algorithm for the
# %             separation of superimposed signals.
# %
# % Call:       x_i = ic(y, Theta, i, ModelType, IC_Mode)
# %
# % Parameters: - y         : Vector with signal samples
# %             - Theta     : Parameter matrix
# %             - i         : Index to component of interest
# %             - ModelType : Model type string ('static_1d',
# %                           'static_2d', 'linear_1d')
# %             - IC_Mode   : 'serial' or 'parallel' mode
# %
# % Output:     - x_i       : Signal component
# %
from numpy import *
from string import *
from model_func import *


def ic(y, Theta, i, IC_Mode):
    # Get Number of components
    # CompCount = Theta.size(1)
    CompCount = array([Theta]).T.shape[1]

    #  Select processing mode
    if str.lower(IC_Mode)=='parallel':
        # Cancel all interferers
        CompVector = arange(1,CompCount+1)
        # Apart from the considered component
        CompVector = delete(CompVector, where(CompVector == i))

    elif str.lower(IC_Mode)=='serial':
        # Cancel dominant interferers
        CompVector = arange(1,i)

    else:
        print('Unknown IC-Mode !')

    # Setup data
    x_i = y

    # Cancel components sequentially
    for Index in CompVector:
        # Cancel current component
        x_i = x_i - model_func(Theta[:,Index], y.shape)
    return x_i


# [[2.0245 + 0.0245i   1.8826 + 0.5985i   1.6330 + 1.1905i],
#    [1.9182 + 0.6367i   1.6223 + 1.1807i   1.1667 + 1.6083i],
#    [1.6463 + 1.2073i   1.2108 + 1.6498i   0.6451 + 1.9191i]]
#
# array([[2.0245 + 0.0245j, 1.8826 + 0.5985j, 1.6330 + 1.1905j],
#        [1.9182 + 0.6367j, 1.6223 + 1.1807j, 1.1667 + 1.6083j],
#        [1.6463 + 1.2073j, 1.2108 + 1.6498j, 0.6451 + 1.9191j ]])