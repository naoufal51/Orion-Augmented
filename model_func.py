# % ########################################################
# % #                  University of Kassel                #
# % # Department for RF-Techniques / Communication Systems #
# % #         Wilhelmshoeher Allee 73, 34121 Kassel        #
# % #              S. Semmelrodt / August 2003             #
# % # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# % #  Comments and bug reports to: semmelr@uni-kassel.de  #
# % ########################################################
# %
# % Function:   Implements several signal models used in the context
# %             with optimization procedures like the EM, SAGE or
# %             RELAX based algorithms, respectively.
# %
# % Call:       y = model_func(Theta, Dim, ModelType)
# %
# % Parameters: - Theta      : Parameter matrix
# %             - Dim        : Signal dimension
# %             - ModelType  : model type string ('static_1d',
# %                            'static_2d', 'linear_1d')
# %
# % Output:     - y          : Synthesized signal
# %
from numpy import *

def model_func(Theta, Dim):

    # Initialise variables
    beta = Theta[0]
    f1 = Theta[1]
    f2 = Theta[2]

    N1 = Dim[0]
    N2 = Dim[1]
    # Setup domains

    if N1 % 2 == 0 and N1 % 2 == 0:
        n1 = array([arange(0, N1)]).transpose()
        n2 = array([arange(0, N2)])
    elif N1 % 2 == 0 and N1 % 2 != 0:
        n1 = array([arange(0, N1)]).transpose()
        n2 = array([arange(-floor(N2 / 2), floor(N2 / 2)+1)])
    elif N1 % 2 != 0 and N1 % 2 == 0:
        n1 = array([arange(-floor(N1 / 2), floor(N1 / 2)+1)]).transpose()
        n2 = array([arange(0, N2)])
    else:
        n1 = array([arange(-floor(N1 / 2), floor(N1 / 2)+1)]).transpose()
        n2 = array([arange(-floor(N2 / 2), floor(N2 / 2)+1)])

    # Calculate signal model for given parameters
    y = beta * exp(2j*pi * f1 * n1) * exp(2j*pi * f2 * n2)

    return y