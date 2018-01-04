# % Call:       c = cost_func(ParamVal, ParamIndex, Theta_i_i, x_i, ModelType, ModeFlag)
# %
# % Parameters: - ParamVal   : Certain parameter value for evaluating the cost function
# %             - ParamIndex : Index within the parameter vector Theta_i_i
# %             - Theta_i    : Parameter matrix
# %             - Theta_i    : Parameter matrix
# %             - x_i        : separated signal
# %             - ModelType  : model type string ('static_1d',
# %                            'static_2d', 'linear_1d')
# %             - ModeFlag   : 0 -> standard, 2 -> initialization

from numpy import *


# noinspection PyRedundantParentheses
def cost_func(ParamVal, ParamIndex, Theta_i, x_i,ModeFlag):
    if (ParamIndex != 0):
        Theta_i[ParamIndex] = ParamVal

    # Get dimensions
    [N1, N2] = size(x_i, 0), size(x_i, 1)

    # Initialise variables
    C = zeros([N1, N2])

    f1 = Theta_i[1]
    f2 = Theta_i[2]

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
    C = exp(-2j * pi * f1 * n1) * exp(-2j * pi * f2 * n2) * x_i

    # Calculate discrete integral
    if ModeFlag == 0:
        # Return standard value
        c = sum(sum(C))

    elif ModeFlag == 1:

        if ParamIndex == 2:
            # Return inverted squared absolute value for initialisation
            c = -sum(abs(sum(C)) ** 2)
        else:
            # Return inverted absolute value
            c = -abs(sum(sum(C))) ** 2
    else:
        # Return inverted absolute value
        c = -abs(sum(sum(C))) ** 2
    return c
