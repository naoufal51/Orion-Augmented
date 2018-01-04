"""
# % Function:   Estimates weights and frequency parameters of a
# %             multi-component harmonic two-dimensional signal
# %             using the SAGE algorithm.
# %
# % Call:       [beta, f1, f2] = sage_2d(y, p, <fs>, <opts>)
# %
# % Parameters: - y    : Matrix with signal samples
# %             - p    : Model order (number of significant components)
# %             - fs   : Sampling frequency in Hz in both dimensions (default: fs = 1 Hz)
# %             - opts : Option structure with parameter settings
# %                      (see function em_opts)
# %
# % Output:     - beta : Estimated weights of the signal components
# %             - f1   : Estimated frequencies (1st dimension)
# %             - f2   : Estimated frequencies (2nd dimension)
# %
# % Example:    y = sg_cissoid_2d([1 1], [0.2 0.3] , [0.1 0.1], [21 5], [5 2], 30);
# %             options = em_opts
# %             [beta, f1, f2] = sage_2d(y, 2, [5 2], options)
# %
"""
from numpy import *
from scipy.optimize import fminbound
from ic import *
from cost_func import *
from nextpow2 import *


def sage(y, p):
    # Get dimensions

    [N1, N2] = y.shape

    # Setup internal parameters
    ModeFlag = 2
    ParamCount = 3
    if p > 1:
        Theta = zeros(ParamCount, 2)
    else:
        Theta = zeros(ParamCount, 1)

    DispMode = opts['dispmode']
    PlotMode = opts['plotmode']
    IC_Mode = opts['ic_mode']
    MaxCycleCnt = abs(opts['maxcyclecnt'])
    DynamicRange = abs(opts['dynamicrange'])
    N1_FFT = abs(opts['n1_fft'])

    if opts['n1_fft'] != 0:
        N2_FFT = abs(opts['n2_fft'])
    else:
        N2_FFT = 2 ** nextpow2(1024 * N2 / N1)

    if opts['tol_em'] != 0:
        Tol = abs(opts['tol_em'])
    else:
        Tol = [1 / (2 * N1_FFT), 1 / (2 * N1_FFT), 1 / (2 * N1_FFT)]

    # Initialise parameter matrix with initial guess
    Theta_Old = zeros([ParamCount, p])
    Theta = zeros([ParamCount, p])
    Theta_Threshold = Tol * ones([1, p])
    CycleCtr = 1

    # SAGE Initialisation using FFT evaluation
    for Comp_Index in range(0, p):

        # Expectation Step: cancel all interferers
        x_i = ic(y, Theta, Comp_Index, 'serial')

        # Calculate initialisation cost function

        CostFunction1 = fft.fftshift((abs(fft.fft(x_i, N1_FFT, 1)) ** 2).sum(axis=1))

        # Estimate initial parameter (Maximization step)

        Dummy, Index = CostFunction1.max(0), CostFunction1.argmax(0)

        Theta[1, Comp_Index] = (1 / N1_FFT * (Index - 1)) - 0.5

        # Calculate initialisation cost function

        CostFunction2 = fft.fftshift(abs(fft.fft((exp(-2j * pi * (array([arange(0, N1)]).transpose() * ones([1, N2]) * Theta[1, Comp_Index])) * x_i),N2_FFT, 1).sum(axis=0)) ** 2)

        # Estimate initial parameter (Maximization step)

        Dummy, Index = CostFunction2.max(0), CostFunction2.argmax(0)

        Theta[2, Comp_Index] = (1 / N2_FFT * (Index - 1)) - 0.5

        # Initial guess of the path weight
        Theta[0, Comp_Index] = 1 / (N1 * N2) * cost_func([], 0, Theta[:, Comp_Index], x_i, 0)

        # Check component count
        if abs(Theta[0, Comp_Index]) < 10 ^ (-DynamicRange / 20) * (abs(Theta[0, 0])):
            Theta[:, Comp_Index] = zeros[ParamCount, 1]
            print('Component %d Out of dynamic range' % Comp_Index)
            sys.exit()
        # Process until parameter converges
        CycleCtr = 0

    while 1:

        # Display EM-Cycle

        print("EM-Cycle: %d" % CycleCtr)

        # Process EM step for each component
        for Comp_Index in range(0, p):
            # Expectation Step: cancel all interferers
            x_i = ic(y, Theta, Comp_Index, IC_Mode)

            # Coordinate wise updating of the parameter vector
            Theta[1, Comp_Index] = fminbound(cost_func, Theta[1, Comp_Index] - 1 / (CycleCtr * N1),
                                             Theta[1, Comp_Index] + 1 / (CycleCtr * N1),
                                             (1, Theta[:, Comp_Index], x_i, ModeFlag),
                                             opts.tol_ls, opts.maxitercnt)
            # Coordinate wise updating of the parameter vector
            Theta[2, Comp_Index] = fminbound(cost_func, Theta[2, Comp_Index] - 1 / (CycleCtr * N2),
                                             Theta[2, Comp_Index] + 1 / (CycleCtr * N2),
                                             (1, Theta[:, Comp_Index], x_i, ModeFlag),
                                             opts.tol_ls, opts.maxitercnt)
            # Updating of the path weight
            Theta[1, Comp_Index] = 1 / (N1 * N2) * cost_func([], 0, Theta[:, Comp_Index], x_i, 0)

            # Check component count
            if abs(Theta(0, Comp_Index)) < 10 ^ (-DynamicRange / 20) * (abs(Theta(0, 0))):
                Theta[:, Comp_Index] = zeros(ParamCount, 1)
                # ddisp(DispMode, ['Component[',num2str(Comp_Index),']: Out of dynamic range'])
            else:
                # Display results for current iteration

                #    ddisp(DispMode, ['Component[',num2str(Comp_Index),']: ',num2str(abs(Theta(1,Comp_Index))),' ', num2str(Theta[2:end,Comp_Index].transpose)])
                print('ok')

                # Check for convergence (Parameter of two iterations are nearly constant)
            if all(abs(Theta - Theta_Old) < Theta_Threshold):
                #    ddisp(DispMode, ['Algorithm converged after 1+', num2str(CycleCtr),' EM-Cycles']);
                #    ddisp(DispMode, ['Mean convergence tolerance = ', num2str(sum(abs(Theta - Theta_Old)./p,2)')]);
                sys.exit()

            else:
                # Store current parameter vector
                Theta_Old = Theta

            # Increment EM cycle counter
            CycleCtr += 1
            ModeFlag = 2

            # Check cycle counter
            if CycleCtr > MaxCycleCnt:
                # ddisp(DispMode, ['Maximum number of EM-Cycles exceeded -> No convergence']);
                return

                # Return nonzeros components only
    Index = find(Theta[0, :])
    Theta = Theta[:, Index]

    # Return parameters in common form
    beta = Theta[0, :].transpose()
    f1 = (Theta[1, :]).transpose()
    f2 = (Theta[2, :]).transpose()

    return beta, f1, f2, CostFunction1, CostFunction2


"""
Options to use
"""

opts = {'dispmode': 0,
        'plotmode': 0,
        'diagmode': 'off',
        'ic_mode': 'parallel',
        'maxfunceval': 50,
        'maxitercnt': 20,
        'maxcyclecnt': 250,
        'n1_fft': 1024,
        'n2_fft': 0,
        'tol_em': 0,
        'tol_ls': 1e-6,
        'dynamicrange': 30}

"""
Find 2^n that is equal to or greater than.
"""

# This is internal function used by fft(), because the FFT routine
# requires that the data size be a power of 2.

# function [beta, f1, f2,CostFunction1,CostFunction2] = sage_2d(y, p, varargin);
