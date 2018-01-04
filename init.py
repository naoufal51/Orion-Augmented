# % ########################################################
# % #                  University of Kassel                #
# % # Department for RF-Techniques / Communication Systems #
# % #         Wilhelmshoeher Allee 73, 34121 Kassel        #
# % #              S. Semmelrodt / August 2003             #
# % # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# % #  Comments and bug reports to: semmelr@uni-kassel.de  #
# % ########################################################
# %
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

from numpy import *
from scipy.optimize import fminbound
from ic import *
from cost_func import *
"""
Options to use
"""
opts = {'dispmode': 0,
        'plotmode': 0,
        'diagmode': 'off',
        'ic_mode':'parallel',
        'maxfunceval':50,
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
def nextpow2(i):
    n = 1
    while n < i: n *= 2
    return n

   # This is internal function used by fft(), because the FFT routine
   # requires that the data size be a power of 2.

# function [beta, f1, f2,CostFunction1,CostFunction2] = sage_2d(y, p, varargin);

def sage_2d(y, p):
    # Get dimensions
    [N1, N2] = y.shape

     # Setup internal parameters
    ModeFlag     = 2
    ParamCount   = 3
    if p > 1:
       Theta = zeros(ParamCount,2)
    else:
       Theta = zeros(ParamCount,1)

    DispMode     = opts['dispmode']
    PlotMode     = opts['plotmode']
    IC_Mode      = opts['ic_mode']
    MaxCycleCnt  = abs(opts['maxcyclecnt'])
    DynamicRange = abs(opts['dynamicrange'])
    N1_FFT       = abs(opts['n1_fft'])

    if opts['n1_fft'] != 0:
       N2_FFT = abs(opts['n2_fft'])
    else:
       N2_FFT = 2**nextpow2(1024*N2/N1)

    if opts['tol_em'] != 0:
       Tol = abs(opts['tol_em'])
    else:
       Tol = [1/(2*N1_FFT), 1/(2*N1_FFT), 1/(2*N1_FFT)];


     # Initialise parameter matrix with initial guess
    Theta     = zeros(ParamCount, p)

    # SAGE Initialisation using FFT evaluation
    for Comp_Index in range(1,p):

        # Expectation Step: cancel all interferers
        x_i = ic(y, Theta, Comp_Index, 'static_2d', 'serial')

        # Calculate initialisation cost function

        CostFunction1 = fft.fftshift(sum(abs(fft(x_i,N1_FFT,1))**2,2));

        # Estimate initial parameter (Maximization step)

        # [Dummy, Index] = max(CostFunction1);

        Dummy, Index = CostFunction1.max(0), CostFunction1.argmax(0)

        # Theta(2,Comp_Index) = (1/N1_FFT * (Index-1)) - 0.5

        Theta[1, Comp_Index] = (1 / N1_FFT * (Index - 1)) - 0.5

        # Calculate initialisation cost function

        CostFunction2 = fft.fftshift(abs(fft.fft((exp(-2j * pi * (arange(0, N1).transpose() * ones([1, N2]) * 0.5))), N2_FFT, 1).sum(axis=0)) ** 2)

        # Estimate initial parameter (Maximization step)

        Dummy, Index = CostFunction2.max(0), CostFunction2.argmax(0)
        Theta[2,Comp_Index] = (1/N2_FFT * (Index-1)) - 0.5

        # Initial guess of the path weight
        Theta[0, Comp_Index] = 1/(N1*N2) * cost_func([], 0, Theta[:, Comp_Index], x_i, 'static_2d', 0)

        # Check component count
        if abs(Theta[1, Comp_Index]) < 10^(-DynamicRange/20) * max(abs(Theta[1, 1])):
            Theta[:,Comp_Index] = zeros[ParamCount,1]
            print('Component %d Out of dynamic range' %Comp_Index)
            return

    return Theta

