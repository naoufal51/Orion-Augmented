"""
sage2d.py
Apply SAGE algotithm to estimate parameters
S. Semmelrodt / 2003 
M Mahfoudi / 2017
"""
import numpy as np
from scipy.optimize import fminbound
import sys

"""
Options to use
"""

opts = {
'ic_mode': 'parallel',
'maxfunceval': 50,
'maxitercnt': 20,
'maxcyclecnt': 500,
'n1_fft': 1024,
'n2_fft': 1024,
'tol_em': 0,
'tol_ls': 1e-6,
'dynamicrange': 30}


class sage2d:
    def __init__(self):
        ...
    def sage(self,y, p):
        # Get dimensions
        y=np.squeeze(y)
        [N1, N2] = y.shape

        # Setup internal parameters
        ModeFlag = 2
        ParamCount = 3
        if p > 1:
            Theta = np.zeros((ParamCount, 2),dtype=complex)
        else:
            Theta = np.zeros((ParamCount, 1),dtype=complex)
        IC_Mode = opts['ic_mode']
        MaxCycleCnt = abs(opts['maxcyclecnt'])
        DynamicRange = abs(opts['dynamicrange'])
        N1_FFT = abs(opts['n1_fft'])
        if opts['n1_fft'] != 0:
            N2_FFT = abs(opts['n2_fft'])
        else:
            N2_FFT = 2 ** self.nextpow2(1024 * N2 / N1)

        if opts['tol_em'] != 0:
            Tol = abs(opts['tol_em'])
        else:
            Tol = [1 / (2 * N1_FFT), 1 / (2 * N1_FFT), 1 / (2 * N1_FFT)]


        # Initialise parameter matrix with initial guess
        Theta_Old = np.zeros([ParamCount, p],dtype=complex)
        Theta = np.zeros([ParamCount, p],dtype=complex)
        Theta_Threshold = np.outer(Tol , np.ones((1, p),dtype=complex))
        CycleCtr = 1

        # SAGE Initialisation using FFT evaluation
        for Comp_Index in range(0, p):
            # Expectation Step: cancel all interferers
            x_i = self.ic(y, Theta, Comp_Index, 'serial')

            # Calculate initialisation cost function
            CostFunction1 = np.fft.fftshift((abs(np.fft.fft(x_i, N1_FFT, 0)) ** 2).sum(axis=1))

            # Estimate initial parameter (Maximization step)

            Dummy, Index = CostFunction1.max(0), CostFunction1.argmax(0)

            Theta[1, Comp_Index] = (1 / N1_FFT * (Index - 1)) - 0.5

            # Calculate initialisation cost function

            CostFunction2 = np.fft.fftshift(abs(np.fft.fft((np.exp(-2j * np.pi * (np.array([np.arange(0, N1)]).transpose() * np.ones([1, N2]) * Theta[1, Comp_Index])) * x_i),N2_FFT, 1).sum(axis=0)) ** 2)

            # Estimate initial parameter (Maximization step)

            Dummy, Index = CostFunction2.max(0), CostFunction2.argmax(0)

            Theta[2, Comp_Index] = (1 / N2_FFT * (Index - 1)) - 0.5

            # Initial guess of the path weight
            Theta[0, Comp_Index] = 1 / (N1 * N2) * self.cost_func([], 0, Theta[:, Comp_Index], x_i, 0)

            # Check component count
            if abs(Theta[0, Comp_Index]) < 10 ** (-DynamicRange / 20) * (abs(Theta[0, 0])):
                Theta[:, Comp_Index] = np.squeeze(np.zeros([ParamCount, 1],dtype=complex))
                print('Component %d Out of dynamic range' % Comp_Index)
        # Process until parameter converges
        CycleCtr = 1
        done=False
        while not done:

            # Process EM step for each component
            for Comp_Index in range(0, p):

                # Expectation Step: cancel all interferers
                x_i = self.ic(y, Theta, Comp_Index, IC_Mode)

                # Coordinate wise updating of the parameter vector
                Theta[1, Comp_Index] = fminbound(self.cost_func, Theta[1, Comp_Index] - 1 / (CycleCtr * N1),Theta[1, Comp_Index] + 1 / (CycleCtr * N1),(1, Theta[:, Comp_Index], x_i, ModeFlag),opts['tol_ls'], opts['maxitercnt'])

                # Coordinate wise updating of the parameter vector
                Theta[2, Comp_Index] = fminbound(self.cost_func, Theta[2, Comp_Index] - 1 / (CycleCtr * N2),Theta[2, Comp_Index] + 1 / (CycleCtr * N2),(2, Theta[:, Comp_Index], x_i, ModeFlag),opts['tol_ls'], opts['maxitercnt'])

                # Updating of the path weight
                Theta[0, Comp_Index] = 1 / (N1 * N2) * self.cost_func([], 0, Theta[:, Comp_Index], x_i, 0)

                # Check component count
                if abs(Theta[0, Comp_Index]) < 10 **(-DynamicRange / 20) * (abs(Theta[0, 0])):
                    Theta[:, Comp_Index] = np.squeeze(np.zeros([ParamCount, 1], dtype=complex))
                    done = True

                # Check for convergence (Parameter of two iterations are nearly constant)
                if np.all(abs(Theta - Theta_Old) < Theta_Threshold):
                    done = True

                else:
                    # Store current parameter vector
                    Theta_Old = Theta

                # Increment EM cycle counter
                CycleCtr += 1
                ModeFlag = 2

                # Check cycle counter
                if CycleCtr > MaxCycleCnt:
                    done=True
                    return

        # Return nonzeros components only
        Index = np.nonzero(Theta[0, :])
        Theta = Theta[:, Index]

        # Return parameters in common form
        beta = Theta[0, :].transpose()
        f1 = (Theta[1, :]).transpose()
        f2 = (Theta[2, :]).transpose()

        return beta, f1, f2, CostFunction1, CostFunction2

    def domain(self,N1,N2):
    # Setup domains

        if N1 % 2 == 0 and N2 % 2 == 0:
            n1 = np.array([np.arange(0, N1)]).transpose()
            n2 = np.array([np.arange(0, N2)])
        elif N1 % 2 == 0 and N2 % 2 != 0:
            n1 = np.array([np.arange(0, N1)]).transpose()
            n2 = np.array([np.arange(-np.floor(N2 / 2), np.floor(N2 / 2)+1)])
        elif N1 % 2 != 0 and N2 % 2 == 0:
            n1 = np.array([np.arange(-np.floor(N1 / 2), np.floor(N1 / 2)+1)]).transpose()
            n2 = np.array([np.arange(0, N2)])
        else:
            n1 = np.array([np.arange(-np.floor(N1 / 2), np.floor(N1 / 2)+1)]).transpose()
            n2 = np.array([np.arange(-np.floor(N2 / 2), np.floor(N2 / 2)+1)])
        self.n1= n1
        self.n2= n2
        return self.n1, self.n2


    def cost_func(self,ParamVal, ParamIndex, Theta_i, x_i,ModeFlag):
        if (ParamIndex != 0):
            Theta_i[ParamIndex] = ParamVal

        # Get dimensions
        [N1, N2] = np.size(x_i, 0), np.size(x_i, 1)

        # Initialise variables
        C = np.zeros((N1, N2),dtype=complex)

        f1 = Theta_i[1]
        f2 = Theta_i[2]

        # Setup domains

        [n1,n2]=self.domain(N1,N2)

        # Calculate signal model for given parameters
        C = np.exp(-2j * np.pi * f1 *n1) * np.exp(-2j * np.pi * f2 *n2) * x_i

        # Calculate discrete integral
        if ModeFlag == 0:
            # Return standard value
            c = sum(sum(C))

        elif ModeFlag == 1:

            if ParamIndex == 2:
                # Return inverted squared absolute value for initialisation
                c = -abs(sum(sum(C))) ** 2
            else:
                # Return inverted absolute value
                c = -abs(sum(sum(C))) ** 2
        else:
            # Return inverted absolute value
            c = -abs(sum(sum(C))) ** 2
        return c


    def ic(self,y, Theta, i, IC_Mode):
        # Get Number of component
        CompCount = Theta.shape
        CompCount=CompCount[1]

        # Select processing mode
        if str.lower(IC_Mode) == 'parallel':
            # Cancel all interferers
            CompVector = np.arange(0, CompCount)
            # Apart from the considered component
            CompVector = np.delete(CompVector, np.where(CompVector == i))

        elif str.lower(IC_Mode) == 'serial':
            # Cancel dominant interferers
            CompVector = np.arange(0, i)
        else:
            print('Unknown IC-Mode !')

        # Setup data
        x_i = y

        # Cancel components sequentially
        for Index in CompVector:
            # Cancel current component
            x_i = x_i - self.model_func(Theta[:, Index], y.shape)
        return x_i


    def model_func(self,Theta, Dim):

        # Initialise variables
        beta = Theta[0]
        f1 = Theta[1]
        f2 = Theta[2]

        N1 = Dim[0]
        N2 = Dim[1]
        # Setup domains

        [n1,n2]=self.domain(N1,N2)

        # Calculate signal model for given parameters
        y = beta * np.exp(2j * np.pi * f1 * n1) * np.exp(2j * np.pi * f2 * n2)

        return y

    def nextpow2(self,i):
        n = 1
        while n < i: n *= 2
        return n
