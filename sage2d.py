from numpy import *
from scipy.optimize import fminbound

"""
sage2d.py
This class has methods that extracts CSI information from NETLINK data. Can also parse
data from file, stored using log_to_file.c tool.

"""
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


class Sage2d:
    def __init__(self):
        ...

    def cost_func(self,ParamVal, ParamIndex, Theta_i, x_i,ModeFlag):
        if (ParamIndex != 0):
            Theta_i[ParamIndex] = ParamVal

        # Get dimensions
        [N1, N2] = size(x_i, 0), size(x_i, 1)

        # Initialise variables
        C = zeros([N1, N2])

        f1 = Theta_i[1]
        f2 = Theta_i[2]

        # Setup domains

        [n1,n2]=self.domain(N1,N2)

        # Calculate signal model for given parameters
        C = exp(-2j * pi * f1 *n1) * exp(-2j * pi * f2 *n2) * x_i

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


    def ic(self,y, Theta, i, IC_Mode):
        # Get Number of components
        # CompCount = Theta.size(1)
        CompCount = array([Theta]).T.shape[1]

        # Select processing mode
        if str.lower(IC_Mode) == 'parallel':
            # Cancel all interferers
            CompVector = arange(1, CompCount + 1)
            # Apart from the considered component
            CompVector = delete(CompVector, where(CompVector == i))

        elif str.lower(IC_Mode) == 'serial':
            # Cancel dominant interferers
            CompVector = arange(1, i)

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
        y = beta * exp(2j * pi * f1 * n1) * exp(2j * pi * f2 * n2)

        return y

    def nextpow2(self,i):
        n = 1
        while n < i: n *= 2
        return n

    def sage(self,y, p):
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
            N2_FFT = 2 ** self.nextpow2(1024 * N2 / N1)

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
            x_i = self.ic(y, Theta, Comp_Index, 'serial')

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
            Theta[0, Comp_Index] = 1 / (N1 * N2) * self.cost_func([], 0, Theta[:, Comp_Index], x_i, 0)

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
                x_i = self.ic(y, Theta, Comp_Index, IC_Mode)

                # Coordinate wise updating of the parameter vector
                Theta[1, Comp_Index] = fminbound(self.cost_func, Theta[1, Comp_Index] - 1 / (CycleCtr * N1),
                Theta[1, Comp_Index] + 1 / (CycleCtr * N1),
                (1, Theta[:, Comp_Index], x_i, ModeFlag),
                opts.tol_ls, opts.maxitercnt)
                # Coordinate wise updating of the parameter vector
                Theta[2, Comp_Index] = fminbound(self.cost_func, Theta[2, Comp_Index] - 1 / (CycleCtr * N2),
                Theta[2, Comp_Index] + 1 / (CycleCtr * N2),
                (1, Theta[:, Comp_Index], x_i, ModeFlag),
                opts.tol_ls, opts.maxitercnt)
                # Updating of the path weight
                Theta[1, Comp_Index] = 1 / (N1 * N2) * self.cost_func([], 0, Theta[:, Comp_Index], x_i, 0)

                # Check component count
                if abs(Theta(0, Comp_Index)) < 10 ^ (-DynamicRange / 20) * (abs(Theta(0, 0))):
                    Theta[:, Comp_Index] = zeros(ParamCount, 1)
                    # ddisp(DispMode, ['Component[',num2str(Comp_Index),']: Out of dynamic range'])
                else:
                    # Display results for current iteration

                    # ddisp(DispMode, ['Component[',num2str(Comp_Index),']: ',num2str(abs(Theta(1,Comp_Index))),' ', num2str(Theta[2:end,Comp_Index].transpose)])
                    print('ok')

                # Check for convergence (Parameter of two iterations are nearly constant)
                if all(abs(Theta - Theta_Old) < Theta_Threshold):
                    # ddisp(DispMode, ['Algorithm converged after 1+', num2str(CycleCtr),' EM-Cycles']);
                    # ddisp(DispMode, ['Mean convergence tolerance = ', num2str(sum(abs(Theta - Theta_Old)./p,2)')]);
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

    def domain(self,N1,N2):
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
        self.n1= n1
        self.n2= n2
        return self.n1, self.n2

