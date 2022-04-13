# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the rotor/stator class for the turbomachinery desing analisys 
#       

# importing libraries
import numpy as np 
import matplotlib.pyplot as plt
from geometry import bladeGenerator
from turboClass.bladeSection import section
from scipy import integrate 
from scipy import interpolate 
from scipy.misc import derivative
from turboCoeff import losses
import warnings

class blade:
    '''
    Blade object, it is used in the stage object.
        AIM:
            --- blade geometry description 
            --- blade thermodynamics 
    '''

    def __init__(self, ID, turboType, nSection, omega, inletBladeHeight, outletBladeHeight, inletHubRadius, outletHubRadius):
        '''
        Rotor object declaration: 
            variables:
                ID                  -- blade identifier
                turboType           -- blade type: stator/rotor
                nSection            -- # of blade sections 
                inletBladeHeight    -- blade inlet height
                outletBladeHeight   -- blade outlet height 
                inletHubRadius      -- inlet hub radius 
                outletHubRadius     -- outlet hub radius 
                airfoilPath         -- path where are stored the airfoil properties 
        '''

        self.ID             = ID
        self.turboType      = turboType
        self.nSection       = nSection
        self.omega          = omega
        
        # section objects allocation 
        self.inletSection = self.allocateSection(hubRadius=inletHubRadius, bladeHeight=inletBladeHeight, nSection=nSection)
        self.outletSection = self.allocateSection(hubRadius=outletHubRadius, bladeHeight=outletBladeHeight, nSection=nSection)

    def allocateSection(self, hubRadius=0, tipRadius=0, bladeHeight=0, nSection=0, plot=False):
        '''
        This function allocates the section vector in the blade object.
            inputs:
                hubRadius   -- radius of the hub 
                tipRadius   -- radius of the tip 
                bladeHeight -- blade height 
                nSection    -- # of sections the blade is composed of 
        '''

        # computing main quantities 
        sectionVec = []
        height = bladeHeight/nSection

        # loop generation for the section generation
        for ii in range(nSection):
            if hubRadius != 0:
                midpoint = hubRadius + ii * height + height / 2 
                bottom = hubRadius + ii * height
                tip = hubRadius + (ii+1) * height
            elif tipRadius != 0:
                midpoint = tipRadius - ii * height - height/2
                bottom = tipRadius - (ii+1) * height
                tip = tipRadius - ii * height 

            # appending section object to the vector 
            sectionVec.append(section(midpoint, bottom, tip, height))

        if plot:
            fig = plt.figure(figsize=(8,8))
            for ii in range(len(sectionVec)):
                plt.plot(0, sectionVec[ii].midpoint, 'r*')
                plt.plot(0, sectionVec[ii].tip, 'ob')
                plt.plot(0, sectionVec[ii].bottom, 'ok')
                plt.grid(linestyle='--')
            plt.show()

        return sectionVec 

    def plotMeridional(self):
        '''
        This function plots the blade sections in the meridional plane.
        '''

        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        import matplotlib.colors as mcolors

        # getting max and minimum axial velocity
        Vmax = 0 
        Vmin = self.outletSection[0].Va
        for ii in range(self.nSection):
            if Vmax < self.outletSection[ii].Va: 
                Vmax = self.outletSection[ii].Va
            if Vmin > self.outletSection[ii].Va:
                Vmin = self.outletSection[ii].Va

        nSection = self.nSection

        def sectionFill(inletSection, outletSection, color):
            '''
            This function plots a streamtube section.
            '''

            x = [0, 0, 0, 1, 1, 1]
            y = [inletSection.tip, inletSection.midpoint, inletSection.bottom, outletSection.bottom, outletSection.midpoint, outletSection.tip]

            plt.fill(x, y, color=color)

            if nSection < 15:
                plt.plot(0, self.inletSection[ii].midpoint, 'r*')
                plt.plot(0, self.inletSection[ii].tip, 'ob')
                plt.plot(0, self.inletSection[ii].bottom, 'ob')
                plt.plot(1, self.outletSection[ii].midpoint, 'r*')
                plt.plot(1, self.outletSection[ii].tip, 'ok')
                plt.plot(1, self.outletSection[ii].bottom, 'ok')

        fig = plt.figure(figsize=(8,8))
        for ii in range(nSection):
            sectionFill(self.inletSection[ii], self.outletSection[ii], [self.outletSection[ii].Va/Vmax, 0.2, 0.35])

        color_list = [(self.outletSection[nSection-1-ii].Va/Vmax, 0.2, 0.35) for ii in range(nSection)]
        cmap = mcolors.LinearSegmentedColormap.from_list("my_colormap", color_list)
        cmappable = ScalarMappable(norm=Normalize(Vmin, Vmax), cmap=cmap)

        plt.grid(linestyle='--')
        plt.title('Axial velocity')
        plt.xlabel('chord')
        plt.ylabel('r')
        plt.colorbar(cmappable, label=r'$\frac{m}{s}$')
        plt.show()

    def printMeridional(self):
        '''
        This function plots the main inlet and outlet quantities of the blade.
        '''

        fig, ax = plt.subplots(ncols=2, nrows=1)
        fig.set_figwidth(18)
        fig.set_figheight(9.5)

        # inlet section plot 
        twiny1 = ax[0].twiny()
        twiny2 = ax[0].twiny() 
        twiny3 = ax[0].twiny()
        # moving axis position
        twiny1.spines["top"].set_position(("axes", 1))
        twiny2.spines["top"].set_position(("axes", 1.2))
        twiny3.spines["top"].set_position(("axes", 1.1))
        
        # vector declaration
        midpoint = np.zeros(self.nSection)
        Va = np.zeros(self.nSection)
        Vt = np.zeros(self.nSection)
        Pt = np.zeros(self.nSection)
        s = np.zeros(self.nSection)

        # vector allocation
        for ii in range(self.nSection):
            midpoint[ii] = self.inletSection[ii].midpoint 
            Va[ii]       = self.inletSection[ii].Va 
            Vt[ii]       = self.inletSection[ii].Vt 
            Pt[ii]       = self.inletSection[ii].Pt / 1e+5
            s[ii]        = self.inletSection[ii].s 

        p0, = ax[0].plot(Va, midpoint, 'k', label=r'$V_a$')
        p1, = twiny1.plot(Vt, midpoint, 'm', label=r'$V_t$')
        p2, = twiny2.plot(Pt, midpoint, 'b--', label=r'$P_t$')
        p3, = twiny3.plot(s, midpoint, 'r--', label=r'$s$')
        ax[0].set_ylim(self.inletSection[0].midpoint, self.inletSection[-1].midpoint)
        ax[0].set_ylabel('r')
        ax[0].set_xlabel(r'$V_a$')
        twiny1.set_xlabel(r'$V_t$')
        twiny2.set_xlabel(r'$P_t$')
        twiny3.set_xlabel(r'$s$')

        ax[0].legend(handles=[p0,p1,p2,p3])
        ax[0].set_title('Inlet')

        # outlet section plot 
        twiny21 = ax[1].twiny()
        twiny22 = ax[1].twiny() 
        twiny23 = ax[1].twiny()
        # moving axis position
        twiny21.spines["top"].set_position(("axes", 1))
        twiny22.spines["top"].set_position(("axes", 1.2))
        twiny23.spines["top"].set_position(("axes", 1.1))

        # vector allocation
        for ii in range(self.nSection):
            midpoint[ii] = self.outletSection[ii].midpoint 
            Va[ii]       = self.outletSection[ii].Va 
            Vt[ii]       = self.outletSection[ii].Vt 
            Pt[ii]       = self.outletSection[ii].Pt / 1e+5
            s[ii]        = self.outletSection[ii].s 

        p20, = ax[1].plot(Va, midpoint, 'k', label=r'$V_a$')
        p21, = twiny21.plot(Vt, midpoint, 'm', label=r'$V_t$')
        p22, = twiny22.plot(Pt, midpoint, 'b--', label=r'$P_t$')
        p23, = twiny23.plot(s, midpoint, 'r--', label=r'$s$')
        ax[1].set_ylim(self.outletSection[0].midpoint, self.outletSection[-1].midpoint)
        ax[1].set_ylabel('r')
        ax[1].set_xlabel(r'$V_a$')
        twiny21.set_xlabel(r'$V_t$')
        twiny22.set_xlabel(r'$P_t$')
        twiny23.set_xlabel(r'$s$')

        ax[1].legend(handles=[p20,p21,p22,p23])
        ax[1].set_title('Outlet')

        plt.tight_layout()
        plt.show()

    def allocateDynamics(self, rMean, VtMean, VaMean, omega, section='inlet'):
        '''
        This function allocates the velocity vectors at each section points using the FREE VORTEX model.
            inputs:
                rMean   -- mean radius
                VtMean  -- tangential mean velocity 
                VaMean  -- axial mean velocity 
                omega   -- angular velocity
                section -- inlet/outlet section 
        '''

        if section == 'inlet':
            for ii in range(self.nSection):
                # tangential speed computation with respect to the FREE VORTEX model 
                Vt = VtMean * rMean / self.inletSection[ii].midpoint
                # rotation speed 
                U = self.inletSection[ii].midpoint * omega 

                # data allocation in section object
                self.inletSection[ii].allocateDynamics(VaMean, Vt, U)
        elif section == 'outlet':
            for ii in range(self.nSection):
                # tangential speed computation with respect to the FREE VORTEX model 
                Vt = VtMean * rMean / self.outletSection[ii].midpoint
                # rotation speed 
                U = self.outletSection[ii].midpoint * omega 
                
                # data allocation in section object
                self.outletSection[ii].allocateDynamics(VaMean, Vt, U)

    def allocateThermodynamics(self, Tt0, Pt0, eta, R=287.06, gamma=1.4, printout=False):
        '''
        This function allocates the thermodynamic properties to each section.
            inputs:
                Tt0     -- inlet total temperature 
                Pt0     -- inlet total pressure 
                Leu     -- real euler work  
                eta     -- stage efficiency   
        '''

        # cP computation
        cP = gamma / (gamma - 1) * R

        for ii in range(self.nSection):
            # temperature computation 
            T0 = Tt0 - self.inletSection[ii].V**2 / (2 * cP)

            # pressure computation
            P0 = Pt0 * (T0/Tt0)**(gamma/(gamma-1))

            # density computation 
            rho0 = P0 / (R * T0)

            # total density computation
            rhot0 = Pt0 / (R * Tt0)

            # Tt0r computation
            Tt0r = Tt0 - self.inletSection[ii].V**2 / (2 * cP) + self.inletSection[ii].W**2 / (2 * cP)

            # relative total pressure computation
            Pt0r = Pt0 * (Tt0r / Tt0)**(gamma/gamma-1)

            # euler work computation 
            if self.turboType == 'stator':
                Leu = 0 
            else: 
                Leu = self.outletSection[ii].U * self.outletSection[ii].Vt - self.inletSection[ii].U * self.inletSection[ii].Vt

            # total temperature computation
            Tt1 = Leu / cP + Tt0

            # ideal temperature computation if the process is completely 
            # isentropic without losses but the work produced is related 
            # to a process that takes into account losses in the stage  
            T1 = Tt1 - self.outletSection[ii].V**2 / (2 * cP)

            # T1/Tt1 isentropic computation 
            #   this correction activates only if eta != 1
            T1iso = T0 + eta * (T1 - T0)
            Tt1iso = Tt0 + eta * (Tt1 - Tt0)

            # pressure computation
            P1 = Pt0 * (T1iso/Tt0)**(gamma/(gamma-1))

            # total pressure computation
            Pt1 = Pt0 * (Tt1iso/Tt0)**(gamma/(gamma-1))

            # Ttr1 computation
            Tt1r = Tt1 - self.outletSection[ii].V**2 / (2 * cP) + self.outletSection[ii].W**2 / (2 * cP)

            # relative total pressure computation
            Pt1r = Pt1 * (Tt1r / Tt1)**(gamma/gamma-1)

            # density computation 
            rho1 = P1 / (R * T1)

            # total density computation
            rhot1 = Pt1 / (R * Tt1)

            if printout:
                print('midpont = ', self.inletSection[ii].midpoint)
                print('Uin  = ', self.inletSection[ii].U)
                print('Vain = ', self.inletSection[ii].Va)
                print('Vtin  = ', self.inletSection[ii].Vt)
                print('Vin = ', self.inletSection[ii].V)
                print('Uout = ', self.outletSection[ii].U)
                print('Vaout = ', self.outletSection[ii].Va)
                print('Vtout = ', self.outletSection[ii].Vt)
                print('Vout = ', self.outletSection[ii].V)
                print('Leu  = ', Leu)
                print('T0 = ', T0)
                print('P0 = ', P0)
                print('Tt0 = ', Tt0)
                print('Pt0 = ', Pt0)
                print('Tt0r = ', Tt0r)
                print('Pt0r = ', Pt0r)
                print('T1 = ', T1)
                print('P1 = ', P1)
                print('Tt1 = ', Tt1)
                print('Pt1 = ', Pt1)
                print('Tt1r = ', Tt1r)
                print('Pt1r = ', Pt1r)
                print('')

            # variable allocation in seection objects            
            self.inletSection[ii].allocateThermodynamics(Tt=Tt0, Pt=Pt0, T=T0, P=P0, Ttr=Tt0r, Ptr=Pt0r, rho=rho0, rhot=rhot0, s=0)
            self.outletSection[ii].allocateThermodynamics(Tt=Tt1, Pt=Pt1, T=T1, P=P1, Ttr=Tt1r, Ptr=Pt1r, rho=rho1, rhot=rhot1, s=0)

    def radialEquilibrium(self, Pt0, Tt0, eta, mFlux, plot=False, R=287.06, gamma=1.4):
        '''
        This function computes the radial equilibrium of the section taking into account losses. 
            inputs:
        ''' 
    
        # cP computation
        cP = gamma / (gamma - 1) * R

        # INLET VARIABLES INTERPOLATION -> they do not change so there is not the necessity to change them in the while loop
        # inlet midpoint generation 
        midpointInlet = [self.inletSection[ii].midpoint for ii in range(self.nSection)]

        # s1 function generation --> blade inlet 
        s1 = [self.inletSection[ii].s for ii in range(self.nSection)]
        s1 = interpolate.interp1d(midpointInlet, s1, kind='linear', bounds_error=False, fill_value='extrapolate')

        # Va1 function generation --> blade inlet 
        Va1 = [self.inletSection[ii].Va for ii in range(self.nSection)]
        Va1 = interpolate.interp1d(midpointInlet, Va1, kind='linear', bounds_error=False, fill_value='extrapolate')

        # Vt1 function generation --> blade inlet 
        rVt1 = [self.inletSection[ii].Vt * self.inletSection[ii].midpoint for ii in range(self.nSection)]
        rVt1 = interpolate.interp1d(midpointInlet, rVt1, kind='linear', bounds_error=False, fill_value='extrapolate')

        # Tt1 function generation --> blade inlet 
        Tt1 = [self.inletSection[ii].Tt for ii in range(self.nSection)]
        Tt1 = interpolate.interp1d(midpointInlet, Tt1, kind='linear', bounds_error=False, fill_value='extrapolate')

        # OUTLET VARIABLES -> they do not change in the while loop so there is not the necessity to change them
        # outlet midpoint generation 
        midpointOutlet = [self.outletSection[ii].midpoint for ii in range(self.nSection)]

        # rVt2 function generation --> blade outlet 
        rVt2 = [self.outletSection[ii].Vt * self.outletSection[ii].midpoint for ii in range(self.nSection)]
        rVt2 = interpolate.interp1d(midpointOutlet, rVt2, kind='linear', bounds_error=False, fill_value='extrapolate')

        # Vt2 function generation --> blade outlet 
        Vt2 = [self.outletSection[ii].Vt for ii in range(self.nSection)]
        Vt2 = interpolate.interp1d(midpointOutlet, Vt2, kind='linear', bounds_error=False, fill_value='extrapolate')

        # ODE function necessary tools
        # omega 
        omega = self.omega 

        # Va2 function generation 
        def radialFunc(y, t, s2):
            '''
            Radial equilibrium ODE. 
                y = Va2**2
                t = r 
            '''

            # derivative step setup
            dx = 1e-3
            # derivative computation 
            dTt1 = derivative(Tt1, t, dx=dx, n=1, order=3)
            drVt2 = derivative(rVt2, t, dx=dx, n=1, order=3)
            drVt1 = derivative(rVt1, t, dx=dx, n=1, order=3)
            ds2 = derivative(s2, t, dx=dx, n=1, order=3)

            # y derivative computation
            dydt = - 2 * (- y/(2*cP) * ds2 - cP * dTt1 - omega * drVt2 + omega * drVt1 + Tt1(t) * ds2 + omega / cP * rVt2(t) * ds2 - omega / cP * rVt1(t) * ds2 - Vt2(t)**2 / (2*cP) * ds2 + Vt2(t) / t * drVt2)

            return dydt 
        
        # radius vector allocation 
        t = midpointOutlet #np.linspace(midpointOutlet[0], midpointOutlet[-1], 1000)

        # setting up boundary condition
        y0 = self.outletSection[0].Va**2

        # setting up tolerances
        relError = 1 
        tol      = 1e-2
        nMax     = 100 
        counter  = 0

        # initialize print
        lineLenght = 80
        iterativeLenght = np.int16((lineLenght - len(' RADIAL EQ. '))/2)
        print('\n' + '*' * iterativeLenght + ' RADIAL EQ. ' + '*' * iterativeLenght)
        
        # LOOP
        while relError > tol and counter < nMax:
            # updating counter 
            counter = counter + 1
              
            # OUTLET VARIABLES INTERPOLATION 
            # s2 function generation --> blade outlet 
            s2 = [self.outletSection[ii].s for ii in range(self.nSection)]
            s2 = interpolate.interp1d(midpointOutlet, s2, kind='linear', bounds_error=False, fill_value='extrapolate')

            # computing solution 
            Va2_squared = integrate.odeint(radialFunc, y0, t, args=(s2,))

            # setting up Va2 vector 
            Va2 = np.zeros(self.nSection)
            for ii in range(self.nSection):
                Va2[ii] = np.sqrt(Va2_squared[ii])

            # outlet section dynamics allocation
            for ii in range(self.nSection):
                self.outletSection[ii].allocateDynamics(Va2[ii], self.outletSection[ii].Vt, self.outletSection[ii].U)

            # getting info on s 
            s2Vec = [self.outletSection[ii].s for ii in range(self.nSection)]
            s2Min = np.min(s2Vec)
            s2Max = np.max(s2Vec)
            s2Ave = np.sum(s2Vec)/self.nSection

            # getting info on Pt
            Pt2Vec = [self.outletSection[ii].Pt for ii in range(self.nSection)]
            Pt2Min = np.min(Pt2Vec)
            Pt2Max = np.max(Pt2Vec)
            Pt2Ave = np.sum(Pt2Vec)/self.nSection

            # getting info on losses 
            try:
                lossMax = np.max(lossVec)
                lossMin = np.min(lossVec)
                lossAve = np.sum(lossVec)/self.nSection
            except:
                pass

            # computing rho for the mass flux 
            self.allocateThermodynamics(Tt0=Tt0, Pt0=Pt0, eta=1)

            # check mass flux 
            newFlux = 0
            for ii in range(self.nSection):
                newFlux = newFlux + self.outletSection[ii].mFlux()

            # relative mass flux error
            relError = np.abs(newFlux - mFlux)/mFlux

            # setting up initial velocity for the new ODE 
            if newFlux > mFlux:
                y0 = y0 * (1 - relError)
            else: 
                y0 = y0 * (1 + relError)

            # printing main values
            iterativeLenght = np.int16((lineLenght - len(' ITERATION   '))/2)
            print('*' * iterativeLenght + ' ITERATION {0:d} '.format(counter) + '*' * iterativeLenght)
            print('-- mFlux   {0:>8.2f} kg/s -- iterFlux {1:>8.2f} kg/s -- rel. error {2:>8.2f}'.format(mFlux,newFlux,relError))
            try:
                print('-- lossMin {0:>8.2f}      -- lossAve  {1:>8.2f}      -- lossMax    {2:>8.2f}'.format(lossMin,lossAve,lossMax))
            except: 
                pass
            print('-- s2_min  {0:>8.2f} J/kg -- s2_ave   {1:>8.2f} J/kg -- s2_max     {2:>8.2f} J/kg'.format(s2Min,s2Ave,s2Max))
            print('-- Pt2_min {0:>8.2f} bar  -- Pt2_ave  {1:>8.2f} bar  -- Pt2_max    {2:>8.2f} bar'.format(Pt2Min/1e+5,Pt2Ave/1e+5,Pt2Max/1e+5))

            # computing losses for the new loop cycle (used if relError > tol) 
            # target s2
            # LOSSES COMPUTATION 
            # computing pressure losses 
            lossVec = np.zeros(self.nSection)

            for ii in range(self.nSection):
                # variables allocation -> these variables are already updated with the previous self.allocateDynamics() function call 
                # velocity 
                W1 = self.inletSection[ii].W
                W2 = self.outletSection[ii].W

                # angles
                beta1 = self.inletSection[ii].beta
                beta2 =  self.outletSection[ii].beta
                
                # angles check -> the Leiblein model treats with positive angle for beta1 
                # this implies changing the angle sign with respect to beta1 sign angle  
                if beta1 < 0:
                    beta1 = - beta1 
                    beta2 = - beta2

                # solidity allocation    
                solidity = 1

                # loss computation -> Lieblein model 
                lossVec[ii], _ = losses.lossCoeff(W1, W2, beta1, beta2, solidity)

            # computing entropy generation through losses 
            for ii in range(self.nSection):
                # new outlet pressure computation -> using losses
                self.outletSection[ii].Ptr = self.inletSection[ii].Ptr - lossVec[ii] * (self.inletSection[ii].Ptr - self.inletSection[ii].P)
                # entropy computation
                self.outletSection[ii].s = self.inletSection[ii].s - R * np.log(self.outletSection[ii].Ptr / self.inletSection[ii].Ptr)                        

        print('*' * lineLenght)

        # plotting interpolated functions 
        if plot:
            midpointInlet = [self.inletSection[ii].midpoint for ii in range(self.nSection)]
            midpointOutlet = [self.outletSection[ii].midpoint for ii in range(self.nSection)]    

            fig, ax = plt.subplots(ncols=2, nrows=1)
            fig.set_figwidth(18)
            fig.set_figheight(9.5)

            # inlet section plot 
            twiny1 = ax[0].twiny()
            twiny2 = ax[0].twiny() 
            # moving axis position
            twiny1.spines["top"].set_position(("axes", 1))
            twiny2.spines["top"].set_position(("axes", 1.1))

            # Vt1 allocation 
            Vt1 = [self.inletSection[ii].Vt for ii in range(self.nSection)]

            p0, = ax[0].plot(Va1(midpointInlet), midpointInlet, 'ko-', label=r'$V_a$')
            p1, = twiny1.plot(Vt1, midpointInlet, 'm*-', label=r'$V_t$')
            p2, = twiny2.plot(s1(midpointInlet), midpointInlet, 'r--', label=r'$s$')
            ax[0].set_ylim(self.inletSection[0].midpoint, self.inletSection[-1].midpoint)
            ax[0].set_ylabel('r')
            ax[0].set_xlabel(r'$V_a$')
            twiny1.set_xlabel(r'$V_t$')
            twiny2.set_xlabel(r'$s$')

            ax[0].legend(handles=[p0,p1,p2])
            ax[0].set_title('Inlet')

            # inlet section plot 
            twiny1 = ax[1].twiny()
            twiny2 = ax[1].twiny() 
            twiny3 = ax[1].twiny()
            # moving axis position
            twiny1.spines["top"].set_position(("axes", 1))
            twiny2.spines["top"].set_position(("axes", 1.1))
            twiny3.spines["top"].set_position(("axes", 1.2))

            # Vt2 allocation 
            Vt2 = [self.outletSection[ii].Vt for ii in range(self.nSection)]
            # Pt2 allocation 
            Pt2 = [self.outletSection[ii].Pt for ii in range(self.nSection)]

            p0, = ax[1].plot(Va2, midpointInlet, 'ko-', label=r'$V_a$')
            p1, = twiny1.plot(Vt2, midpointInlet, 'm*-', label=r'$V_t$')
            p2, = twiny2.plot(s2(midpointOutlet), midpointInlet, 'r--', label=r'$s$')
            p3, = twiny3.plot(Pt2, midpointOutlet, 'b-o', label=r'$P_t$')
            ax[1].set_ylim(self.inletSection[0].midpoint, self.inletSection[-1].midpoint)
            ax[1].set_ylabel('r')
            ax[1].set_xlabel(r'$V_a$')
            twiny1.set_xlabel(r'$V_t$')
            twiny2.set_xlabel(r'$s$')
            twiny3.set_xlabel(r'$P_t$')

            ax[1].legend(handles=[p0,p1,p2], loc='lower left')
            ax[1].set_title('Outlet')

            plt.tight_layout()
            plt.show()

#    # outlet thermodynamics computation
#    self.allocateThermodynamics(Tt0, Pt0, eta)
#    # LOSSES COMPUTATION 
#    # computing pressure losses 
#    lossVec = np.zeros(self.nSection)
#
#    for ii in range(self.nSection):
#        # variables allocation
#        # velocity 
#        W1 = self.inletSection[ii].W
#        W2 = self.outletSection[ii].W
#
#        # angles
#        beta1 = self.inletSection[ii].beta
#        beta2 =  self.outletSection[ii].beta
#        
#        # angles check -> the Leiblein model treats with positive angle for beta1 
#        # this implies changing the angle sign with respect to beta1 sign angle  
#        if beta1 < 0:
#            beta1 = - beta1 
#            beta2 = - beta2
#
#        # solidity allocation    
#        solidity = 1
#
#        # loss computation -> Lieblein model 
#        lossVec[ii], _ = losses.lossCoeff(W1, W2, beta1, beta2, solidity)
#
#    # computing change in total pressure due to losses 
#    for ii in range(self.nSection):
#        # new outlet pressure computation -> using losses
#        self.outletSection[ii].Ptr = self.inletSection[ii].Ptr - lossVec[ii] * (self.inletSection[ii].Ptr - self.inletSection[ii].P)
#        # entropy computation
#        self.outletSection[ii].s = self.inletSection[ii].s - R * np.log(self.outletSection[ii].Ptr / self.inletSection[ii].Ptr)
#        # total pressure compuation 
#        print('delta s / R =', (self.inletSection[ii].s - self.outletSection[ii].s) / R)
#        print(- self.outletSection[ii].s)
#        print('inlet ptr',self.inletSection[ii].Ptr)
#        print('outlet ptr',self.outletSection[ii].Ptr)
#        self.outletSection[ii].Pt = self.inletSection[ii].Pt * np.exp((self.inletSection[ii].s - self.outletSection[ii].s) / R) 
#        # static temperature computation 
#        self.outletSection[ii].T = self.outletSection[ii].Tt - self.outletSection[ii].V**2 / (2 * cP)
#        # total relative temperature computation
#        self.outletSection[ii].Ttr = self.outletSection[ii].Tt - self.outletSection[ii].W**2 / (2 * cP)
#        # static pressure computation
#        self.outletSection[ii].P = self.outletSection[ii].Pt * (self.outletSection[ii].Tt/self.outletSection[ii].T)**(gamma/(gamma-1))
#        # static density computation
#        self.outletSection[ii].rho = self.outletSection[ii].P / (R * self.outletSection[ii].T)
#        # total density computation
#        self.outletSection[ii].rhot = self.outletSection[ii].Pt / (R * self.outletSection[ii].Tt)
#
#