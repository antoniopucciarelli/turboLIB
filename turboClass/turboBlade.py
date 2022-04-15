# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the rotor/stator class for the turbomachinery desing analisys 
#       

# importing libraries
import numpy                 as np 
import matplotlib.pyplot     as plt
from scipy                   import integrate, interpolate
from scipy.misc              import derivative
from turboCoeff              import lieblein, losses
from turboClass.bladeSection import section

class blade:
    '''
    Blade object, it is used in the stage object.
        AIM:
            --- blade geometry description 
            --- blade thermodynamics 
    '''

    def __init__(self, ID, turboType, nSection, omega, nBlade, inletBladeHeight, outletBladeHeight, inletHubRadius, outletHubRadius):
        '''
        Rotor object declaration: 
            variables:
                ID                  -- blade identifier
                turboType           -- blade type: stator/rotor
                nSection            -- # of blade sections 
                omega               -- angular velocity 
                nBlade              -- # of blades 
                inletBladeHeight    -- blade inlet height
                outletBladeHeight   -- blade outlet height 
                inletHubRadius      -- inlet hub radius 
                outletHubRadius     -- outlet hub radius 
                airfoilPath         -- path where are stored the airfoil properties 
        '''

        self.ID                 = ID
        self.turboType          = turboType
        self.nSection           = nSection
        self.omega              = omega
        self.nBlade             = nBlade
        self.inletBladeHeight   = inletBladeHeight
        self.outletBladeHeight  = outletBladeHeight
        
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
                pitch = 2 * np.pi * midpoint / self.nBlade
            elif tipRadius != 0:
                midpoint = tipRadius - ii * height - height/2
                bottom = tipRadius - (ii+1) * height
                tip = tipRadius - ii * height 
                pitch = 2 * np.pi * midpoint / self.nBlade

            # appending section object to the vector 
            sectionVec.append(section(midpoint, bottom, tip, height, pitch))

        if plot:
            fig = plt.figure(figsize=(8,8))
            for ii in range(len(sectionVec)):
                plt.plot(0, sectionVec[ii].midpoint, 'r*')
                plt.plot(0, sectionVec[ii].tip, 'ob')
                plt.plot(0, sectionVec[ii].bottom, 'ok')
                plt.plot(0, sectionVec[ii].pitch, 'og')
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

    def allocateKinetics(self, rMean, VtMean, VaMean, omega, section='inlet'):
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
                self.inletSection[ii].allocateKinetics(VaMean, Vt, U)
        elif section == 'outlet':
            for ii in range(self.nSection):
                # tangential speed computation with respect to the FREE VORTEX model 
                Vt = VtMean * rMean / self.outletSection[ii].midpoint
                # rotation speed 
                U = self.outletSection[ii].midpoint * omega 
                
                # data allocation in section object
                self.outletSection[ii].allocateKinetics(VaMean, Vt, U)

    def allocateThermodynamics(self, Tt0, Pt0, eta, R=287.06, gamma=1.4):
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
            Pt0r = Pt0 * (Tt0r / Tt0)**(gamma/(gamma-1))

            # computing total relative density 
            rhot0r = Pt0r / (R * Tt0r)

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
            Pt1r = Pt1 * (Tt1r / Tt1)**(gamma/(gamma-1))

            # density computation 
            rho1 = P1 / (R * T1)

            # total density computation
            rhot1 = Pt1 / (R * Tt1)

            # computing total relative density 
            rhot1r = Pt1r / (R * Tt1r)

            # variable allocation in section objects            
            self.inletSection[ii].allocateThermodynamics(Tt=Tt0, Pt=Pt0, T=T0, P=P0, Ttr=Tt0r, Ptr=Pt0r, rho=rho0, rhot=rhot0, rhotr=rhot0r, s=self.inletSection[ii].s)
            self.outletSection[ii].allocateThermodynamics(Tt=Tt1, Pt=Pt1, T=T1, P=P1, Ttr=Tt1r, Ptr=Pt1r, rho=rho1, rhot=rhot1, rhotr=rhot1r, s=self.outletSection[ii].s)

    def radialEquilibrium(self, Pt0, Tt0, mFlux, plot=False, save=False, nMaxS=100, nMaxFlux=100, tolS=1e-2, tolFlux=1e-2, position0='inOut.pgf', position1='betaInbetaOut.pgf', R=287.06, gamma=1.4):
        '''
        This function computes the radial equilibrium of the section taking into account losses. 
            inputs:
                Pt0     -- inlet total pressure 
                Tt0     -- inlet total temperature 
                eta     -- stage efficiency 
                plot    -- boolean value for the plotting of the results
                save    -- boolean value for the saving of the results
                R       -- gas constant 
                gamma   -- specific heat ratio 
        ''' 
    
        # cP computation
        cP = gamma / (gamma - 1) * R

        # omega 
        omega = self.omega 

        # saving initial beta2 angle for plotting 
        beta2start = [self.outletSection[ii].beta for ii in range(self.nSection)]

        # INLET VARIABLES INTERPOLATION 
        # these variables do not change so they are set once in all the process
        # midpoint --> blade inlet  
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

        # OUTLET VARIABLES 
        # these variables do not change so they are set once in all the process 
        # midpoint --> blade outlet
        midpointOutlet = [self.outletSection[ii].midpoint for ii in range(self.nSection)]

        # rVt2 function generation --> blade outlet 
        rVt2 = [self.outletSection[ii].Vt * self.outletSection[ii].midpoint for ii in range(self.nSection)]
        rVt2 = interpolate.interp1d(midpointOutlet, rVt2, kind='linear', bounds_error=False, fill_value='extrapolate')

        # Vt2 function generation --> blade outlet 
        Vt2 = [self.outletSection[ii].Vt for ii in range(self.nSection)]
        Vt2 = interpolate.interp1d(midpointOutlet, Vt2, kind='linear', bounds_error=False, fill_value='extrapolate')

        # setting up entropy tolerances and counters
        # iteration counter
        counterS = 0
        # mass flux 
        relErrorS = 1 

        # initialize print
        lineLenght = 80
        iterativeLenght = np.int16((lineLenght - len(' RADIAL EQ. '))/2)
        print('' + '*' * iterativeLenght + ' RADIAL EQ. ' + '*' * iterativeLenght)

        # entropy outer loop 
        while relErrorS > tolS and counterS < nMaxS: 
            # updating entropy loop counter 
            counterS = counterS + 1

            # print
            outerIterativeLenght = np.int16((lineLenght - len(' OUTER ITERATION   '))/2)
            print('\n' + '*' * outerIterativeLenght + ' OUTER ITERATION {0:d} '.format(counterS) + '*' * outerIterativeLenght)

            # ODE function necessary tools
            # OUTLET VARIABLES INTERPOLATION => this function changes with respect to the Va2 changes 
            # s2 is updated with as soon as Va2 reached convergence
            # s2 function generation --> blade outlet 
            s2 = [self.outletSection[ii].s for ii in range(self.nSection)]
            s2 = interpolate.interp1d(midpointOutlet, s2, kind='linear', bounds_error=False, fill_value='extrapolate')
            
            # Va2 function generation 
            def radialFunc(y, t):
                '''
                Radial equilibrium ODE. 
                    y = Va2**2
                    t = r 
                '''

                # derivative step setup
                dx = 1e-3

                # derivative computation
                # dTt1 / dr  
                dTt1 = derivative(Tt1, t, dx=dx, n=1, order=3)
                # drVt2 / dr 
                drVt2 = derivative(rVt2, t, dx=dx, n=1, order=3)
                # drVt1 / dr
                drVt1 = derivative(rVt1, t, dx=dx, n=1, order=3)
                # ds2 / dr
                ds2 = derivative(s2, t, dx=dx, n=1, order=3)

                # y derivative computation
                dydt = - 2 * (- y / (2*cP) * ds2 - cP * dTt1 - omega * drVt2 + omega * drVt1 + Tt1(t) * ds2 + omega / cP * rVt2(t) * ds2 - omega / cP * rVt1(t) * ds2 - Vt2(t)**2 / (2*cP) * ds2 + Vt2(t) / t * drVt2)

                return dydt 
            
            # radius vector allocation 
            t = midpointOutlet 

            # setting up boundary condition
            y0 = self.outletSection[0].Va**2

            # setting up tolerances and counters 
            # iterations counter
            counterFlux = 0
            # mass flux 
            relErrorFlux = 1 

            # computing losses for the new loop cycle (used if relError > tol) 
            # LOSSES COMPUTATION => target s2
            # computing pressure losses 
            lossVec = self.computeLosses()

            # computing entropy generation through losses 
            for ii in range(self.nSection):
                # storing s2 old 
                s2old = [self.outletSection[ii].s for ii in range(self.nSection)]
                
                # COMPUTING NEW s VALUES
                # new outlet pressure computation -> using losses
                self.outletSection[ii].Ptr = self.inletSection[ii].Ptr - lossVec[ii] * (self.inletSection[ii].Ptr - self.inletSection[ii].P)
                # entropy computation
                self.outletSection[ii].s = self.inletSection[ii].s - R * np.log(self.outletSection[ii].Ptr / self.inletSection[ii].Ptr)                        

                # storing new s2 
                s2new = [self.outletSection[ii].s for ii in range(self.nSection)]

            # loss error computation 
            relErrorS = np.abs(np.sum(s2old) - np.sum(s2new)) / np.abs(sum(s2old))
        
            # getting info on s 
            s2Vec = [self.outletSection[ii].s for ii in range(self.nSection)]
            s2Min = np.min(s2Vec)
            s2Max = np.max(s2Vec)

            # getting info on losses 
            try:
                lossMax = np.max(lossVec)
                lossMin = np.min(lossVec)
                lossAve = np.sum(lossVec)/self.nSection
            except:
                pass

            # continuity inner loop
            while relErrorFlux > tolFlux and counterFlux < nMaxFlux:
                # updating counter 
                counterFlux = counterFlux + 1

                # computing solution 
                Va2_squared = integrate.odeint(radialFunc, y0, t)

                # setting up Va2 vector 
                Va2 = np.zeros(self.nSection)
                for ii in range(self.nSection):
                    Va2[ii] = np.sqrt(Va2_squared[ii])

                # outlet section dynamics allocation
                for ii in range(self.nSection):
                    self.outletSection[ii].allocateKinetics(Va2[ii], self.outletSection[ii].Vt, self.outletSection[ii].U)

                # computing rho for the mass flux 
                self.allocateThermodynamics(Tt0=Tt0, Pt0=Pt0, eta=1)

                # check mass flux 
                newFlux = 0
                for ii in range(self.nSection):
                    newFlux = newFlux + self.outletSection[ii].mFlux()

                # relative mass flux error
                relErrorFlux = np.abs(newFlux - mFlux) / mFlux

                # setting up initial velocity for the new ODE 
                if newFlux > mFlux:
                    y0 = y0 * (1 - relErrorFlux)
                else: 
                    y0 = y0 * (1 + relErrorFlux)

                # printing main values
                innerIterativeLenght = np.int16((lineLenght - len(' INNER ITERATION   '))/2)
                print('*' * innerIterativeLenght + ' INNER ITERATION {0:d} '.format(counterFlux) + '*' * innerIterativeLenght)
                print('-- mFlux   {0:>8.2f} kg/s -- iterFlux {1:>8.2f} kg/s -- rel. error flux {2:>2.4f}'.format(mFlux,newFlux,relErrorFlux))
                try:
                    print('-- lossMin {0:>8.2f}      -- lossAve  {1:>8.2f}      -- lossMax         {2:>2.2f}'.format(lossMin,lossAve,lossMax))
                    print('-- s2_min  {0:>8.2f} J/kg -- s2_max   {1:>8.2f} J/kg -- rel. error s    {2:>2.4f}'.format(s2Min,s2Max,relErrorS))
                except: 
                    pass
            
            print('*' * lineLenght)
            # reallocating the entropy for the next iteration 
            for ii in range(self.nSection):
                self.outletSection[ii].s = s2new[ii]

        # computing all the new thermodynamic quantities after the radial equilibrium is satisfied
        # computing change in total pressure due to losses 
        for ii in range(self.nSection):
            # total pressure compuation -- wrong
            self.outletSection[ii].Pt = self.outletSection[ii].Pt * np.exp((self.inletSection[ii].s - self.outletSection[ii].s) / R) 
            # static temperature computation 
            self.outletSection[ii].T = self.outletSection[ii].Tt - self.outletSection[ii].V**2 / (2 * cP)
            # total relative temperature computation
            self.outletSection[ii].Ttr = self.outletSection[ii].Tt - (self.outletSection[ii].V**2 - self.outletSection[ii].W**2) / (2 * cP)
            # static pressure computation
            self.outletSection[ii].P = self.outletSection[ii].Pt * (self.outletSection[ii].T/self.outletSection[ii].Tt)**(gamma/(gamma-1))
            # static density computation
            self.outletSection[ii].rho = self.outletSection[ii].P / (R * self.outletSection[ii].T)
            # total density computation
            self.outletSection[ii].rhot = self.outletSection[ii].Pt / (R * self.outletSection[ii].Tt)
            # total relative density computation 
            self.outletSection[ii].rhotr = self.outletSection[ii].Ptr / (R * self.outletSection[ii].Ttr)

        # plotting interpolated functions 
        if plot or save:
            if save:
                # setting matplotlib LaTeX export 
                import matplotlib
                matplotlib.use("pgf")
                matplotlib.rcParams.update({
                    "pgf.texsystem": "pdflatex",
                    'font.family': 'serif',
                    'text.usetex': True,
                    'pgf.rcfonts': False,
                })

            midpointInlet = [self.inletSection[ii].midpoint for ii in range(self.nSection)]
            midpointOutlet = [self.outletSection[ii].midpoint for ii in range(self.nSection)]    

            fig0, ax = plt.subplots(ncols=2, nrows=1)
            fig0.set_figwidth(18)
            fig0.set_figheight(9.5)

            # inlet section plot 
            twiny1 = ax[0].twiny()
            twiny2 = ax[0].twiny() 
            twiny3 = ax[0].twiny()
            # moving axis position
            twiny1.spines["top"].set_position(("axes", 1))
            twiny2.spines["top"].set_position(("axes", 1.1))
            twiny3.spines["top"].set_position(("axes", 1.2))

            # Vt1 allocation 
            Vt1 = [self.inletSection[ii].Vt for ii in range(self.nSection)]
            # Pt1 allocation 
            Pt1 = [self.inletSection[ii].Pt/1e+5 for ii in range(self.nSection)]

            # plotting data
            p0, = ax[0].plot(Va1(midpointInlet), midpointInlet, 'ko-', label=r'$V_a$')
            p1, = twiny1.plot(Vt1, midpointInlet, 'm*-', label=r'$V_t$')
            p2, = twiny2.plot(s1(midpointInlet), midpointInlet, 'r--', label=r'$s$')
            p3, = twiny3.plot(Pt1, midpointInlet, 'b-o', label=r'$P_t$')
            ax[0].set_ylim(self.inletSection[0].midpoint, self.inletSection[-1].midpoint)
            ax[0].set_ylabel(r'$r$')
            ax[0].set_xlabel(r'$V_a \ \frac{{m}}{{s}}$')
            twiny1.set_xlabel(r'$V_t \ \frac{{m}}{{s}}$')
            twiny2.set_xlabel(r'$s \ \frac{{J}}{{kg}}$')
            twiny3.set_xlabel(r'$P_t \ bar$')

            ax[0].legend(handles=[p0,p1,p2,p3], loc='center right')
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
            Pt2 = [self.outletSection[ii].Pt/1e+5 for ii in range(self.nSection)]

            # plotting data
            p0, = ax[1].plot(Va2, midpointOutlet, 'ko-', label=r'$V_a$')
            p1, = twiny1.plot(Vt2, midpointOutlet, 'm*-', label=r'$V_t$')
            p2, = twiny2.plot(s2(midpointOutlet), midpointOutlet, 'r--', label=r'$s$')
            p3, = twiny3.plot(Pt2, midpointOutlet, 'b-o', label=r'$P_t$')
            ax[1].set_ylim(self.outletSection[0].midpoint, self.outletSection[-1].midpoint)
            ax[1].set_ylabel(r'$r$')
            ax[1].set_xlabel(r'$V_a \ \frac{{m}}{{s}}$')
            twiny1.set_xlabel(r'$V_t \ \frac{{m}}{{s}}$')
            twiny2.set_xlabel(r'$s \ \frac{{J}}{{kg}}$')
            twiny3.set_xlabel(r'$P_t \ bar$')

            ax[1].legend(handles=[p0,p1,p2,p3], loc='center right')
            ax[1].set_title('Outlet')

            plt.tight_layout()

            # new figure 
            fig1, ax1 = plt.subplots(ncols=1, nrows=1)
            fig1.set_figwidth(10)
            fig1.set_figheight(8)
            
            # collecting new beta2 
            beta2new = [self.outletSection[ii].beta for ii in range(self.nSection)]
            # collecting beta1 
            beta1 = [self.inletSection[ii].beta for ii in range(self.nSection)]
            # collecting pressure 
            P2 = [self.outletSection[ii].P/1e+5 for ii in range(self.nSection)]

            # subdividing data 
            twiny1 = ax1.twiny()
            twiny2 = ax1.twiny() 
            twiny3 = ax1.twiny()
            # moving axis position
            twiny1.spines["top"].set_position(("axes", 1))
            twiny2.spines["top"].set_position(("axes", 1.1))
            twiny3.spines["top"].set_position(("axes", 1.2))
            
            # plotting data
            p20, = ax1.plot(Va2, midpointOutlet, 'k-*', label=r'$V_{{a \ 2}}$')
            p21, = twiny1.plot(Vt2, midpointOutlet, 'm-*', label=r'$V_{{t \ 2}}$')
            p22, = twiny2.plot(beta2start, midpointOutlet, 'r--', linewidth=2, label=r'$\beta_{{2 \ start}}$')
            p23, = twiny2.plot(beta2new, midpointOutlet, 'c', linewidth=2, label=r'$\beta_{{2 \ new}}$')
            p24, = twiny2.plot(beta1, midpointOutlet, 'g', linewidth=2, label=r'$\beta_{{1 }}$')
            p25, = twiny3.plot(P2, midpointOutlet, 'y-o', label=r'$P_{{2 }}$')
            p26, = twiny3.plot(Pt2, midpointOutlet, 'b-o', label=r'$P_{{t \ 2}}$')
            ax1.set_ylim(self.outletSection[0].midpoint, self.outletSection[-1].midpoint)
            ax1.set_ylabel(r'$r$')
            ax1.set_xlabel(r'$V_{{a \ 2}} \ \frac{{m}}{{s}}$')
            twiny1.set_xlabel(r'$V_{{t \ 2}} \ \frac{{m}}{{s}}$')
            twiny2.set_xlabel(r'$\beta_{{2 }} \ - \ \beta_{{1 }} \ ^\circ$')
            twiny3.set_xlabel(r'$P_{{2 }} \ - \ P_{{t \ 2}} \ bar$')

            ax1.legend(handles=[p20,p21,p22,p23,p24,p25,p26], loc='center left')
            
            plt.tight_layout()

            if save:
                fig0.savefig(position0)
                fig1.savefig(position1)
            else:
                plt.show()

    def generateGeometry(self, pos='/data/airfoils/naca65.txt', STLname='cad', plot=False, printout=False):
        '''
        This function generates the blade shape given already computed flow angles.
            * the geometry sections will be the midsections relative to the streamtubes. 
            * the only exception made is relative to the hub and tip streamtubes; in this case
                the section considered are no more the midsections but the tip section (tip streamtube)
                and the bottom section (hub streamtube).
        '''

        # importing libraries
        from turboClass import bladeStudy
        from geometry import bladeGenerator

        # plotting definition
        if plot:
            import random
            fig = plt.figure()
            ax = plt.axes(projection ='3d')
            number_of_colors = self.nSection + 1
            color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(number_of_colors)]

        # allocation of space for the blade geometry 
        self.blade = [] 

        for jj in range(self.nSection+1):
            # converting counter into needed counter 
            ii = jj - 1

            # this modeling approach works with beta0 >= 0 and beta1 >= 0
            # -- in order to adopt this model to each blade modeling a boolean 
            #       value for the change is generated and then used later for the 
            #       blade geoemetry generation
            checkAngle = False  

            if self.turboType == 'rotor':
                # angle allocation
                if ii == - 1: 
                    beta0 = self.inletSection[0].beta
                    beta1 = self.outletSection[0].beta
                elif ii == self.nSection:
                    beta0 = self.inletSection[-1].beta
                    beta1 = self.outletSection[-1].beta
                else: 
                    beta0 = self.inletSection[ii].beta
                    beta1 = self.outletSection[ii].beta
                # chech on angle sign for adopting Lieblein model
                if beta0 < 0:
                    beta0 = - beta0
                    beta1 = - beta1
                    checkAngle = True 

            elif self.turboType == 'stator':
                # angle allocation 
                if ii == - 1: 
                    beta0 = self.inletSection[0].alpha
                    beta1 = self.outletSection[0].alpha
                elif ii == self.nSection:
                    beta0 = self.inletSection[-1].alpha
                    beta1 = self.outletSection[-1].alpha
                else: 
                    beta0 = self.inletSection[ii].alpha
                    beta1 = self.outletSection[ii].alpha
                # chech on angle sign for adopting Lieblein model
                if beta0 < 0:
                    beta0 = - beta0
                    beta1 = - beta1
                    checkAngle = True 

            # optimal angles computation
            i, delta, theta, solidity, tbc = bladeStudy.optimalAngles(beta0, beta1, printout=printout)

            # computing optimal alpha angle
            ac = 0.5 # this is valid only for NACA-65
            alpha = bladeStudy.alphaFunc(ac, solidity, theta, tbc)

            # computing stagger angle gamma 
            gamma = beta0 - alpha

            # checking angles sign 
            if checkAngle:
                i = - i 
                delta = - delta 
                theta = - theta
                alpha = - alpha 
                gamma = - gamma 

            # computing Cl 
            Cl = ac * np.tan(np.deg2rad(theta)/4) / 0.0551515

            # pitch computation with nBlades 
            if ii == -1:
                pitch = 2 * np.pi * self.inletSection[0].bottom / self.nBlade
            elif ii == self.nSection:
                pitch = 2 * np.pi * self.inletSection[-1].tip / self.nBlade
            else: 
                pitch = 2 * np.pi * self.inletSection[ii].midpoint / self.nBlade
     
            # chord computation from solidity
            chord = pitch * solidity 

            # computing blade inclination 
            if ii == -1:
                zeta = np.rad2deg(np.arcsin((self.inletSection[0].bottom-self.outletSection[0].bottom)/chord))
            elif ii == self.nSection:
                zeta = np.rad2deg(np.arcsin((self.inletSection[-1].tip-self.outletSection[-1].tip)/chord))
            else:
                zeta = np.rad2deg(np.arcsin((self.inletSection[ii].midpoint-self.outletSection[ii].midpoint)/chord))
            
            # airfoil object generation 
            airfoil = bladeGenerator.geometryData(pos)

            # geometry fitting airfoil -> shape and chord dimensions
            airfoil.geometryFitting(Cl=Cl, chord=chord, plot=False)

            # airfoil 3D rotation
            airfoil.geometryRotation(gamma, zeta, plot=False)

            # translation of the airfoil with respect to the hub position 
            if ii == -1: 
                airfoilHub = airfoil 
                translationHub = airfoilHub.middleChord()
                translationHub[2] = 0.0
            
            # translation of profiles
            # translation vector 
            if ii == -1:
                height = self.inletSection[0].bottom - self.inletSection[0].bottom
            elif ii == self.nSection:
                height = self.inletSection[-1].tip - self.inletSection[0].bottom
            else:
                height = self.inletSection[ii].midpoint - self.inletSection[0].bottom

            # airfoil section center
            translationSection = airfoil.middleChord()
            translationSection[0] = 0 
            translationSection[1] = 0
            # translation vector  
            translationVec = translationHub + translationSection
            # translation procedure 
            airfoil.geometryTranslation(translationVec, height, plot=False)

            if plot:
                ax.plot3D(airfoil.upper[:,0], airfoil.upper[:,1], airfoil.upper[:,2], color=color[ii], label=str(ii))
                ax.plot3D(airfoil.lower[:,0], airfoil.lower[:,1], airfoil.lower[:,2], color=color[ii])

            # adding airfoil to blade object array
            self.blade.append(airfoil)

            # allocate section data 
            if ii != -1 and ii != self.nSection:
                self.inletSection[ii].allocateQuantities(i, delta, solidity, chord, pitch, gamma, Cl, tbc)
                self.outletSection[ii].allocateQuantities(i, delta, solidity, chord, pitch, gamma, Cl, tbc)

            # printout 
            if printout:
                starDim = 32
                bladeDim = np.int16((starDim - len(' BLADE ANGLES '))/2)
                print('*' * bladeDim + ' BLADE ANGLES ' + '*' * bladeDim)
                if ii == 0:
                    print('-- r             = {0:>7.3f} cm'.format(self.inletSection[ii].bottom*1e+2))
                elif ii == self.nSection - 1:
                    print('-- r             = {0:>7.3f} cm'.format(self.inletSection[ii].tip*1e+2))
                else:
                    print('-- r             = {0:>7.3f} cm'.format(self.inletSection[ii].midpoint*1e+2))
                print('-- i             = {0:>7.3f} deg'.format(i))
                print('-- delta         = {0:>7.3f} deg'.format(delta))
                print('-- alpha         = {0:>7.3f} deg'.format(alpha))
                print('-- beta0 - beta1 = {0:>7.3f} deg'.format(beta0 - beta1))
                print('-- gamma         = {0:>7.3f} deg'.format(gamma))
                print('-- theta         = {0:>7.3f} deg'.format(theta))
                print('-- zeta          = {0:>7.3f} deg'.format(zeta))
                print('-- Cl            = {0:>7.3f}'.format(np.abs(Cl)))
                print('-- tbc           = {0:>7.3f}'.format(tbc))
                print('-- s             = {0:>7.3f} cm'.format(pitch*1e+2))
                print('-- c             = {0:>7.3f} cm'.format(chord*1e+2))
                print('*' * starDim + '\n')

        if plot:
            if self.nSection < 10:
                plt.legend()
            plt.title('Blade')
            plt.show()

        # STL file generation 
        bladeGenerator.STLsaving(self.blade, STLname=STLname)

    def computeLosses(self, variableSpeed=True):
        '''
        This function computes the losses for the blade.
            inputs:
                variableSpeed   -- determines the loss coefficient with respect to Leiblein model that relates to variable axial speed for each streamtube

        '''

        # losse vector allocation
        lossVec = np.zeros(self.nSection)

        # computing losses at each section
        for ii in range(self.nSection):
            # PROFILE LOSSES
            # variables allocation -> these variables are already updated with the previous self.allocateDynamics() function call 
            # for stators V == W and alpha == beta
            # velocity 
            W1 = self.inletSection[ii].W
            W2 = self.outletSection[ii].W

            # angles
            beta1 = self.inletSection[ii].beta
            beta2 = self.outletSection[ii].beta
            
            # angles check -> the Leiblein model treats with positive angle for beta1 
            # this implies changing the angle sign with respect to beta1 sign angle  
            if beta1 < 0:
                beta1 = - beta1 
                beta2 = - beta2

            # solidity allocation    
            solidity = self.inletSection[ii].solidity

            # the Lieblein losses model dependes also on the change of axial speed along the blade section
            r1  = self.inletSection[ii].midpoint
            r2  = self.outletSection[ii].midpoint
            Vt1 = self.inletSection[ii].Vt
            Vt2 = self.outletSection[ii].Vt
            Va1 = self.inletSection[ii].Va

            # loss computation -> Lieblein model 
            if variableSpeed:
                lossVec[ii], _ = losses.lossCoeff(W1=W1, W2=W2, beta1=beta1, beta2=beta2, solidity=solidity, D=0)
            else:
                lossVec[ii], _ = losses.lossCoeff(W1=W1, W2=W2, beta1=beta1, beta2=beta2, r1=r1, r2=r2, Vt1=Vt1, Vt2=Vt2, Va1=Va1, solidity=solidity, D=0)

            # TIP LOSSES & SECONDARY LOSSES
            # computing if the rotor blade is close to the tip
            if np.abs(self.inletSection[-1].midpoint - self.inletSection[ii].midpoint) / self.inletBladeHeight < 0.2:
                endWall = True 
            elif np.abs(self.inletSection[ii].midpoint - self.inletSection[0].midpoint) / self.inletBladeHeight < 0.1:
                endWall = True
            else:
                endWall = False 
            
            # !!!!
            endWall = True
            # loss computation -> Howel additional loss model 
            if endWall:
                lossVec[ii] = lossVec[ii] + losses.lossHowell(beta1=beta1, beta2=beta2, solidity=solidity, pitch=self.inletSection[ii].pitch, bladeHeight=self.inletBladeHeight, endWall=endWall)
                
        return lossVec

    def bladeGenerator(self, Pt0, Tt0, mFlux, STLname='cad', relTolShape=1e-3, nMaxShape=100, plot=False):
        '''
        This function computes the final shape of a blade given blade number and total inlet quantites.
            inputs:
                Pt0     -- total inlet pressure 
                Tt0     -- total inlet temperature
                mFlu    -- mass flux
        '''
        
        # getting deviatoric angle vector 
        ClVecOld = np.zeros(self.nSection)
        ClVecNew = np.zeros(self.nSection)

        # loop on the geometry variation of the blade with radial equilibrium and entropy equilibrium 
        errorShape = 1 
        counterShape = 0

        # printing properties 
        starDim = 80
        geometryDim = np.int16(np.floor(starDim - len(' SHAPE LOOP   '))/2)

        while errorShape > relTolShape and counterShape < nMaxShape:
            # updating counter 
            counterShape = counterShape + 1

            # printing 
            print('-' * geometryDim + ' SHAPE LOOP {0:d} '.format(counterShape) + '-' * geometryDim)

            # blade design through iterative process on radial equilibrium 
            self.radialEquilibrium(Pt0=Pt0, Tt0=Tt0, mFlux=mFlux ,plot=plot)

            # rotor blade geometry allocation
            self.generateGeometry(pos='data/airfoils/naca65.txt', STLname=STLname, plot=False, printout=False)

            # allocating Cl vector 
            for ii in range(self.nSection):
                ClVecNew[ii] = self.inletSection[ii].Cl

            # computing relative error
            errorShape = np.abs((np.sum(ClVecNew) - np.sum(ClVecOld)) / np.max(ClVecNew))

            # reallocation Old vector to New vector 
            ClVecOld = ClVecNew 

            # printing 
            print('-- rel. error shape = {0:.4f}'.format(errorShape))

        print('-' * starDim + '\n') 

        # check flow chockin in flow passages
        self.checkChoking()

    def checkChoking(self, R=287.06, gamma=1.4):
        '''
        This function checks the choking of each section of the blade given flow properties and blade geometry.
        '''
        
        # initializing choking descriptor
        chokedFlow = False

        # initializing throat vector 
        oVec = np.zeros(self.nSection)
        rhotrVec = np.zeros(self.nSection)
        WrVec = np.zeros(self.nSection)

        # check on all the blade sections 
        for ii in range(self.nSection):
            # allocating quantities 
            rho   = self.inletSection[ii].rho 
            W     = self.inletSection[ii].W
            pitch = self.inletSection[ii].pitch 
            beta  = self.inletSection[ii].beta 
            # allocating total relative quantities
            rhotrVec[ii] = self.inletSection[ii].rhot
            WrVec[ii]    = np.sqrt(2 * gamma / (gamma + 1) * R * self.inletSection[ii].Ttr)

            # computing minumum throat 
            oVec[ii] = lieblein.oFunc(gamma=self.inletSection[ii].gamma, tbc=0.1, solidity=self.inletSection[ii].solidity, pitch=self.inletSection[ii].pitch, Cl=self.inletSection[ii].Cl)

            # chocking check condition
            if rho * W * pitch * np.cos(np.deg2rad(beta)) > rhotrVec[ii] * oVec[ii] * WrVec[ii]:
                chokedFlow = True

        # printing results
        starDim = 45
        chokingDim = np.int16(np.floor(starDim - len(' CHOKING CHECK '))/2)
        print('*' * chokingDim + ' CHOKING CHECK ' + '*' * chokingDim)
        print('-- rho1*W1*pitch*cos(beta1) < rho* W* o')
        print('-- choked                  = {0}'.format(chokedFlow))
        print('-- minimum throat position = {0:d}'.format(np.argmin(oVec)))
        print('-- minimum throat diameter = {0:>8.3f} m'.format(np.min(oVec)))
        print('-- minimum rho*            = {0:>8.3f} kg/m3'.format(np.min(rhotrVec)))
        print('-- minimum W*              = {0:>8.3f} m/s'.format(np.min(WrVec)))
        print('*' * starDim)

    def copySection(self, blade, fromSection='outlet', toSection='inlet'):
        '''
        This function copies the properties of blade another section. 
        The blades must have the same number of sections. The interpolation of data is not available (till now).
            input:
                blade       -- blade source data
                            -- turboBlade object
                fromSection -- name of the section where to copy data
                            -- inlet/outlet
                toSection   -- name of the section where to store data
                            -- inlet/outlet
        '''

        if fromSection == 'outlet':
            if toSection == 'inlet':
                print('ok')
                for ii in range(self.nSection):
                    self.inletSection[ii] = blade.outletSection[ii]
            else:
                for ii in range(self.nSection):
                    self.outletSection[ii] = blade.outletSection[ii]
        else:
            if toSection == 'inlet':
                for ii in range(self.nSection):
                    self.inletSection[ii] = blade.inletSection[ii]
            else:
                for ii in range(self.nSection):
                    self.outletSection[ii] = blade.inletSection[ii]
        
