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

        self.ID                = ID
        self.turboType         = turboType
        self.nSection          = nSection
        self.omega             = omega
        self.nBlade            = nBlade
        self.inletBladeHeight  = inletBladeHeight
        self.outletBladeHeight = outletBladeHeight
        
        # section objects allocation 
        self.inletSection = self.allocateSection(hubRadius=inletHubRadius, bladeHeight=inletBladeHeight, nSection=nSection)
        self.outletSection = self.allocateSection(hubRadius=outletHubRadius, bladeHeight=outletBladeHeight, nSection=nSection)

    def plotMeridional(self):
        '''
        This function plots the blade sections in the meridional plane.
        '''

        import matplotlib.colors as mcolors
        from matplotlib.colors   import Normalize
        from matplotlib.cm       import ScalarMappable

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

    def printMeridional(self, save=False, position0='entropyFlow.pgf', position1='betaThermo.pgf'):
        '''
        This function plots the main inlet and outlet quantities of the blade.
        '''

        # allocating marker descriptors
        markersize      = 6 
        markeredgecolor = 'k'
        markeredgewidth = 1.5   

        # setting up saving options
        if save:
            import matplotlib
            matplotlib.use("pgf")
            matplotlib.rcParams.update({
                "pgf.texsystem": "pdflatex",
                'font.family': 'serif',
                'text.usetex': True,
                'pgf.rcfonts': False,
            })

        # fig0 generation
        fig0, ax = plt.subplots(ncols=2, nrows=1)
        fig0.set_figwidth(18)
        fig0.set_figheight(9.5)

        # inlet section plot 
        twiny1 = ax[0].twiny()
        twiny2 = ax[0].twiny() 
        twiny3 = ax[0].twiny()
        twiny4 = ax[0].twiny()
        # moving axis position
        twiny1.spines["top"].set_position(("axes", 1))
        twiny2.spines["top"].set_position(("axes", 1.2))
        twiny3.spines["top"].set_position(("axes", 1.1))
        twiny4.spines["top"].set_position(("axes", 1.3))
        
        # vector declaration
        midpoint = np.zeros(self.nSection)
        Va       = np.zeros(self.nSection)
        Vt       = np.zeros(self.nSection)
        Pt       = np.zeros(self.nSection)
        s        = np.zeros(self.nSection)
        Tt       = np.zeros(self.nSection)

        # vector allocation
        for ii in range(self.nSection):
            midpoint[ii] = self.inletSection[ii].midpoint 
            Va[ii]       = self.inletSection[ii].Va 
            Vt[ii]       = self.inletSection[ii].Vt 
            Pt[ii]       = self.inletSection[ii].Pt / 1e+5
            s[ii]        = self.inletSection[ii].s 
            Tt[ii]       = self.inletSection[ii].Tt

        # plotting ax0
        ax[0].plot(Va,  midpoint, linestyle='-',  linewidth=2, marker='o', markersize=markersize, markeredgecolor=markeredgecolor, color='g', markeredgewidth=markeredgewidth, label=r'$V_a$')
        twiny1.plot(Vt, midpoint, linestyle='-',  linewidth=2, marker='s', markersize=markersize, markeredgecolor=markeredgecolor, color='m', markeredgewidth=markeredgewidth, label=r'$V_t$')
        twiny2.plot(Pt, midpoint, linestyle='--', linewidth=2, marker='^', markersize=markersize, markeredgecolor=markeredgecolor, color='b', markeredgewidth=markeredgewidth, label=r'$P_t$')
        twiny3.plot(s,  midpoint, linestyle='--', linewidth=2, marker='D', markersize=markersize, markeredgecolor=markeredgecolor, color='r', markeredgewidth=markeredgewidth, label=r'$s$')
        twiny4.plot(Tt,  midpoint, linestyle='--', linewidth=2, marker='p', markersize=markersize, markeredgecolor=markeredgecolor, color='gold', markeredgewidth=markeredgewidth, label=r'$T_t$')
        
        # ax0 setup
        ax[0].set_ylim(self.outletSection[0].midpoint, self.outletSection[-1].midpoint)
        ax[0].set_xlabel(r'$V_a \ [\frac{{m}}{{s}}]$')
        twiny1.set_xlabel(r'$V_t \ [\frac{{m}}{{s}}]$')
        twiny2.set_xlabel(r'$P_t \ [bar]$')
        twiny3.set_xlabel(r'$s \ [\frac{{J}}{{kg }}]$')
        twiny4.set_xlabel(r'$T \ [K]$')

        # setting up axes 
        ax[0].set_title('Inlet')
        ax[0].spines["left"].set_position(("axes", 1))
        ax[0].yaxis.set_ticks_position('right')

        # outlet section plot 
        twiny21 = ax[1].twiny()
        twiny22 = ax[1].twiny() 
        twiny23 = ax[1].twiny()
        twiny24 = ax[1].twiny()
        # moving axis position
        twiny21.spines["top"].set_position(("axes", 1))
        twiny22.spines["top"].set_position(("axes", 1.2))
        twiny23.spines["top"].set_position(("axes", 1.1))
        twiny24.spines["top"].set_position(("axes", 1.3))

        # vector allocation
        for ii in range(self.nSection):
            midpoint[ii] = self.outletSection[ii].midpoint 
            Va[ii]       = self.outletSection[ii].Va 
            Vt[ii]       = self.outletSection[ii].Vt 
            Pt[ii]       = self.outletSection[ii].Pt / 1e+5
            s[ii]        = self.outletSection[ii].s 
            Tt[ii]       = self.outletSection[ii].Tt

        # plotting ax1
        ax[1].plot(Va,   midpoint, linestyle='-',  marker='o', color='g', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$V_a$')
        twiny21.plot(Vt, midpoint, linestyle='-',  marker='s', color='m', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$V_t$')
        twiny22.plot(Pt, midpoint, linestyle='--', marker='^', color='b', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$P_t$')
        twiny23.plot(s,  midpoint, linestyle='--', marker='D', color='r', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$s$')
        twiny24.plot(Tt,  midpoint, linestyle='--', marker='p', color='gold', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$T_t$')
        
        # ax1 setup
        ax[1].set_title('Outlet')
        ax[1].set_ylim(self.outletSection[0].midpoint, self.outletSection[-1].midpoint)
        ax[1].set_ylabel(r'$r \ [m]$')
        ax[1].set_xlabel(r'$V_a \ [\frac{{m}}{{s}}]$')
        twiny21.set_xlabel(r'$V_t \ [\frac{{m}}{{s}}]$')
        twiny22.set_xlabel(r'$P_t \ [bar]$')
        twiny23.set_xlabel(r'$s \ [\frac{{J}}{{kg }}]$')
        twiny24.set_xlabel(r'$T \ [K]$')
    
        # getting labels information for the legend plot
        lines  = []
        labels = []
        axVec = [ax[1],twiny21, twiny22, twiny23, twiny24]
        for ax in axVec:
            Line, Label = ax.get_legend_handles_labels()
            lines.extend(Line)
            labels.extend(Label)

        # legend plot 
        fig0.legend(lines, labels, loc='lower center', ncol=len(lines))

        plt.tight_layout()

        # fig1 generation 
        fig1, ax1 = plt.subplots(ncols=2, nrows=1)
        fig1.set_figwidth(18)
        fig1.set_figheight(9.5)

        # inlet section plot 
        twiny31 = ax1[0].twiny()
        twiny32 = ax1[0].twiny()
        # moving axis position
        twiny31.spines["top"].set_position(("axes", 1))
        twiny32.spines["top"].set_position(("axes", 1.1))
        if self.turboType == 'rotor':
            # this axes is related to reaction degree that for a stator is put as 0 in the computation
            #   so it is avoided to print it in the final plot
            twiny33 = ax1[0].twiny() 
            twiny33.spines["top"].set_position(("axes", 1.2))
            twiny33.set_xlim(0,1)
        
        # vector declaration
        midpoint = np.zeros(self.nSection)
        betaIn   = np.zeros(self.nSection)
        betaOut  = np.zeros(self.nSection)
        rD       = np.zeros(self.nSection)
        gamma    = np.zeros(self.nSection)
        chord    = np.zeros(self.nSection)
        solidity = np.zeros(self.nSection)
        s        = np.zeros(self.nSection)
        Pin      = np.zeros(self.nSection)
        Pout     = np.zeros(self.nSection)
        Ptin     = np.zeros(self.nSection)
        Ptout    = np.zeros(self.nSection)
        Tin      = np.zeros(self.nSection)
        Tout     = np.zeros(self.nSection)
        Ttin     = np.zeros(self.nSection)
        Ttout    = np.zeros(self.nSection)
        MrIn     = np.zeros(self.nSection)
        MrOut    = np.zeros(self.nSection)

        # vector allocation
        for ii in range(self.nSection):
            midpoint[ii] = self.inletSection[ii].midpoint 
            gamma[ii]    = self.inletSection[ii].gamma
            rD[ii]       = self.inletSection[ii].rD
            chord[ii]    = self.inletSection[ii].chord 
            solidity[ii] = self.inletSection[ii].solidity
            s[ii]        = self.inletSection[ii].s 
            Pin[ii]      = self.inletSection[ii].P / 1e+5
            Pout[ii]     = self.outletSection[ii].P / 1e+5
            Ptin[ii]     = self.inletSection[ii].Pt / 1e+5
            Ptout[ii]    = self.outletSection[ii].Pt / 1e+5
            Tin[ii]      = self.inletSection[ii].T 
            Tout[ii]     = self.outletSection[ii].T 
            Ttin[ii]     = self.inletSection[ii].Tt 
            Ttout[ii]    = self.outletSection[ii].Tt 
            if self.turboType == 'stator':
                MrIn[ii]     = self.inletSection[ii].M
                MrOut[ii]    = self.outletSection[ii].M
                betaIn[ii]   = self.inletSection[ii].alpha
                betaOut[ii]  = self.outletSection[ii].alpha 
            elif self.turboType == 'rotor':
                MrIn[ii]     = self.inletSection[ii].Mr
                MrOut[ii]    = self.outletSection[ii].Mr
                betaIn[ii]   = self.inletSection[ii].beta
                betaOut[ii]  = self.outletSection[ii].beta 

        # plotting ax0
        p30, = ax1[0].plot(betaIn,    midpoint, linestyle='-',  linewidth=2, marker='o', markersize=markersize, markeredgecolor=markeredgecolor, color='g', markeredgewidth=markeredgewidth, label=r'$\beta_{{1 }}$')
        p31, = ax1[0].plot(betaOut,   midpoint, linestyle='-',  linewidth=2, marker='s', markersize=markersize, markeredgecolor=markeredgecolor, color='m', markeredgewidth=markeredgewidth, label=r'$\beta_{{2 }}$')
        p32, = ax1[0].plot(gamma,     midpoint, linestyle='-',  linewidth=2, marker='p', markersize=markersize, markeredgecolor=markeredgecolor, color='c', markeredgewidth=markeredgewidth, label=r'$\gamma$')
        p33, = twiny31.plot(chord,    midpoint, linestyle='--', linewidth=2, marker='d', markersize=markersize, markeredgecolor=markeredgecolor, color='y', markeredgewidth=markeredgewidth, label=r'$c$')
        p34, = twiny32.plot(solidity, midpoint, linestyle='--', linewidth=2, marker='D', markersize=markersize, markeredgecolor=markeredgecolor, color='r', markeredgewidth=markeredgewidth, label=r'$\sigma$')
        if self.turboType == 'rotor':
            p35, = twiny33.plot(rD,       midpoint, linestyle='--', linewidth=2, marker='^', markersize=markersize, markeredgecolor=markeredgecolor, color='b', markeredgewidth=markeredgewidth, label=r'$\chi$')
        
        # ax0 setup
        ax1[0].set_ylim(self.outletSection[0].midpoint, self.outletSection[-1].midpoint)
        ax1[0].set_xlabel(r'$\beta_{{1 }} \ \beta_{{2 }} \ \gamma \ [^{\circ}]$')
        twiny31.set_xlabel(r'$c \ [m]$')
        twiny32.set_xlabel(r'$\sigma$')
        if self.turboType == 'rotor':
            twiny33.set_xlabel(r'$\chi$')

        # setting up axes 
        ax1[0].set_title('Main geometric/load quantities')
        ax1[0].spines["left"].set_position(("axes", 1))
        ax1[0].yaxis.set_ticks_position('right')

        # legend generation 
        try:
            ax1[0].legend(handles=[p30,p31,p32,p33,p34,p35], loc='upper right', bbox_to_anchor=[0,1])
        except:
            ax1[0].legend(handles=[p30,p31,p32,p33,p34], loc='upper right', bbox_to_anchor=[0,1])
        
        plt.tight_layout()

        # outlet section plot 
        twiny41 = ax1[1].twiny()
        twiny42 = ax1[1].twiny()
        # moving axis position
        twiny41.spines["top"].set_position(("axes", 1))
        twiny42.spines["top"].set_position(("axes", 1.1))

        # plotting ax1
        p40, = ax1[1].plot(Pin,    midpoint, linestyle='-',        marker='o', color='peru', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$P_{{in }}$')
        p41, = ax1[1].plot(Pout,   midpoint, linestyle='-',        marker='s', color='deeppink', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$P_{{out }}$')
        p42, = twiny41.plot(Tin,   midpoint, linestyle='--',       marker='^', color='lightseagreen', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$T_{{in }}$')
        p43, = twiny41.plot(Tout,  midpoint, linestyle='--',       marker='D', color='indianred', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$T_{{out }}$')
        p44, = ax1[1].plot(Ptin,   midpoint, linestyle='-',        marker='o', color='skyblue', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$P_{t_{in }}$')
        p45, = ax1[1].plot(Ptout,  midpoint, linestyle='-',        marker='s', color='lime', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$P_{t_{out }}$')
        p46, = twiny41.plot(Ttin,  midpoint, linestyle='--',       marker='^', color='salmon', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$T_{t_{in }}$')
        p47, = twiny41.plot(Ttout, midpoint, linestyle='--',       marker='D', color='royalblue', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$T_{t_{out }}$')
        p48, = twiny42.plot(MrIn,  midpoint, linestyle='dashdot',  marker='8', color='red', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$M_{r_{in }}$')
        p49, = twiny42.plot(MrOut, midpoint, linestyle='dashdot',  marker='8', color='gold', markersize=markersize, markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor, label=r'$M_{r_{out }}$')

        # ax1 setup
        ax1[1].set_title('Main thermodynamic quantities')
        ax1[1].set_ylim(self.outletSection[0].midpoint, self.outletSection[-1].midpoint)
        ax1[1].set_ylabel(r'$r \ [m]$')
        ax1[1].set_xlabel(r'$P \ [bar]$')
        twiny41.set_xlabel(r'$T \ [K]$')
        twiny42.set_xlabel(r'$M$')
        twiny42.set_xlim(0,1.1)

        # legend generation 
        ax1[1].legend(handles=[p40,p41,p42,p43,p44,p45,p46,p47,p48,p49], loc='upper left', bbox_to_anchor=[1,1])

        plt.tight_layout()

        # setting up figure supertitle 
        fig1.suptitle(self.turboType)

        if save:
            fig0.savefig(position0, bbox_inches='tight')
            fig1.savefig(position1, bbox_inches='tight')
        else:
            plt.show()

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
        height     = bladeHeight/nSection

        # loop generation for the section generation
        for ii in range(nSection):
            if hubRadius != 0:
                midpoint = hubRadius + ii * height + height / 2 
                bottom   = hubRadius + ii * height
                tip      = hubRadius + (ii+1) * height
                pitch    = 2 * np.pi * midpoint / self.nBlade
            elif tipRadius != 0:
                midpoint = tipRadius - ii * height - height/2
                bottom   = tipRadius - (ii+1) * height
                tip      = tipRadius - ii * height 
                pitch    = 2 * np.pi * midpoint / self.nBlade

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

    def allocateKinetics(self, rMean=0, VtMean=0, VaMean=0, omega=0, Vt1=0, Vt2=0, rMean1=0, rMean2=0, section='inlet', kind='FV', a=0, b=0, n=0, func=0):
        '''
        This function allocates the velocity vectors at each section points using kind model.
            inputs:
                rMean   -- mean radius
                VtMean  -- tangential mean velocity 
                VaMean  -- axial mean velocity 
                omega   -- angular velocity
                section -- inlet/outlet section 
                kind    -- vortex model 
                        -- FV  => free vortex model 
                        -- VVD => variable vortex design
                        -- MVD => mixed vorted design
                        -- C   => constant angle model 

        '''

        if kind == 'FV':
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
        elif kind == 'VVD':
            # variables computation with respect to n 
            a = (Vt2 * rMean2 + Vt1 * rMean1) / (rMean1**(n+1) + rMean2**(n+1))
            b = a * rMean1**(n+1) - Vt1 * rMean1

            if section == 'inlet':
                for ii in range(self.nSection):
                    # computing variable vortex design inlet velocity
                    r = self.inletSection[ii].midpoint 
                    Vt = a * r**n - b / r    
                    # rotation speed 
                    U = omega * r 

                    # data allocation
                    self.inletSection[ii].allocateKinetics(VaMean, Vt, U)    
            elif section == 'outlet':
                for ii in range(self.nSection):
                    # computing variable vortex design inlet velocity
                    r = self.outletSection[ii].midpoint 
                    Vt = a * r**n + b / r    
                    # rotation speed 
                    U = omega * r 

                    # data allocation
                    self.outletSection[ii].allocateKinetics(VaMean, Vt, U)
        elif kind == 'MVD':
            if section == 'inlet':
                # variables computation 
                a = (Vt1 - b * rMean1) * rMean1
                
                for ii in range(self.nSection):
                    # computing mixed vortex design inlet velocity
                    r = self.inletSection[ii].midpoint
                    Vt = a / r + b * r

                    # rotation speed 
                    U = omega * r

                    # data allocation 
                    self.inletSection[ii].allocateKinetics(VaMean, Vt, U)
            elif section == 'outlet':
                # variables computation 
                a = (Vt2 - b * rMean2) * rMean2

                for ii in range(self.nSection):
                    # computing mixed vortex design outlet velocity 
                    r = self.outletSection[ii].midpoint
                    Vt = a / r + b * r 
                    
                    # rotation speed 
                    U = omega * r 

                    # data allocation 
                    self.outletSection[ii].allocateKinetics(VaMean, Vt, U)
        elif kind == 'eqn':
            if section == 'inlet':
                for ii in range(self.nSection):
                    # section position 
                    r = self.inletSection[ii].midpoint
                    # Vt computation 
                    Vt = func(r)

                    # rotation speed 
                    U = omega * r

                    # data allocation
                    self.inletSection[ii].allocateKinetics(VaMean, Vt, U)
            elif section == 'outlet':
                for ii in range(self.nSection):
                    # section position 
                    r = self.outletSection[ii].midpoint
                    # Vt computation 
                    Vt = func(r)

                    # rotation speed 
                    U = omega * r

                    # data allocation
                    self.outletSection[ii].allocateKinetics(VaMean, Vt, U)
        else:
            raise ValueError("Invalid vortex model")

        #if section == 'inlet':
        #    print('Inlet section')
        #    for ii in range(self.nSection):
        #        print('Section = ', ii)
        #        print('alpha   = ', self.inletSection[ii].beta)
        #        print('Vt      = ', self.inletSection[ii].Vt)
        #else:
        #    print('Outlet section')
        #    for ii in range(self.nSection):
        #        print('Section = ', ii)
        #        print('alpha   = ', self.outletSection[ii].beta)
        #        print('Vt      = ', self.outletSection[ii].Vt)

    def allocateThermodynamics(self, Tt0=0, Pt0=0, eta=1, R=287.06, gamma=1.4):
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
            # check if the thermodynamics properites are not already computed for the inlet section
            if Tt0 == 0 and Pt0 == 0:
                Tt0 = self.inletSection[ii].Tt
                Pt0 = self.inletSection[ii].Pt

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
                Leu = 0.0
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

    def allocateShape(self, bladeHeight, AR, nBlade):
        '''
        This function allocates the section geometry. The allocation comes from a first analysis of the profile + secondary losses on the mean line.
            inputs:
                bladeHeight -- inlet blade height
                AR          -- aspect ratio -> bladeHeight / chord
                nBlade      -- blade number 
        '''

        # blade number allocation
        self.nBlade = nBlade

        # computing chord, pitch and solidity
        try:
            for ii in range(self.nSection):
                chord    = bladeHeight / AR[ii]
                pitch    = self.inletSection[ii].midpoint * 2 * np.pi / nBlade
                solidity = chord / pitch
                # allocate data  
                self.inletSection[ii].chord     = chord 
                self.inletSection[ii].pitch     = pitch 
                self.inletSection[ii].solidity  = solidity
                self.outletSection[ii].chord    = chord 
                self.outletSection[ii].pitch    = pitch 
                self.outletSection[ii].solidity = solidity
        except:
            for ii in range(self.nSection):
                chord    = bladeHeight / AR
                pitch    = self.inletSection[ii].midpoint * 2 * np.pi / nBlade
                solidity = chord / pitch
                # allocate data  
                self.inletSection[ii].chord     = chord 
                self.inletSection[ii].pitch     = pitch 
                self.inletSection[ii].solidity  = solidity
                self.outletSection[ii].chord    = chord 
                self.outletSection[ii].pitch    = pitch 
                self.outletSection[ii].solidity = solidity

    def radialEquilibrium(self, mFlux, clearance, nMaxS=100, nMaxFlux=100, tolS=0.2, tolFlux=1e-2, NISRE=True, plot=False, save=False, position0='entropyFlow.pgf', position1='betaThermo.pgf', R=287.06, gamma=1.4):
        '''
        This function computes the radial equilibrium of the section taking into account losses. 
            inputs:
                mFlux       -- mass flux 
                clearance   -- rotor tip clearance
                nMaxS       -- # of entropy loop iterations
                nMaxFlux    -- # of continuity loop iterations
                tolS        -- entropy loop relative tolerance
                tolFlux     -- continuity loop relative tolerance 
                NISRE       -- boolean value that enables pressure losses study
                plot        -- boolean value for the plotting of the results
                save        -- boolean value for the saving of the results
                position0   -- saving path for the first figure
                position1   -- saving path for the last figure 
                R           -- gas constant 
                gamma       -- specific heat ratio 
            
            function steps:
                1. setting up variables
                    1.1. interpolating from blade inlet/outlet allocated variables the needed function that will be used in the radial equilibrium 
                2. setting up loop properties 
                    2.1. computing entropy from geometry and flow data
                    2.2. continuity loop
                        2.2.1. setting up radial equilibrium ODE (radialFunc) that will not change during the inner loop because the 
                               outlet entropy envelop is defined at the beginning of the continuity loop
                        2.2.2. setting up boundary conditions: y0 == Va**2
                               setting up study points: t == r 
                        2.2.3. computing losses -> this vector changes at all the steps because it is dependant of the flow angles that 
                               change with Va (output of the ODE integration)
                        2.2.4. solving ODE problem and getting Va**2 envelop 
                        2.2.5. computing Va from Va**2
                        2.2.6. kinematics allocation -> dependant on Va
                        2.2.7. thermodynamics allocation -> dependant on Va 
                        2.2.8. computing mass flux 
                        2.2.9. computing flux error 
                    2.3. reallocating new entropy data 
                    2.4. computing new flow thermodynamic properties 
        '''         
    
        # cP computation
        cP = gamma / (gamma - 1) * R

        # omega 
        omega = self.omega 

        # INLET VARIABLES INTERPOLATION 
        # these variables do not change so they are set once in all the process
        # midpoint --> blade inlet  
        midpointInlet = [self.inletSection[ii].midpoint for ii in range(self.nSection)]

        if self.turboType == 'stator':
            for ii in range(self.nSection):
                self.outletSection[ii].Pt = self.inletSection[ii].Pt
                self.outletSection[ii].Tt = self.inletSection[ii].Tt
                #lossVec = self.computeLosses(mFlux=mFlux, clearance=clearance)

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
        relErrorS = 1.0

        # initialize print
        lineLenght = 80
        iterativeLenght = np.int16((lineLenght - len(' RADIAL EQ. '))/2)
        print('' + '*' * iterativeLenght + ' RADIAL EQ. ' + '*' * iterativeLenght)

        #if plot:
        #    fig = plt.figure()
        #    plt.plot(s1(midpointInlet),   midpointInlet, 'r')
        #    plt.plot(Tt1(midpointInlet),  midpointInlet, 'b')
        #    plt.plot(rVt1(midpointInlet), midpointInlet, 'g')
        #    plt.plot(Va1(midpointInlet),  midpointInlet, 'black')
        #    plt.plot(Vt2(midpointOutlet), midpointInlet, 'royalblue')
        #    plt.show()

        # entropy outer loop 
        while  counterS < nMaxS and relErrorS > tolS: 
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
                dydt = 2 * (y / (2 * cP) * ds2 - Tt1(t) * ds2 - omega * rVt2(t) / cP * ds2  + omega * rVt1(t) / cP * ds2 + Vt2(t)**2 / (2 * cP) * ds2 - Vt2(t) / t * drVt2 +  dTt1 + omega * drVt2 - omega * drVt1)

                return dydt 
            
            # radius vector allocation 
            t = midpointOutlet 

            def radialFunc1(y, t):
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
                dydt = 2 * (y / (2 * cP) * ds2 - Tt1(t) * ds2 - omega * rVt2(t) / cP * ds2  + omega * rVt1(t) / cP * ds2 + Vt2(t)**2 / (2 * cP) * ds2 - Vt2(t) / t * drVt2 + cP * dTt1 + omega * drVt2 - omega * drVt1)
                
                Tt1ds2 = Tt1(t) * ds2
                Vt22ds2 = Vt2(t)**2 / (2 * cP) * ds2
                Vt2rdrVt2 = Vt2(t) / t * drVt2
                cpdTt1 = cP * dTt1

                return dydt/2, Tt1(t), ds2, rVt2(t), rVt1(t), Vt2(t), drVt2, drVt1, -Tt1ds2, Vt22ds2, -Vt2rdrVt2, cpdTt1

            #if plot:
            #    dydtval = np.zeros(self.nSection) 
            #    Tt1val = np.zeros(self.nSection)
            #    ds2val = np.zeros(self.nSection) 
            #    rVt2val = np.zeros(self.nSection) 
            #    rVt1val = np.zeros(self.nSection) 
            #    Vt2val = np.zeros(self.nSection) 
            #    drVt2val = np.zeros(self.nSection) 
            #    drVt1val = np.zeros(self.nSection)
            #    Tt1ds2val = np.zeros(self.nSection)
            #    Vt22ds2val = np.zeros(self.nSection)
            #    Vt2rdrVt2val = np.zeros(self.nSection)
            #    cpdTt1val = np.zeros(self.nSection)
            #    for ii in range(self.nSection):
            #        dydtval[ii], Tt1val[ii], ds2val[ii], rVt2val[ii], rVt1val[ii], Vt2val[ii], drVt2val[ii], drVt1val[ii], Tt1ds2val[ii], Vt22ds2val[ii], Vt2rdrVt2val[ii], cpdTt1val[ii] = radialFunc1(0, t[ii])
            #
            #    fig2 = plt.figure()
            #    plt.plot(dydtval, t, 'r', label='dydt')
            #    plt.plot(Tt1ds2val, t, 'b', label='Vt2 ds2/sr')
            #    plt.plot(Vt2rdrVt2val, t, 'g', label='Vt2/r d r Vt2 / dr')
            #    plt.plot(cpdTt1val, t, 'k', label='cp dTt1 / dr')
            #    fig3 = plt.figure()
            #    plt.plot(Tt1val, t, 'b', label='Tt1')
            #    plt.plot(ds2val, t, 'k', label='ds2/dr')
            #    plt.plot(rVt2val, t, 'gold', label='r Vt2')
            #    plt.plot(rVt1val, t, 'green', label='r Vt1')
            #    plt.plot(Vt2val, t, 'magenta', label='Vt2')
            #    plt.plot(drVt2val, t, 'pink', label = 'd r Vt2 / dr')
            #    plt.plot(drVt1val, t, 'slategray', label = 'd r Vt1 / dr')
            #    fig3.legend()
            #    fig2.legend()
            #    plt.show()

            # setting up boundary condition
            y0 = self.outletSection[0].Va**2

            # setting up tolerances and counters 
            # iterations counter
            counterFlux = 0
            # mass flux 
            relErrorFlux = 1.0

            # computing losses for the new loop cycle (used if relError > tol) 
            # LOSSES COMPUTATION => target s2
            # computing pressure losses 
            # NISRE allows studying the losses
            if NISRE:
                lossVec = self.computeLosses(mFlux=mFlux, clearance=clearance)
            else:
                lossVec = np.zeros(self.nSection)

            # storing s2 old 
            s2old = [self.outletSection[ii].s for ii in range(self.nSection)]

            # computing entropy generation through losses 
            for ii in range(self.nSection):    
                # COMPUTING NEW s VALUES
                # new outlet pressure computation -> using losses
                if self.turboType == 'rotor':
                    self.outletSection[ii].Ptr = self.inletSection[ii].Ptr - lossVec[ii] * (self.inletSection[ii].Ptr - self.inletSection[ii].P)
                    # entropy computation
                    self.outletSection[ii].s = self.inletSection[ii].s - R * np.log(self.outletSection[ii].Ptr / self.inletSection[ii].Ptr) 
                elif self.turboType == 'stator':
                    self.outletSection[ii].Pt = self.inletSection[ii].Pt - lossVec[ii] * (self.inletSection[ii].Pt - self.inletSection[ii].P)
                    # entropy computation
                    self.outletSection[ii].s = self.inletSection[ii].s - R * np.log(self.outletSection[ii].Pt / self.inletSection[ii].Pt) 

            # storing new s2 
            s2new = [self.outletSection[ii].s for ii in range(self.nSection)]

            # loss error computation 
            relErrorS = 0.0
            for ii in range(self.nSection):
                if np.abs(s2old[ii] - s2new[ii]) > relErrorS:
                    relErrorS = np.abs(s2old[ii] - s2new[ii])/np.max(s2new)
        
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

                # computing all the new thermodynamic quantities after the radial equilibrium is satisfied
                # computing change in total pressure due to losses 
                for ii in range(self.nSection): 
                    # static temperature computation
                    self.outletSection[ii].T = self.outletSection[ii].Tt - self.outletSection[ii].V**2 / (2 * cP)

                    # rotor and stator work with different total pressures
                    # rotor -> total relative pressure
                    # stator -> total pressure 
                    if self.turboType == 'rotor':
                        # total relative temperature computation
                        #print(self.outletSection[ii].V**2 - self.outletSection[ii].W**2)
                        self.outletSection[ii].Ttr = self.outletSection[ii].Tt - (self.outletSection[ii].V**2 - self.outletSection[ii].W**2) / (2 * cP)
                        # total pressure computation 
                        self.outletSection[ii].Pt = self.outletSection[ii].Ptr * (self.outletSection[ii].Tt/self.outletSection[ii].Ttr)**(gamma/(gamma-1))
                        # static pressure computation
                        self.outletSection[ii].P = self.outletSection[ii].Pt * (self.outletSection[ii].T/self.outletSection[ii].Tt)**(gamma/(gamma-1))
                        # static density computation
                        self.outletSection[ii].rho = self.outletSection[ii].P / (R * self.outletSection[ii].T)
                        # total density computation
                        self.outletSection[ii].rhot = self.outletSection[ii].Pt / (R * self.outletSection[ii].Tt)
                        # total relative density computation 
                        self.outletSection[ii].rhotr = self.outletSection[ii].Ptr / (R * self.outletSection[ii].Ttr)
                        # sound speed computation 
                        self.outletSection[ii].a = np.sqrt(gamma * R * self.outletSection[ii].T)
                        # mach number computation
                        self.outletSection[ii].M = self.outletSection[ii].V / self.outletSection[ii].a
                        # relative mach number computation
                        self.outletSection[ii].Mr = self.outletSection[ii].W / self.outletSection[ii].a
                    elif self.turboType == 'stator':
                        # static pressure computation
                        self.outletSection[ii].P = self.outletSection[ii].Pt * (self.outletSection[ii].T/self.outletSection[ii].Tt)**(gamma/(gamma-1))
                        # static density computation
                        self.outletSection[ii].rho = self.outletSection[ii].P / (R * self.outletSection[ii].T)
                        # total density computation
                        self.outletSection[ii].rhot = self.outletSection[ii].Pt / (R * self.outletSection[ii].Tt)
                        # sound speed computation 
                        self.outletSection[ii].a = np.sqrt(gamma * R * self.outletSection[ii].T)
                        # mach number computation
                        self.outletSection[ii].M = self.outletSection[ii].V / self.outletSection[ii].a

                        # printing

                # check mass flux 
                newFlux = 0
                for ii in range(self.nSection):
                    newFlux = newFlux + self.outletSection[ii].mFlux()

                # relative mass flux error
                relErrorFlux = np.abs(newFlux - mFlux) / mFlux

                # setting up initial velocity for the new ODE 
                if newFlux > mFlux:
                    if relErrorFlux < 0.25:
                        y0 = y0 * (1 - relErrorFlux)
                    else:
                        y0 = y0 * (1 - 0.25)
                else: 
                    if relErrorFlux < 0.25:
                        y0 = y0 * (1 + relErrorFlux)
                    else:
                        y0 = y0 * (1 + 0.25)

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
    
                # reaction degree computation 
                if self.turboType == 'rotor':
                    reactionDegree = (self.inletSection[ii].T - self.outletSection[ii].T) / (self.inletSection[ii].Tt - self.outletSection[ii].Tt)
                    self.inletSection[ii].rD = reactionDegree
                    self.outletSection[ii].rD = reactionDegree
            
        # plotting interpolated functions 
        if plot:
            self.printMeridional(save=save, position0=position0, position1=position1)

        return lossVec

    def generateGeometry(self, pos='/data/airfoils/naca65.txt', STLname='cad', plot=False, printout=False):
        '''
        This function generates the blade shape given already computed flow angles.
            * the geometry sections will be the midsections relative to the streamtubes. 
            * the only exception made is relative to the hub and tip streamtubes; in this case
                the section considered are no more the midsections but the tip section (tip streamtube)
                and the bottom section (hub streamtube).
        '''

        # importing libraries
        from geometry   import bladeGenerator
        from turboClass import bladeStudy

        # plotting definition
        if plot:
            import random
            ax = plt.axes(projection ='3d')
            number_of_colors = self.nSection + 1
            color = ["#"+''.join([random.choice('0123456789ABCDEF') for _ in range(6)]) for _ in range(number_of_colors)]

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
                    beta0    = self.inletSection[0].beta
                    beta1    = self.outletSection[0].beta
                    solidity = self.inletSection[0].solidity
                    tbc      = self.inletSection[0].tbc
                elif ii == self.nSection:
                    beta0    = self.inletSection[-1].beta
                    beta1    = self.outletSection[-1].beta
                    solidity = self.inletSection[-1].solidity
                    tbc      = self.inletSection[-1].tbc
                else: 
                    beta0    = self.inletSection[ii].beta
                    beta1    = self.outletSection[ii].beta
                    solidity = self.inletSection[ii].solidity
                    tbc      = self.inletSection[ii].tbc
                # chech on angle sign for adopting Lieblein model
                if beta0 < 0:
                    beta0 = - beta0
                    beta1 = - beta1
                    checkAngle = True 

            elif self.turboType == 'stator':
                # angle allocation 
                if ii == - 1: 
                    beta0    = self.inletSection[0].alpha
                    beta1    = self.outletSection[0].alpha
                    solidity = self.inletSection[0].solidity
                    tbc      = self.inletSection[0].tbc
                elif ii == self.nSection:
                    beta0    = self.inletSection[-1].alpha
                    beta1    = self.outletSection[-1].alpha
                    solidity = self.inletSection[-1].solidity
                    tbc      = self.inletSection[-1].tbc
                else: 
                    beta0    = self.inletSection[ii].alpha
                    beta1    = self.outletSection[ii].alpha
                    solidity = self.inletSection[ii].solidity
                    tbc      = self.inletSection[ii].tbc
                # chech on angle sign for adopting Lieblein model
                if beta0 < 0:
                    beta0 = - beta0
                    beta1 = - beta1
                    checkAngle = True 

            # optimal angles computation
            i, delta, theta = bladeStudy.optimalAngles(beta0, beta1, solidity, tbc, printout=printout)

            # computing optimal alpha angle
            ac    = 0.5 # this is valid only for NACA-65 -> different values of ac need different bladeStudy.alphaFunc()
            alpha = bladeStudy.alphaFunc(ac, solidity, theta, tbc)

            # computing stagger angle gamma 
            gamma = beta0 - alpha

            # checking angles sign 
            if checkAngle:
                i     = - i 
                delta = - delta 
                theta = - theta
                alpha = - alpha 
                gamma = - gamma 

            # allocate main angles
            try:
                self.inletSection[ii].theta  = theta
                self.inletSection[ii].i      = i
                self.outletSection[ii].theta = theta 
                self.outletSection[ii].i     = i 
            except:
                pass

            # computing Cl 
            Cl = ac * np.tan(np.deg2rad(theta)/4) / 0.0551515

            # pitch computation with # of blades
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
                airfoilHub        = airfoil 
                translationHub    = airfoilHub.middleChord()
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
            translationSection    = airfoil.middleChord()
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
        bladeGenerator.STLsaving(self.blade, STLname=STLname, kind=self.turboType)

    def computeLosses(self, mFlux, clearance=3e-3, variableSpeed=False):
        '''
        This function computes the losses for the blade. The blade losses are referred to the Leiblein and Howell loss model.
            inputs:
                clearance       -- rotor tip clearance
                variableSpeed   -- defines the loss coefficient with respect to Leiblein model that relates to variable axial speed for each streamtube
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
            theta = np.abs(self.inletSection[ii].theta)
            i     = np.abs(self.inletSection[ii].i)
            tbc   = self.inletSection[ii].tbc
            c     = self.inletSection[ii].chord
            gamma = self.inletSection[ii].gamma
            pitch = self.inletSection[ii].pitch 
            
            if self.turboType == 'rotor':
                M1    = self.inletSection[ii].Mr
                Ttr   = self.inletSection[ii].Ttr
                Ptr   = self.inletSection[ii].Ptr
                P     = self.inletSection[ii].P
            else:
                M1    = self.inletSection[ii].M 
                Ttr   = self.inletSection[ii].Tt
                Ptr   = self.inletSection[ii].Pt
                P     = self.inletSection[ii].P

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
                lossVec[ii], _ = losses.profileLosses(W1=W1, W2=W2, beta1=beta1, beta2=beta2, solidity=solidity, D=0)
            else:
                lossVec[ii], _ = losses.profileLosses(W1=W1, W2=W2, beta1=beta1, beta2=beta2, r1=r1, r2=r2, Vt1=Vt1, Vt2=Vt2, Va1=Va1, solidity=solidity, D=0)

            # compressibility effects
            if theta != 0 and i != 0:
                lossVec[ii] = losses.machLosses(beta1, beta2, theta, i, W1, M1, Ttr, lossVec[ii], solidity)
            
            # secondary flow losses
            lossVec[ii] = lossVec[ii] + losses.lossHowell(beta1=beta1, beta2=beta2, solidity=solidity, pitch=pitch, bladeHeight=self.inletBladeHeight, endWall=False)

            # end wall losses computation
            lossVec[ii] = lossVec[ii] + losses.lossHowell(beta1=beta1, beta2=beta2, solidity=solidity, pitch=pitch, bladeHeight=self.inletBladeHeight, endWall=True)

            # shock losses computation 
            if M1 >= 1: 
                lossVec[ii] = lossVec[ii] + losses.shockLosses(theta, tbc, c, gamma, pitch, W1, M1, Ptr, P)

            # TIP LOSSES 
            # computing tip losses 
            if self.turboType == 'rotor':
                # data allocation
                mid       = np.int16(self.nSection/2)
                r         = self.inletSection[ii].midpoint
                VaInMid   = self.inletSection[mid].Va
                VaOutMid  = self.outletSection[mid].Va
                VtInMid   = self.inletSection[mid].Vt
                VtOutMid  = self.outletSection[mid].Vt
                rhoInMid  = self.inletSection[mid].rho
                rhoOutMid = self.outletSection[mid].rho
                rhotMid   = self.inletSection[mid].rhot
                r1Mid     = self.inletSection[mid].midpoint
                r2Mid     = self.outletSection[mid].midpoint
                rHub      = self.inletSection[0].bottom
                rTip      = self.inletSection[-1].tip
                nBlade    = self.nBlade 
                chordMid  = self.inletSection[mid].chord
                gammaMid  = self.inletSection[mid].gamma
                Nrow      = self.ID
                deltaPiso = self.inletSection[ii].Pt - self.inletSection[ii].P

                # tip loss computation
                lossVec[ii] = lossVec[ii] + losses.tipLosses(r, clearance, VaInMid, VaOutMid, VtInMid, VtOutMid, rhoInMid, rhoOutMid, rhotMid, r1Mid, r2Mid, rHub, rTip, nBlade, chordMid, gammaMid, Nrow, mFlux, deltaPiso)

        return lossVec

    def bladeGenerator(self, mFlux, clearance=3e-3, NISRE=True, STLname='cad', relTolShape=1e-3, nMaxShape=100, nMaxFlux=100, nMaxS=10, plot=False):
        '''
        This function computes the final shape of a blade given blade number and total inlet quantites.
            inputs:
                mFlu        -- mass flux
                clearance   -- rotor tip clearance 

            function steps:
                1. setting up loop tolerances and storing variables for the error check
                2. external loop based on the change in blade geoometry
                    2.1. internal loop based on continuity and entropy production
                    2.2. blade generation
                    2.3. error computation
                3. check results
                4. check choking 
        '''
        
        # getting deviatoric angle vector 
        ClVecOld = np.zeros(self.nSection)
        ClVecNew = np.zeros(self.nSection)

        # loop on the geometry variation of the blade with radial equilibrium and entropy equilibrium 
        errorShape = 1 
        counterShape = 0

        # printing properties 
        starDim = 80
        geometryDim = np.int16(np.floor(starDim - len(' SHAPE ITERATION   '))/2)

        while errorShape > relTolShape and counterShape < nMaxShape:
            # updating counter 
            counterShape = counterShape + 1

            # printing 
            print('-' * geometryDim + ' SHAPE ITERATION {0:d} '.format(counterShape) + '-' * geometryDim)

            # blade design through iterative process on radial equilibrium 
            lossVec = self.radialEquilibrium(mFlux=mFlux, clearance=clearance, NISRE=NISRE, plot=plot, nMaxFlux=nMaxFlux, nMaxS=nMaxS)

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
            #print('-- rel. error shape = {0:.4f}'.format(errorShape))

        print('-' * starDim + '\n') 

        # check flow chockin in flow passages
        self.checkChoking()

        # computing mean total pressure
        self.computeMeanPressure(kind='total')

        return lossVec

    def checkChoking(self, R=287.06, gamma=1.4):
        '''
        This function checks the choking of each section of the blade given flow properties and blade geometry.
            inputs:
                R       -- gas constant
                gamma   -- specific heat ratio
        '''
        
        # initializing choking descriptor
        chokedFlow = False

        # initializing throat vector 
        oVec     = np.zeros(self.nSection)
        rhotrVec = np.zeros(self.nSection)
        WrVec    = np.zeros(self.nSection)

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
            gamma_    = self.inletSection[ii].gamma     # pay attention with the specific heat ratio gamma
            tbc_      = self.inletSection[ii].tbc
            solidity_ = self.inletSection[ii].solidity
            pitch_    = self.inletSection[ii].pitch
            Cl_       = self.inletSection[ii].Cl
            # minimum throath
            oVec[ii] = lieblein.oFunc(gamma=gamma_, tbc=tbc_, solidity=solidity_, pitch=pitch_, Cl=Cl_)

            # chocking check condition
            if rho * W * pitch * np.cos(np.deg2rad(beta)) > rhotrVec[ii] * oVec[ii] * WrVec[ii]:
                chokedFlow = True

        # printing results
        starDim = 82
        chokingDim = np.int16(np.floor(starDim - len(' CHOKING CHECK '))/2)
        print('*' * chokingDim + ' CHOKING CHECK ' + '*' * chokingDim)
        print('-- ASME DEFINITION: rho1*W1*pitch*cos(beta1) < rho* W* o')
        print('-- choked                  = {0}'.format(chokedFlow))
        print('-- minimum throat position = {0:>4d}           -- minimum throat diam = {1:>8.3f} cm'.format(np.argmin(oVec), np.min(oVec)*1e+2))
        print('-- minimum rho*            = {0:>8.3f} kg/m3 -- minimum W*          = {1:>8.3f} m/s'.format(np.min(rhotrVec), np.min(WrVec)))
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
                for ii in range(self.nSection):
                    # kinetics allocation 
                    Va    = blade.outletSection[ii].Va
                    Vt    = blade.outletSection[ii].Vt
                    V     = blade.outletSection[ii].V
                    alpha = blade.outletSection[ii].alpha
                    Wa    = blade.outletSection[ii].Wa
                    Wt    = blade.outletSection[ii].Wt
                    W     = blade.outletSection[ii].W
                    beta  = blade.outletSection[ii].beta
                    # kinetics data allocation
                    self.inletSection[ii].Va    = Va
                    self.inletSection[ii].Vt    = Vt
                    self.inletSection[ii].V     = V
                    self.inletSection[ii].alpha = alpha
                    self.inletSection[ii].Wa    = Wa
                    self.inletSection[ii].Wt    = Wt
                    self.inletSection[ii].W     = W
                    self.inletSection[ii].beta  = beta
                    # thermodynamics allocation 
                    P     = blade.outletSection[ii].P
                    T     = blade.outletSection[ii].T
                    rho   = blade.outletSection[ii].rho
                    Pt    = blade.outletSection[ii].Pt
                    Tt    = blade.outletSection[ii].Tt
                    rhot  = blade.outletSection[ii].rhot
                    Ptr   = blade.outletSection[ii].Ptr
                    Ttr   = blade.outletSection[ii].Ttr
                    rhotr = blade.outletSection[ii].rhotr
                    a     = blade.outletSection[ii].a
                    M     = blade.outletSection[ii].M
                    Mr    = blade.outletSection[ii].Mr
                    s     = blade.outletSection[ii].s
                    # thermodynamics data allocation
                    self.inletSection[ii].P     = P
                    self.inletSection[ii].T     = T
                    self.inletSection[ii].rho   = rho
                    self.inletSection[ii].Pt    = Pt
                    self.inletSection[ii].Tt    = Tt
                    self.inletSection[ii].rhot  = rhot
                    self.inletSection[ii].Ptr   = Ptr
                    self.inletSection[ii].Ttr   = Ttr
                    self.inletSection[ii].rhotr = rhotr
                    self.inletSection[ii].a     = a
                    self.inletSection[ii].M     = M
                    self.inletSection[ii].Mr    = Mr
                    self.inletSection[ii].s     = s
                    # section shape properties 
                    gamma = blade.outletSection[ii].gamma
                    self.inletSection[ii].gamma = gamma

            elif toSection == 'outlet':
                for ii in range(self.nSection):
                    # kinetics allocation 
                    self.outletSection[ii].Va    = blade.outletSection[ii].Va
                    self.outletSection[ii].Vt    = blade.outletSection[ii].Vt
                    self.outletSection[ii].V     = blade.outletSection[ii].V
                    self.outletSection[ii].alpha = blade.outletSection[ii].alpha
                    self.outletSection[ii].Wa    = blade.outletSection[ii].Wa
                    self.outletSection[ii].Wt    = blade.outletSection[ii].Wt
                    self.outletSection[ii].W     = blade.outletSection[ii].W
                    self.outletSection[ii].beta  = blade.outletSection[ii].beta
                    # thermodynamics allocation 
                    self.outletSection[ii].P     = blade.outletSection[ii].P
                    self.outletSection[ii].T     = blade.outletSection[ii].T
                    self.outletSection[ii].rho   = blade.outletSection[ii].rho
                    self.outletSection[ii].Pt    = blade.outletSection[ii].Pt
                    self.outletSection[ii].Tt    = blade.outletSection[ii].Tt
                    self.outletSection[ii].rhot  = blade.outletSection[ii].rhot
                    self.outletSection[ii].Ptr   = blade.outletSection[ii].Ptr
                    self.outletSection[ii].Ttr   = blade.outletSection[ii].Ttr
                    self.outletSection[ii].rhotr = blade.outletSection[ii].rhotr
                    self.outletSection[ii].a     = blade.outletSection[ii].a
                    self.outletSection[ii].M     = blade.outletSection[ii].M
                    self.outletSection[ii].Mr    = blade.outletSection[ii].Mr
                    self.outletSection[ii].s     = blade.outletSection[ii].s
                    # section shape properties 
                    self.outletSection[ii].gamma = blade.outletSection[ii].gamma
        elif fromSection == 'inlet':
            if toSection == 'inlet':
                for ii in range(self.nSection):
                    # kinetics allocation 
                    self.inletSection[ii].Va    = blade.inletSection[ii].Va
                    self.inletSection[ii].Vt    = blade.inletSection[ii].Vt
                    self.inletSection[ii].V     = blade.inletSection[ii].V
                    self.inletSection[ii].alpha = blade.inletSection[ii].alpha
                    self.inletSection[ii].Wa    = blade.inletSection[ii].Wa
                    self.inletSection[ii].Wt    = blade.inletSection[ii].Wt
                    self.inletSection[ii].W     = blade.inletSection[ii].W
                    self.inletSection[ii].beta  = blade.inletSection[ii].beta
                    # thermodynamics allocation 
                    self.inletSection[ii].P     = blade.inletSection[ii].P
                    self.inletSection[ii].T     = blade.inletSection[ii].T
                    self.inletSection[ii].rho   = blade.inletSection[ii].rho
                    self.inletSection[ii].Pt    = blade.inletSection[ii].Pt
                    self.inletSection[ii].Tt    = blade.inletSection[ii].Tt
                    self.inletSection[ii].rhot  = blade.inletSection[ii].rhot
                    self.inletSection[ii].Ptr   = blade.inletSection[ii].Ptr
                    self.inletSection[ii].Ttr   = blade.inletSection[ii].Ttr
                    self.inletSection[ii].rhotr = blade.inletSection[ii].rhotr
                    self.inletSection[ii].a     = blade.inletSection[ii].a
                    self.inletSection[ii].M     = blade.inletSection[ii].M
                    self.inletSection[ii].Mr    = blade.inletSection[ii].Mr
                    self.inletSection[ii].s     = blade.inletSection[ii].s
                    # section shape properties 
                    self.inletSection[ii].gamma = blade.inletSection[ii].gamma
            elif toSection == 'outlet':
                for ii in range(self.nSection):
                    # kinetics allocation 
                    self.outletSection[ii].Va    = blade.inletSection[ii].Va
                    self.outletSection[ii].Vt    = blade.inletSection[ii].Vt
                    self.outletSection[ii].V     = blade.inletSection[ii].V
                    self.outletSection[ii].alpha = blade.inletSection[ii].alpha
                    self.outletSection[ii].Wa    = blade.inletSection[ii].Wa
                    self.outletSection[ii].Wt    = blade.inletSection[ii].Wt
                    self.outletSection[ii].W     = blade.inletSection[ii].W
                    self.outletSection[ii].beta  = blade.inletSection[ii].beta
                    # thermodynamics allocation 
                    self.outletSection[ii].P     = blade.inletSection[ii].P
                    self.outletSection[ii].T     = blade.inletSection[ii].T
                    self.outletSection[ii].rho   = blade.inletSection[ii].rho
                    self.outletSection[ii].Pt    = blade.inletSection[ii].Pt
                    self.outletSection[ii].Tt    = blade.inletSection[ii].Tt
                    self.outletSection[ii].rhot  = blade.inletSection[ii].rhot
                    self.outletSection[ii].Ptr   = blade.inletSection[ii].Ptr
                    self.outletSection[ii].Ttr   = blade.inletSection[ii].Ttr
                    self.outletSection[ii].rhotr = blade.inletSection[ii].rhotr
                    self.outletSection[ii].a     = blade.inletSection[ii].a
                    self.outletSection[ii].M     = blade.inletSection[ii].M
                    self.outletSection[ii].Mr    = blade.inletSection[ii].Mr
                    self.outletSection[ii].s     = blade.inletSection[ii].s
                    # section shape properties 
                    self.outletSection[ii].gamma = blade.inletSection[ii].gamma

    def velocityTriangles(self, sectionNumber, save=False, position='velocityTriangle.pgf'):
        '''
        This function plots the velocity triangles of a blade section.
            inputs:
                sectionNumber   -- the section vector to analise
                                -- it must be a vector 
        '''

        # generating axes
        if save:
            plt.rcParams['text.usetex'] = True

        fig, ax = plt.subplots(ncols=1, nrows=len(sectionNumber), figsize=(5,10))
        
        for ii,secNum in enumerate(sectionNumber):
            # setting up dimensions
            VaIn = self.inletSection[secNum].Va
            VtIn = self.inletSection[secNum].Vt
            VaOut = self.outletSection[secNum].Va
            VtOut = self.outletSection[secNum].Vt
            if self.turboType == 'rotor':
                WaIn = self.inletSection[secNum].Wa
                WtIn = self.inletSection[secNum].Wt
                WaOut = self.outletSection[secNum].Wa
                WtOut = self.outletSection[secNum].Wt
                U = self.inletSection[secNum].U

            # finding min and max velocities 
            if VaIn < VaOut:
                maxVelAx = VaOut
            else:
                maxVelAx = VaIn
            
            try:
                VtVec = [VtIn, VtOut, WtIn, WtOut]
            except:
                VtVec = [VtIn, VtOut]

            # computing maximum and minimum tangential velocity 
            minVelTan = np.min(VtVec)
            maxVelTan = np.max(VtVec)

            # plotting quiver 
            ax[ii].quiver(0, 0, VaIn, VtIn, color='firebrick', angles='xy', scale_units='xy', scale=1)
            ax[ii].quiver(0, 0, VaOut, VtOut, color='royalblue', angles='xy', scale_units='xy', scale=1)            
            try:
                # if the blade is from a rotor
                ax[ii].quiver(0, 0, WaIn, WtIn, color='forestgreen', angles='xy', scale_units='xy', scale=1)
                ax[ii].quiver(0, 0, WaOut, WtOut, color='darkorange', angles='xy', scale_units='xy', scale=1)
                ax[ii].quiver(WaIn, WtIn, 0, U, color='black', angles='xy', scale_units='xy', scale=1)
                ax[ii].quiver(WaOut, WtOut, 0, U, color='black', angles='xy', scale_units='xy', scale=1)
            except:
                pass

            ax[ii].set_xlim(-10, maxVelAx * 1.1)
            ax[ii].set_ylim(-50 + minVelTan * 1.1, maxVelTan * 1.1 + 50)
            try:
                ax[ii].set_title('section {0:d}\n'.format(secNum+1) + r'$U = {0:.2f} \frac{{m }}{{s }}$'.format(self.inletSection[secNum].U) + '\t' + r'$\Delta V_t = {0:.2f} \frac{{m }}{{s }}$'.format(VtOut - VtIn))
            except:
                ax[ii].set_title('section {0:d}\n'.format(secNum+1) + '\t' + r'$\Delta V_t = {0:.2f} \frac{{m }}{{s }}$'.format(VtOut - VtIn))
            ax[ii].grid(linestyle='--')
            ax[ii].set_xlabel(r'$V_a \ [\frac{{m }}{{s }}]$')   
            ax[ii].set_ylabel(r'$V_t \ [\frac{{m }}{{s }}]$')

        fig.suptitle('\n')
        if self.turboType == 'rotor':
            fig.legend(labels=[r'$V_{in }$', r'$V_{out }$', r'$W_{in }$', r'$W_{out }$', r'$U$'], loc='upper center', ncol=5, borderaxespad=0)
        else:
            fig.legend(labels=[r'$V_{in }$', r'$V_{out }$'], loc='upper center', ncol=3, borderaxespad=0)
        plt.tight_layout()

        if save:
            fig.savefig(position)
        else:
            plt.show()

    def computeBladeEfficiency(self, Va, lossVec, R=287.06, gamma=1.4):
        '''
        This function computes the efficiency of the blade.
            inputs:
                VaOut   -- isentropic axial outlet speed
        '''

        # cP computation
        cP = gamma / (gamma - 1) * R

        # allocating vector 
        eta = np.zeros(self.nSection)

        if self.turboType == 'rotor':
            # allocate outlet velocity 
            VaVec = [self.outletSection[ii].Va for ii in range(self.nSection)]
            Va = np.max(VaVec)

        # computing entropy generation through losses 
        for ii in range(self.nSection):    
            # new outlet pressure computation -> using losses
            if self.turboType == 'rotor':

                # data allocation 
                T1    = self.inletSection[ii].T
                W1    = self.inletSection[ii].W
                U1    = self.inletSection[ii].U
                W2iso = np.sqrt(self.outletSection[ii].Wt**2 + Va**2)
                U2    = self.outletSection[ii].U
                Ptr1  = self.inletSection[ii].Ptr
                P1    = self.inletSection[ii].P
                Ttr   = self.inletSection[ii].Ttr
                P2    = self.outletSection[ii].P
                
                # rothalpy computation
                R = cP * T1 + W1**2 / 2 - U1**2 / 2

                # T2 computation
                T2iso = (R - W2iso**2/2 + U2**2/2) / cP

                # Ttr 
                Ttr = T1 + W1**2/(2 * cP)

                # P2 computation
                P2 = Ptr1 * (T2iso/Ttr)**(gamma/(gamma-1))

                # Ptr2 computation
                Ptr2 = Ptr1 - lossVec[ii] * (Ptr1 - P1)

                # T2 computation
                T2 = Ttr * (P2/Ptr2)**((gamma-1)/gamma)

                # W2 computation
                W2 = np.sqrt(2 * cP * (Ttr - T2))

                # computing efficiency
                eta[ii] = (W1**2 - W2iso**2) / (W1**2 - W2**2)   

            elif self.turboType == 'stator':
                # data allocation 
                Pt1   = self.inletSection[ii].Pt
                Tt1   = self.inletSection[ii].Tt
                Tt2   = Tt1
                P1    = self.inletSection[ii].P
                V1    = np.sqrt(Va**2 + self.inletSection[ii].Vt**2)
                T1    = self.inletSection[ii].T
                V2iso = np.sqrt(Va**2 + self.outletSection[ii].Vt**2)

                # due to losses there is a reduction in the relative total pressure
                Pt2 = Pt1 - lossVec[ii] * (Pt1 - P1)

                # computing T2iso 
                T2iso = T1 + V1**2/(2*cP) - V2iso**2/(2*cP)

                # computing exit pressure 
                P2 = Pt1 * (T2iso/Tt1)**(gamma/(gamma - 1))

                # computing T2 real 
                T2 = Tt2 * (P2/Pt2)**((gamma - 1)/gamma)

                # computing V2 
                V2 = np.sqrt((T1 - T2)*2*cP + V1**2) 

                # computing efficiency
                eta[ii] = (V1**2 - V2iso**2) / (V1**2 - V2**2)   
                
        # efficiency
        etaBlade = np.sum(eta) / self.nSection

        # printing results
        starDim = 82
        efficiencyDim = np.int16(np.floor(starDim - len(' EFFICIENCY '))/2)
        print('*' * efficiencyDim + ' EFFICIENCY ' + '*' * efficiencyDim)
        print('-- etaBlade           = {0:>4.3f}'.format(etaBlade))
        print('-- minimum efficiency = {0:>4.3f} -- minimum efficiency position = {1:d}'.format(np.min(eta), np.argmin(eta)))
        print('-- maximum efficiency = {0:>4.3f} -- maximum efficiency position = {1:d}'.format(np.max(eta), np.argmax(eta)))
        print('*' * starDim)

    def computeMeanPressure(self, kind='total'):
        '''
        This function computes the mean pressure static/total of the outlet section of a blade
            inputs:
                kind    -- static/total
        '''          

        # field points allocation 
        if kind == 'total':
            # vector allocation 
            PVec = [self.outletSection[ii].Pt for ii in range(self.nSection)]
        elif kind == 'static':
            # vector allocation 
            PVec = [self.outletSection[ii].P for ii in range(self.nSection)]
        
        # field interpolation 
        outletMidpoint = [self.outletSection[ii].midpoint for ii in range(self.nSection)]
        PressureFunc   = interpolate.interp1d(outletMidpoint, PVec, kind='linear', bounds_error=False, fill_value='extrapolate')

        # setting up integration extremes 
        hub = self.outletSection[0].bottom
        tip = self.outletSection[-1].tip

        # computing mean pressure 
        Pmean = integrate.quad(PressureFunc, hub, tip, epsrel=1e-3)[0]/(tip - hub)

        # printing results
        starDim = 82
        pressureDim = np.int16(np.floor(starDim - len(' MEAN PRESSURE '))/2)
        print('*' * pressureDim + ' MEAN PRESSURE ' + '*' * pressureDim)
        print('-- mean ', kind, ' pressure = {0:>4.3f}'.format(Pmean))
        print('*' * starDim)

        return Pmean