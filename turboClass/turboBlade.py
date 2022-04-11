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

    def allocateThermodynamics(self, Tt0, Pt0, Leu, eta, R=287.06, gamma=1.4):
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

            # total temperature computation
            Tt1 = Leu / cP + Tt0

            # ideal temperature computation if the process is completely 
            # isentropic without losses but the work produced is related 
            # to a process that takes into account losses in the stage  
            T1 = Tt1 - self.outletSection[ii].V**2 / (2 * cP)

            # T1 isoentropic computation 
            #   this correction activates only is eta != 1
            T1 = T0 + eta * (T1 - T0)

            # pressure computation
            P1 = P0 * (T1/T0)**(gamma/(gamma-1))

            # total pressure computation
            Pt1 = P1 * (Tt1/T1)**(gamma/(gamma-1))

            # density computation 
            rho1 = P1 / (R * T1)

            # total density computation
            rhot1 = Pt1 / (R * Tt1)

            # variable allocation in seection objects            
            self.inletSection[ii].allocateThermodynamics(T=T0, P=P0, Tt=Tt0, Pt=Pt0, rho=rho0, rhot=rhot0, s=0)
            self.outletSection[ii].allocateThermodynamics(T=T1, P=P1, Tt=Tt1, Pt=Pt1, rho=rho1, rhot=rhot1, s=0)

    def radialEquilibrium(self, R=287.06, gamma=1.4):
        '''
        This function computes the radial equilibrium of the section taking into account losses. 
            inputs:
        ''' 
    
        # cP computation
        cP = gamma / (gamma - 1) * R

        for ii in range(self.nSection):
            # variable allocation
            ht0 = self.inletSection[ii].ht
            Vt0 = self.inletSection[ii].Vt
            U0 = self.inletSection[ii].U
            s0 = self.inletSection[ii].s
            Vt1 = self.outletSection[ii].Vt
            U1 = self.outletSection[ii].U 
            s1 = self.outletSection[ii].s 

            # total enthalpy computation 
            ht1 = ht0 + U1 * Vt1 - U0 * Vt0
            self.outletSection[ii].ht = ht1 
            
            # total temperature computation
            self.outletSection[ii].Tt = (ht1 - ht0) / cP + self.inletSection[ii].Tt



