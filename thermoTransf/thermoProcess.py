# TURBOMACHINERY -- LIBRARY FOR THE PLOTTING OF THERMODYNAMIC PROCESSES
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   THERMODYNAMICS PROCESSES: plot and computation
#       -- isothermal process
#       -- isentropic process
#       -- isobaric process
#       -- polytropic process
#       -- adiabatic real process 
#
#   FLUID: air
#  

# importing libraries
import pyromat as pm 
import numpy as np
import matplotlib.pyplot as plt 

# setting up dimensions
unit_pressure = pm.config['unit_pressure'] = 'bar'
unit_temperature = pm.config['unit_temperature'] = 'K'
unit_energy = pm.config['unit_energy'] = 'kJ'
unit_mass = pm.config['unit_mass'] = 'kg'

# air object generation 
air = pm.get('ig.air')

def isothermal(T=0,pIn=0,pOut=0,dim=100,name='isothermal',color='k',plot=True):
    '''
    Isothermal transformation plot

        function inputs:
            T       -- constant temperature during the process
            Pin     -- initial pressure 
            Pout    -- final pressure 
            dim     -- number of points for the transformation plot 
                    -- default set to 100
            name    -- transformation name
            color   -- plotting color
    '''

    # getting air data 
    global air

    # getting dimensions
    global unit_energy, unit_mass, unit_pressure, unit_temperature

    if pIn != 0 and pOut != 0:
        # generating array of value for the pressure 
        pVec = np.linspace(pIn, pOut, dim)

        # computing entropy values for the isothermal transformation
        sVec = air.s(T=T,p=pVec)

        # computing enthalpy values for the isothermal transformation
        hVec = air.h(T=T,p=pVec)
    
    else:
        print('input error in:', name)

    # initial values
    hIn = hVec[0]
    sIn = sVec[0]
    pIn = pVec[0]
    Tin = T

    # final values
    hOut = hVec[-1]
    sOut = sVec[-1]
    pOut = pVec[-1]
    Tout = T

    # printing of useful values
    print('+++++++++++++++++++++++++++++++++++++++++++')
    print('Isothermal transformation: ', name)
    print('-- Initial values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sIn, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hIn, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pIn, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tin, unit_temperature))
    print('-- Final values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sOut, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hOut, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pOut, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tout, unit_temperature))
    print('+++++++++++++++++++++++++++++++++++++++++++\n')

    # plotting transformation
    if plot:
        plt.plot(sVec, hVec, color=color, label=name)

    return [hIn, sIn, pIn, Tin], [hOut, sOut, pOut, Tout]

def isentropic(s=0,pIn=0,pOut=0,Tin=0,Tout=0,dim=100,name='isentropic',color='k',plot=True):
    '''
    Isentropic transformation plot

        function inputs:
            s       -- constant entropy during the process
            Pin     -- initial pressure | -> pressure based 
            Pout    -- final pressure   |
            Tin     -- intial temperature | -> temperature based
            Tout    -- final temperature  | 
            dim     -- number of points for the transformation plot 
                    -- default set to 100
            name    -- transformation name
            color   -- plotting color
            plot    -- boolean for plotting
    '''

    # getting air data 
    global air

    # getting dimensions
    global unit_energy, unit_mass, unit_pressure, unit_temperature

    # generation of variables arrays
    pVec = np.linspace(pIn, pOut, dim)
    Tvec = np.linspace(Tin, Tout, dim)

    # computing transformation with respect to the input data 
    if pIn != 0 and pOut != 0:
        # temperature vector allocation 
        Tvec = np.zeros([dim])
        # temperature computation
        for ii, p in enumerate(pVec):
            Tvec[ii] = air.T_s(s=s,p=p)

        # enthalpy computation
        hVec = air.h(T=Tvec, p=pVec)

    elif Tin != 0 and Tout != 0:
        # pressure vector allocation
        pVec = np.zeros([dim])
        # pressure computation
        for ii, T in enumerate(Tvec):
            pVec[ii] = air.p_s(s=s, T=T)
        
        # enthalpy computation
        hVec = air.h(T=Tvec, p=pVec)

    else:
        print('input error in:', name)

    # initial values
    hIn = hVec[0]
    sIn = s
    pIn = pVec[0]
    Tin = Tvec[0]

    # final values
    hOut = hVec[-1]
    sOut = s
    pOut = pVec[-1]
    Tout = Tvec[-1]

    # printing of useful values
    print('+++++++++++++++++++++++++++++++++++++++++++')
    print('Isentropic transformation: ', name)
    print('-- Initial values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sIn, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hIn, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pIn, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tin, unit_temperature))
    print('-- Final values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sOut, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hOut, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pOut, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tout, unit_temperature))
    print('+++++++++++++++++++++++++++++++++++++++++++\n')

    # plotting transformation
    if plot:
        plt.plot(s*np.ones([dim]), hVec, color=color, label=name)

    return [hIn, sIn, pIn, Tin], [hOut, sOut, pOut, Tout]

def isobaric(p=0,Tin=0,Tout=0,dim=100,name='isobaric',color='k',plot=True):
    '''
    Isobaric transformation plot

        function inputs:
            p       -- constant pressure during the process
            Tin     -- intial temperature | -> temperature based
            Tout    -- final temperature  | 
            dim     -- number of points for the transformation plot 
                    -- default set to 100
            name    -- transformation name
            color   -- plotting color
            plot    -- boolean for plotting
    '''

    # getting air data 
    global air

    if Tin != 0 and Tout != 0:
        # array allocation 
        sVec = np.zeros([dim])
        hVec = np.zeros([dim])
        
        # temperature vector generation 
        Tvec = np.linspace(Tin, Tout, dim)

        # values computation
        for ii, T in enumerate(Tvec):
            sVec[ii] = air.s(T=T,p=p)
            hVec[ii] = air.h(T=T,p=p)
    
    else:
        print('input error in: ', name)
    
    # initial values
    hIn = hVec[0]
    sIn = sVec[0]
    pIn = p
    Tin = Tvec[0]

    # final values
    hOut = hVec[-1]
    sOut = sVec[-1]
    pOut = p
    Tout = Tvec[-1]

    # printing of useful values
    print('+++++++++++++++++++++++++++++++++++++++++++')
    print('Isobaric transformation: ', name)
    print('-- Initial values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sIn, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hIn, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pIn, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tin, unit_temperature))
    print('-- Final values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sOut, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hOut, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pOut, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tout, unit_temperature))
    print('+++++++++++++++++++++++++++++++++++++++++++\n')

    # plotting transformation
    if plot:
        plt.plot(sVec, hVec, color=color, label=name)

    return [hIn, sIn, pIn, Tin], [hOut, sOut, pOut, Tout]

def polytropic(n=0,pIn=0,pOut=0,Tin=0,Tout=0,R=8.3145e-2,dim=100,name='polytropic',color='k',plot=True):
    '''
    Isobaric transformation plot

        function inputs:
            n       -- polytropic thermodynamics process index
            Tin     -- intial temperature 
            Tout    -- final temperature   
            pIn     -- initial pressure   
            pOut    -- final pressure  
            R       -- perfect gas constant
                    -- !!! relate it to the engineering unit of measure adopted
            dim     -- number of points for the transformation plot 
                    -- default set to 100
            name    -- transformation name
            color   -- plotting color
            plot    -- boolean for plotting
    '''

    # getting air data 
    global air

    # getting dimensions
    global unit_energy, unit_mass, unit_pressure, unit_temperature
    
    if n != 0 and pIn != 0 and pOut != 0:
        # pressure vector generation 
        pVec = np.linspace(pIn, pOut, dim)

        # different cases with respect to the known flow properties
        if Tin != 0:
            dIn = air.d(T=Tin, p=pIn)
            dVec = (pVec / pIn)**(1/n) * dIn
        elif Tout != 0:
            dOut = air.d(T=Tout, p=pOut)
            dVec = (pVec / pOut)**(1/n) * dOut
        else:
            print('input error on temperature in: ', name)

        # temperature vector allocation 
        Tvec = np.zeros([dim])
        for ii in range(dim):
            # temperature computation
            # !!! PYroMat doesn't allow the inverse computation of T through the density 
            # !!! in order to solve this issue ideal gas relation has been used for the computation of T
            Tvec[ii] = pVec[ii] / ( dVec[ii] * R / air.mw() )

        # enthaply and entropy vector allocation
        hVec = np.zeros([dim])
        sVec = np.zeros([dim])
        
        # enthalpy and entropy computation
        for ii in range(dim):
            hVec[ii] = air.h(T=Tvec[ii], p=pVec[ii])
            sVec[ii] = air.s(T=Tvec[ii], p=pVec[ii])

    elif n != 0 and Tin != 0 and Tout != 0:
        # temperature vector generation 
        Tvec = np.array(Tin, Tout, dim)

        # different cases with respect to the known flow properties
        if pIn != 0:
            dIn = air.d(T=Tin, p=pIn)
            dVec = (pVec / pIn)**(1/n) * dIn
        elif pOut != 0:
            dOut = air.d(T=Tout, p=pOut)
            dVec = (pVec / pOut)**(1/n) * dOut

        # pressure vector allocation 
        pVec = np.zeros([dim])
        for ii in range(dim):
            # pressure computation
            # !!! PYroMat does not allow the inverse computation of p through the density 
            # !!! in order to solve this issue ideal gas relation has been used for the computation of p
            pVec[ii] = Tvec[ii] * (R / air.mw()) * dVec[ii] 

        # enthaply and entropy vector allocation
        hVec = np.zeros([dim])
        sVec = np.zeros([dim])
        
        # enthalpy and entropy computation
        for ii in range(dim):
            hVec[ii] = air.h(T=Tvec[ii], p=pVec[ii])
            sVec[ii] = air.s(T=Tvec[ii], p=pVec[ii])

    else:
        print('input error in: ', name)

    # plotting transformation 
    if plot:
        plt.plot(sVec, hVec, color=color, label=name)

    # initial values
    hIn = hVec[0]
    sIn = sVec[0]
    pIn = pVec[0]
    Tin = Tvec[0]

    # final values
    hOut = hVec[-1]
    sOut = sVec[-1]
    pOut = pVec[-1]
    Tout = Tvec[-1]

    # printing of useful values
    print('+++++++++++++++++++++++++++++++++++++++++++')
    print('Polytropic transformation: ', name)
    print('-- n = {0}'.format(n))
    print('-- R = {0:f}'.format(R))
    print('-- Initial values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sIn, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hIn, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pIn, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tin, unit_temperature))
    print('-- Final values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sOut, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hOut, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pOut, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tout, unit_temperature))
    print('+++++++++++++++++++++++++++++++++++++++++++\n')

    return [hIn, sIn, pIn, Tin], [hOut, sOut, pOut, Tout]

def adiabatic(eta=1,pIn=0,pOut=0,Tin=0,s=0,dim=100,name='adiabatic',color='k',plot=True):
    '''
    Adiabatic transformation plot

        function inputs:
            eta     -- transformation efficiency 
                    -- default set to 1 
            pIn     -- initial pressure   
            pOut    -- final pressure  
            Tin     -- intial temperature  
            s       -- constant entropy for the ideal process 
            dim     -- number of points for the transformation plot 
            R       -- perfect gas constant
                    -- !!! relate it to the engineering unit of measure adopted
                    -- default set to 100
            name    -- transformation name
            color   -- plotting color
            plot    -- boolean for plotting
    '''

    # getting air data
    global air 

    # getting dimensions
    global unit_energy, unit_mass, unit_pressure, unit_temperature

    # checking if the transformation is for a compressor or a turbine
    if pOut > pIn: 
        # compressor
        process = 'compression'
    elif pOut < pIn:
        # turbine 
        process = 'expansion'
    else:
        print('input error in: ', name)

    # generation of pressure vector 
    pVec = np.linspace(pIn, pOut, dim)

    # enthalpy computation 
    if Tin != 0 or s != 0:
        if Tin !=0 :    
            # entropy computation -- isentropic transformation 
            s = air.s(T=Tin, p=pIn)
        
        # temperature vector allocation -- isentropic transformation
        Tvec_is = np.zeros([dim])

        # temperature computation -- isentropic transformation
        for ii, p in enumerate(pVec):
            Tvec_is[ii] = air.T_s(s=s, p=p)

        # enthalpy vector allocation -- isentropic transformation
        hVec_is = np.zeros([dim])

        # enthalpy computation -- isentropic transformation
        for ii in range(dim):
            hVec_is[ii] = air.h(T=Tvec_is[ii], p=pVec[ii])

        # enthaply vector allocation -- real transformation 
        hVec = np.zeros([dim])

        # enthalpy computation -- real transformation 
        if pOut < pIn:
            # if the trasformation is relative to a turbine
            # eta = (h1 - h2) / (h1 - h2_is)
            # h2 = eta * h2_is + (1 - eta) * h1
            for ii in range(dim):
                hVec[ii] = eta * hVec_is[ii] + (1 - eta) * hVec_is[0]
        elif pIn < pOut: 
            # if the transformation is relative to a compressor
            # eta = (h2_is - h1) / (h2 - h1)
            # h2 = (h2_is - h1) / eta + h1
            for ii in range(dim):
                hVec[ii] = ( hVec_is[ii] - hVec_is[0] ) / eta + hVec_is[0]
        else:
            print('input error in: ', name)

        # temperature vector allocation -- real transformation 
        Tvec = np.zeros([dim])

        # temperature computation -- real transformation 
        for ii in range(dim):
            Tvec[ii] = air.T_h(h=hVec[ii],p=pVec[ii])

        # entropy vector allocation -- real transformation 
        sVec = np.zeros([dim])

        # entropy computation -- real transformation 
        for ii in range(dim):
            sVec[ii] = air.s(T=Tvec[ii], p=pVec[ii])
  
    else:
        print('input error in: ', name)

    # plotting transformation 
    if plot:
        plt.plot(sVec, hVec, color=color, label=name)

    # initial values
    hIn = hVec[0]
    sIn = sVec[0]
    pIn = pVec[0]
    Tin = Tvec[0]

    # final values
    hOut = hVec[-1]
    sOut = sVec[-1]
    pOut = pVec[-1]
    Tout = Tvec[-1]

    # printing of useful values
    print('+++++++++++++++++++++++++++++++++++++++++++')
    print('Adiabatic transformation: ', name)
    print('-- eta     = {0}'.format(eta))
    print('-- process = {0:s}'.format(process))
    print('-- Initial values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sIn, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hIn, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pIn, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tin, unit_temperature))
    print('-- Final values:')
    print('\ts = {0:>10.3f} {1:s}/{2:s}'.format(sOut, unit_energy, unit_mass))
    print('\th = {0:>10.3f} {1:s}/{2:s}'.format(hOut, unit_energy, unit_mass))
    print('\tp = {0:>10.3f} {1:s}'.format(pOut, unit_pressure))
    print('\tT = {0:>10.3f} {1:s}'.format(Tout, unit_temperature))
    print('+++++++++++++++++++++++++++++++++++++++++++\n')

    return [hIn, sIn, pIn, Tin], [hOut, sOut, pOut, Tout]
