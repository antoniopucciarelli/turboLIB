# TURBOMACHINERY -- ENGINEERING DIMENSION AND COEFFICIENT FUNCTIONS LIBRARY
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY: compressor/turbine
#   CONTENT: engineering coefficient functions for the computation of turbomachinery performace
#  

# WORK COMPUTATION 

def L_is(Tin=0, Pin=0, Pout=0, beta=0, kind='', gamma=1.4, R=287.06):
    '''
    Compressor work computation
        working conditions: isentropic transformation

        function inputs:
            Tin         [K]         -- initial temperature 
            Pin         [Pa]        -- initial pressure
            Pout        [Pa]        -- final pressure 
            beta        [--]        -- pressure ratio 
            kind                    -- turbomachinery type 
                                    -- values:
                                        -- compressor
                                        -- turbine 
            gamma       [--]        -- gas heat capacity ratio
            R = 287.06  [J/kg K]    -- air gas constant Ru/Mair 
    '''

    # there are 2 ways for L computation
    # 1. using directly beta
    # 2. computing directly beta from Pin and Pout 
    if Pin != 0 and Pout !=0: 
        if kind == 'compressor':
            beta = Pout/Pin
        elif kind == 'turbine':
            beta = Pin/Pout 

    # work computation
    if kind == 'compressor':
        Lis = (gamma * R)/(gamma - 1) * Tin * (beta**((gamma - 1)/gamma) - 1)
    elif kind == 'turbine':
        Lis = (gamma * R)/(gamma - 1) * Tin * (1 - beta**((1 - gamma)/gamma))

    return Lis

def L_poly(Tin=0, Pin=0, Pout=0, beta=0, kind='', n=0, R=287.06):
    '''
    Compressor work computation
        working conditions: polytropic transformation

        function inputs:
            Tin         [K]         -- initial temperature 
            Pin         [Pa]        -- initial pressure
            Pout        [Pa]        -- final pressure 
            beta        [--]        -- pressure ratio 
            kind                    -- turbomachinery type 
                                    -- values:
                                        -- compressor
                                        -- turbine 
            n           [--]        -- gas polytropic index
            R = 287.06  [J/kg K]    -- air gas constant Ru/Mair 
    '''

    # there are 2 ways for L computation
    # 1. using directly beta
    # 2. computing directly beta from Pin and Pout 
    if Pin != 0 and Pout !=0: 
        if kind == 'compressor':
            beta = Pout/Pin
        elif kind == 'turbine':
            beta = Pin/Pout 

    # work computation
    if kind == 'compressor':
        Lpoly = (n * R * Tin)/(n - 1) * (beta**((n - 1)/n) - 1) 
    elif kind == 'turbine':
        Lpoly = (n * R * Tin)/(n - 1) * (1 - beta**((1 - n)/n))
    
    return Lpoly

def L_adR(Tin=0, Pin=0, Pout=0, n=0, beta=0, kind='', gamma=1.4, R=287.06):
    '''
    Compressor work computation
        working conditions: polytropic transformation

        function inputs:
            Tin         -- initial temperature 
            Pin         -- initial pressure
            Pout        -- final pressure 
            beta        -- pressure ratio 
            kind        -- turbomachinery type 
                        -- values:
                            -- compressor
                            -- turbine 
            n           -- polytropic index
            gamma       -- gas heat capacity ratio 
            R = 287.06  -- air gas constant Ru/Mair 
    '''

    # there are 2 ways for L computation
    # 1. using directly beta
    # 2. computing directly beta from Pin and Pout 
    if Pin != 0 and Pout !=0: 
        if kind == 'compressor':
            beta = Pout/Pin
        elif kind == 'turbine':
            beta = Pin/Pout 

    # work computation
    if kind == 'compressor':
        LadR = (gamma * R * Tin)/(gamma - 1) * (beta**((n - 1)/n) - 1)
    elif kind == 'turbine':
        LadR = (gamma * R * Tin)/(gamma - 1) * (1 - beta**((1 - n)/n))
    
    return LadR

# EFFICIENCY COMPUTATION

def stageEfficiency(rotorBlade, statorBlade, R=287.06, gamma=1.4):
    '''
    This function computes the efficiency of a compressor stage.
        inputs:
            rotorBlade  -- turboBlade object of rotor type
            statroBlade -- turboBlade object of stator type 

    '''

    # compute stage efficiency 
    etaStage = 0 

    # cP computation 
    cP = gamma / (gamma - 1) * R 

    for ii in range(rotorBlade.nSection):
        # allocating total pressure 
        Pt1 = rotorBlade.inletSection[ii].Pt
        Pt2 = statorBlade.outletSection[ii].Pt 
        
        # allocation total temperature 
        Tt1 = rotorBlade.inletSection[ii].Tt
        
        # computing total temperature if the transformation Pt1 -> Pt2 were isentropic 
        Tt2 = Tt1 * (Pt2/Pt1)**((gamma-1)/gamma)
        
        # computing work if the transformation were isentropic 
        Lis = cP * (Tt2 - Tt1)
        
        # euler work computation 
        Leu = cP * (rotorBlade.outletSection[ii].Tt - rotorBlade.inletSection[ii].Tt)
        
        # etaSection 
        etaSection = Lis / Leu
        
        # etaStage 
        etaStage = etaStage + etaSection

    # computing total stage efficiency
    etaStage = etaStage / rotorBlade.nSection

    print('\n-- STAGE EFFICIENCY: etaStage = {0:>4.4f}'.format(etaStage))

    return etaStage

def ise_eff(hout_is=0, hin=0, hout=0, Tout_is=0, Tin=0, Tout=0, beta=0, kind='', gamma=0, n=0, L_is=0, L_adR=0):
    '''
    Isentropic efficiency computation
        function inputs:
            hout_is     -- final enthalpy for an isentropic transformation
            hin         -- initial enthaply 
            hout        -- final enthalpy for the real transformation
            Tout_is     -- final temperature for an isentropic transformation
            Tin         -- initial temperature
            Tout        -- final temperature for the real transformation
            beta        -- pressure ratio
            kind        -- turbomachinery type 
                        -- values:
                            -- compressor
                            -- turbine 
            gamma       -- gas heat capacity ratio
            n           -- polytropic index
            L_is        -- isentropic compressor work
            L_adR       -- adiabatic real compressor work 
        
        computation -- there are different ways for the computation of the efficiency
            enthalpy based 
            temperature based 
            pressure ratio based 
            work based
    '''

    if kind == 'compressor':
        # enthalpy based computation
        if hout_is != 0 and hin != 0 and hout != 0:
            eta = (hout_is - hin) / (hout - hin)

        # temperature based 
        if Tout_is != 0 and Tin != 0 and Tout != 0:
            eta = (Tout_is - Tin) / (Tout - Tin)

        # pressure ratio based 
        if beta != 0 and gamma != 0 and n != 0: 
            eta = (beta**((gamma - 1)/gamma) - 1) / (beta**((n - 1)/n) - 1)

        # work based 
        if L_is != 0 and L_adR != 0:
            eta = L_is / L_adR 
   
    elif kind == 'turbine':
        # enthalpy based computation
        if hout_is != 0 and hin != 0 and hout != 0:
            eta = (hout - hin) / (hout_is - hin)

        # temperature based 
        if Tout_is != 0 and Tin != 0 and Tout != 0:
            eta = (Tout - Tin) / (Tout_is - Tin)

        # pressure ratio based 
        if beta != 0 and gamma != 0 and n != 0: 
            eta = (1 - beta**((1 - n)/n)) / (1 - beta**((1 - gamma)/gamma))

        # work based 
        if L_is != 0 and L_adR != 0:
            eta = L_adR / L_is 

    return eta    

def poly_eff(gamma=0, n=0, L_p=0, L_adR=0, kind=''):
    '''
    Polytropic efficiency computation
        function inputs:
            gamma       -- gas heat capacity ratio
            n           -- polytropic index
            L_p         -- polytropic transformation compressor work
            L_adR       -- adiabatic real compressor work 
            kind        -- turbomachinery type 
                        -- values:
                            -- compressor
                            -- turbine 
        
        computation -- there are different ways for the computation of the efficiency
            work based
            transformation factor based 
    '''

    if kind == 'compressor':
        # work based 
        if L_p != 0 and L_adR != 0:
            eta = L_p / L_adR 

        # transformation factor based
        if gamma != 0 and n != 0:
            eta = (n * (gamma - 1)) / (gamma * (n - 1))
    elif kind == 'turbine':
        # work based 
        if L_p != 0 and L_adR != 0:
            eta = L_adR / L_p

        # transformation factor based
        if gamma != 0 and n != 0:
            eta = (gamma * (n - 1)) / (n * (gamma - 1))

    return eta    

def reHeat_factor(eta_poly, eta_is, kind=''):
    '''
    Re-heat factor computation
        function inputs:
            eta_poly    -- polytropic transformation efficiency
            eta_is      -- insentropic trnasformation efficiency
            kind        -- turbomachinery type 
                        -- values:
                            -- compressor
                            -- turbine 
    '''

    if kind == 'compressor':
        reHeatFactor = eta_poly / eta_is - 1
    elif kind == 'turbine':
        reHeatFactor = eta_is / eta_poly - 1
    
    return reHeatFactor