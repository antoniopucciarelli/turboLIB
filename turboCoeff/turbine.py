# TURBOMACHINERY -- ENGINEERING DIMENSION AND COEFFICIENT FUNCTIONS LIBRARY
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY: turbine
#   CONTENT: engineering coefficient functions for the computation of turbomachinery performace
#  

# WORK COMPUTATION 

def L_is(Tin,Pin,Pout,gamma,R=8.3145):
    '''
    Turbine work computation
        working conditions: isentropic transformation

        function inputs:
            Tin        -- initial temperature 
            Pin        -- initial pressure
            Pout       -- final pressure 
            gamma      -- gas heat capacity ratio
            R = 8.3145 -- perfect gas universal constant 
    '''
    return (gamma * R * Tin)/(gamma - 1) * (1 - (Pin/Pout)**((1 - gamma)/gamma))

def L_poly(Tin,Pin,Pout,n,R=8.3145):
    '''
    Turbine work computation
        working conditions: polytropic transformation

        function inputs:
            Tin        -- initial temperature 
            Pin        -- initial pressure
            Pout       -- final pressure 
            n          -- polytropic index 
            R = 8.3145 -- perfect gas universal constant 
    '''
    return (n * R * Tin)/(n - 1) * (1 - (Pin/Pout)**((1 - n)/n))

def L_adR(Tin,Pin,Pout,n,gamma,R=8.3145):
    '''
    Turbine work computation
        working conditions: polytropic transformation

        function inputs:
            Tin        -- initial temperature 
            Pin        -- initial pressure
            Pout       -- final pressure 
            n          -- polytropic index
            gamma      -- gas heat capacity ratio 
            R = 8.3145 -- perfect gas universal constant 
    '''
    return (gamma * R * Tin)/(gamma - 1) * (1 - (Pin/Pout)**((1 - n)/n))

# EFFICIENCY COMPUTATION

def ise_eff(hout_is=0,hin=0,hout=0,Tout_is=0,Tin=0,Tout=0,beta=0,gamma=0,n=0,L_is=0,L_adR=0):
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

def poly_eff(gamma=0,n=0,L_p=0,L_adR=0):
    '''
    Polytropic efficiency computation
        function inputs:
            gamma       -- gas heat capacity ratio
            n           -- polytropic index
            L_p         -- polytropic transformation compressor work
            L_adR       -- adiabatic real compressor work 
        
        computation -- there are different ways for the computation of the efficiency
            work based
            transformation factor based 
    '''

    # work based 
    if L_p != 0 and L_adR != 0:
        eta = L_adR / L_p

    # transformation factor based
    if gamma != 0 and n != 0:
        eta = (gamma * (n - 1)) / (n * (gamma - 1))

    return eta    

def reHeat_factor(eta_poly,eta_is):
    '''
    Re-heat factor computation
        function inputs:
            eta_poly    -- polytropic transformation efficiency
            eta_is      -- insentropic trnasformation efficiency
    '''

    return eta_is / eta_poly - 1