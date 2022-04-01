# AXIAL TURBOMACHINERY DESIGN
# importing libraries 
import numpy as np 
from thermoTransf.thermoProcess import *
from turboCoeff.similarity import *
import matplotlib.pyplot as plt 

# data
# inlet values
Pt0 = 1e+5   # inlet total pressure     [Pa]
Tt0 = 300    # inlet total temperature  [K]
# constraints
mFlux = 100  # mass flux                [kg/s]
betaP = 1.45 # compression ratio        [--]
maxR  = 0.9  # maximum tip radius       [m]

# air properties allocation 
R     = 287.06                  # air gas constant Ru / Mm      [J/kg K]
gamma = 1.4                     # specific heat ratio           [--]
cP    = gamma / (gamma - 1) * R # specific heat ratio @ P cost  [J/kg K]
cV    = cP - R                  # specific heat ratio @ V cost  [J/kg K]

# stage hypothesis
eta = 0.90 # stage efficiency 
rD  = 0.5  # reaction degree 

# inlet velocity --> the velocity is straight -- no guide vanes are present 
# 
