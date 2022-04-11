# importing libraries
from turboClass import turboBlade

# objects generation 
rotorBlade = turboBlade.blade(ID=1, turboType='rotor', nSection=10, inletBladeHeight=1, outletBladeHeight=1, inletHubRadius=1, outletHubRadius=1, airfoilPath='data/airfoils/naca65.txt')
statorBlade = turboBlade.blade(ID=1, turboType='stator', nSection=10, inletBladeHeight=1, outletBladeHeight=1, inletHubRadius=1, outletHubRadius=1, airfoilPath='data/airfoils/naca65.txt')

# 

