# importing libraries
from turboClass import turboBlade

# objects generation 
rotorBlade = turboBlade.blade(ID=1, turboType='rotor', nSection=100, inletBladeHeight=1, outletBladeHeight=1, inletHubRadius=1, outletHubRadius=1, omega=1)
statorBlade = turboBlade.blade(ID=1, turboType='stator', nSection=10, inletBladeHeight=1, outletBladeHeight=1, inletHubRadius=1, outletHubRadius=1, omega=1)

rotorBlade.allocateDynamics(rMean=1.5, VtMean=100, VaMean=100, omega=10, section='outlet')
rotorBlade.plotMeridional()

# 

