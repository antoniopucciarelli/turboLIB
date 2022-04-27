# importing libraries
from turboClass import turboBlade
from turboCoeff import similarity
from turboCoeff import coeff
from turboClass import bladeStudy
from geometry import bladeGenerator
import numpy as np
import contextlib

# data
# constraints
mFlux = 100     # mass flux                [kg/s]
betaP = 1.45    # compression ratio        [--]
maxD  = 0.9     # maximum tip diameter     [m]
maxR  = maxD/2  # maximum tip radius       [m]

# inlet values
Pt0 = 1e+5      # inlet total pressure     [Pa]
Tt0 = 300       # inlet total temperature  [K]

# stage hypothesis
# reaction degree
rD = 0.55
# stage mean radius -> radius @ inlet blade midspan
rMean = 0.325
# rotor inlet tangential velocity
Vt0Umean = 0

# Vt0 = (1 - rD - lam/4) * Umean -> lam = 4 * (1 - rD - Vt0/Umean) 
# psi = lam / 2 
lam = (1 - rD - Vt0Umean) * 4 
psiTarget = lam / 2

# declaring blades
#nRotorBlades  = 28
#nStatorBlades = 28
#nSection      = 50

# output file generation
#file_path = 'compressor_' + str(rD) + '_' + str(rMean) + '_' + str(nRotorBlades) + '_' + str(nStatorBlades) + '.txt'
#with open(file_path, "w") as file:
#    with contextlib.redirect_stdout(file):
# generation of mean line properties to be used for the blade assembly 
adimVec, bladeVec, rotationVec, V0vec, V1vec, V2vec, W0vec, W1vec, W2vec, thermo0, thermo1, _, work = similarity.stageProperties(rD, psiTarget, rMean, mFlux, Tt0, Pt0, betaP, T1real=True, printout=True, save=False)

# values allocation 
Leu                = work[0]
Lis                = work[1]
eta                = adimVec[2]
omega              = rotationVec[1] * 1.05
b0                 = bladeVec[0]
hubRadius          = rMean - b0/2
rotorVaMeanInlet   = V0vec[0]
rotorVtMeanInlet   = V0vec[1] 
rotorVtMeanOutlet  = Lis/0.82 / (rMean * omega) + rotorVtMeanInlet
rotorVaMeanOutlet  = V0vec[0]
#rotorVtMeanOutlet  = V1vec[1]
statorVaMeanInlet  = V0vec[0] 
statorVtMeanInlet  = V1vec[1]
statorVaMeanOutlet = V2vec[0]
statorVtMeanOutlet = V2vec[1]
Tt1                = thermo1[3]
Pt1                = thermo1[4]

# blade shape 
# rotor 
V1a = V0vec[0]
V1t = V0vec[1]
V1  = np.sqrt(V1a**2 + V1t**2)

W1a   = W0vec[0]
W1t   = W0vec[1]
W1    = np.sqrt(W1a**2 + W1t**2)
beta1 = W0vec[2]

V2a = V1vec[0]
V2t = V1vec[1]
V2  = np.sqrt(V2a**2 + V2t**2)

W2a   = W1vec[0]
W2t   = W1vec[1]
W2    = np.sqrt(W2a**2 + W2t**2)
beta2 = W1vec[2]
bladeHeight = b0

bladeStudy.optimalBladeNumber(W1, W2, beta1, beta2, rMean, bladeHeight, r1=rMean, r2=rMean, Vt1=V1t, Vt2=V2t, Va1=V1a, bladeInterval=[25,50], kind='std', save=True, position='latex/figures/rotorBlades.pdf', title='Rotor')

# stator
V3a = V1vec[0]
V3t = V1vec[1]
V3  = np.sqrt(V3a**2 + V3t**2)

W3a   = V1vec[0]
W3t   = V1vec[1]
W3    = np.sqrt(W3a**2 + W3t**2)
beta3 = V1vec[2]
bladeHeight = b0

V4a = V2vec[0]
V4t = V2vec[1]
V4  = np.sqrt(V4a**2 + V4t**2)

W4a   = V2vec[0]
W4t   = V2vec[1]
W4    = np.sqrt(W4a**2 + W4t**2)
beta4 = V2vec[2]
bladeHeight = b0

bladeStudy.optimalBladeNumber(W3, W4, beta3, beta4, rMean, bladeHeight, r1=rMean, r2=rMean, Vt1=V3t, Vt2=V4t, Va1=V3a, bladeInterval=[25,50], kind='std', save=True, position='latex/figures/statorBlades.pdf', title='Stator')