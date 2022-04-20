# importing libraries 
from numpy import False_
from geometry import bladeGenerator
from turboCoeff import similarity
from turboClass import bladeStudy

# data
# constraints
mFlux = 100     # mass flux                [kg/s]
betaP = 1.45    # compression ratio        [--]
maxD  = 0.9     # maximum tip diameter     [m]
maxR  = maxD/2  # maximum tip radius       [m]

# inlet values
Pt0 = 1e+5      # inlet total pressure     [Pa]
Tt0 = 300       # inlet total temperature  [K]

# air properties allocation 
R     = 287.06                  # air gas constant Ru / Mm      [J/kg K]
gamma = 1.4                     # specific heat ratio           [--]
cP    = gamma / (gamma - 1) * R # specific heat ratio @ P cost  [J/kg K]
cV    = cP - R                  # specific heat ratio @ V cost  [J/kg K]

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

# plotting charts 
phi, psi = similarity.stagePerf(psi=psiTarget, rD=rD, perc=0.98, save=True, position='latex/figures/stagePerf.png')
eta = similarity.efficiency(phi=phi, rD=rD, plot=False, save=True, position='latex/figures/efficiency.png')
#similarity.stageStudy(mFlux, betaP, rMean, Pt0, Tt0, rDmin=0.5, rDmax=0.75, Vt0UmeanMin=0, Vt0UmeanMax=0.25, R=287.06, gamma=1.4)
#similarity.reactionStudy(mFlux, betaP, rMean, Pt0, Tt0, rDmin=0.5, rDmax=0.73, Vt0UmeanMin=0, Vt0UmeanMax=0.1, save=True, position0='latex/figures/reactionStudy0.png', position1='latex/figures/reactionStudy1.png')
#bladeStudy.optimalPlot(beta1, beta2, error=1e-4)

exit()

# generation of mean line properties to be used for the blade assembly 
adimVec, bladeVec, rotationVec, V0vec, V1vec, V2vec, _, _, _, thermo0, _, _, work = similarity.stageProperties(rD, psiTarget, rMean, mFlux, Tt0, Pt0, betaP, T1real=False, R=R, gamma=gamma)

# ROTOR BLADE ASSEMBLY 
# values allocation 
meanValues   = [rMean, rotationVec[0], V0vec[0], V0vec[1], V1vec[1]]
b0           = bladeVec[0]
b1           = bladeVec[1]
Leu          = work[0]
inletValues  = [thermo0[3], thermo0[4], thermo0[0], thermo0[1]]
nSections    = 50
nRotorBlades = 40

# rotor blade assembly 
rotorHub, rotorTip =  bladeStudy.bladeGenerator('rotor', meanValues, b0, b1, Leu, inletValues, nSections, STLname='rotor', etaVec=1, hubChord=0, nBlades=nRotorBlades, pos='data/airfoils/naca65.txt')

# STATOR BLADE ASSEMBLY 
# values allocation 
meanValues    = [rMean, 0, V0vec[0], V1vec[1], V2vec[1]]
b1            = bladeVec[1]
b2            = bladeVec[2]
Leu           = 0
inletValues   = [thermo0[3], thermo0[4], thermo0[0], thermo0[1]]
nSections     = 50
nStatorBlades = 40

# stator blade assembly 
statorHub, statorTip, = bladeStudy.bladeGenerator('stator', meanValues, b1, b2, Leu, inletValues, nSections, STLname='stator', etaVec=1, hubChord=0, nBlades=nStatorBlades, pos='data/airfoils/naca65.txt')

# saving data 
bladeGenerator.SCADsaving(nRotorBlades, nStatorBlades, rotorHub, statorHub, rMean, b0, b1, rotorPath='../container/rotor.stl', statorPath='../container/stator.stl', geometryPath='geometry/')


