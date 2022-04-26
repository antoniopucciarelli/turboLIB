from turboCoeff import similarity

# data
# constraints
mFlux = 100     # mass flux                [kg/s]
betaP = 1.45    # compression ratio        [--]
maxD  = 0.9     # maximum tip diameter     [m]
maxR  = maxD/2  # maximum tip radius       [m]

# inlet values
Pt0 = 1e+5      # inlet total pressure     [Pa]
Tt0 = 300       # inlet total temperature  [K]

# choosed quantities
rD    = 0.58
rMean = 0.32

# mean line properties study 
V0vec, V1vec, V2vec, bladeVec, rotationVec = similarity.propertiesStudy(mFlux, betaP, Pt0, Tt0, input=[rD, rMean], rDmin=0.5, rDmax=0.73, rMeanMin=0.25, rMeanMax=0.4, Vt0Umean=0, plot=False, save=True, position='latex/rMeanChi.png')

# data allocation 
Va0 = V0vec[0]
Vt0 = V0vec[1]
Va1 = V1vec[0]
Vt1 = V1vec[1]
omega = rotationVec[1]
bladeHeight = bladeVec[0]
hubRadius = rMean - bladeHeight/2
rMean = [rMean, rMean]
VtMean = [Vt0, Vt1]
VaMean = [Va0, Va1]

WtTarget = 0
bVal = (WtTarget + omega * hubRadius - Vt1 * rMean[1]/hubRadius)/(hubRadius - rMean[1]**2/hubRadius)
b = [0, bVal]

def func(r):
    Vt2new = Vt1 * (rMean[1] / r)**-0.6

    return Vt2new

# angle and velocity study along the blade span
similarity.deltaAngleStudy(hubRadius, bladeHeight, rMean, VtMean, VaMean, omega, kind=['FV', 'MVD'], a=[0,0], b=b, n=[0,0], nSection=50, save=True, position='latex/figures/betaAngles.png')

# data allocation 
Va1 = V1vec[0]
Vt1 = V1vec[1]
Va2 = V2vec[0]
Vt2 = V2vec[1]
VtMean = [Vt1, Vt2]
VaMean = [Va1, Va2]

bVal = (omega - Vt1 * rMean[1]/hubRadius**2)/(1 - rMean[1]**2/hubRadius**2)
WtTarget = -omega*hubRadius*1.5*0
bVal1 = (WtTarget + omega * hubRadius - Vt1 * rMean[1]/hubRadius)/(hubRadius - rMean[1]**2/hubRadius)
b = [bVal1, 0]
omega = 0

VtTarget = 30
def func(r):
    a = 4 * VtTarget / bladeHeight**2
    Vt = a * (r - rMean[1])**2

    return Vt

# stator angles and velocity study 
similarity.deltaAngleStudy(hubRadius, bladeHeight, rMean, VtMean, VaMean, omega, kind=['MVD', 'eqn'], a=[0,0], b=b, n=[0,0], nSection=50, func2=func, save=True, position='latex/figures/alphaAngles.png')
