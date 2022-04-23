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

similarity.propertiesStudy(mFlux, betaP, Pt0, Tt0, rDmin=0.5, rDmax=0.73, rMeanMin=0.25, rMeanMax=0.4, Vt0UmeanMin=0, Vt0UmeanMax=0.25, save=False, position0='reactionStudy0.pgf', position1='reactionStudy1.pgf', gamma=1.4, R=287.06)