import matplotlib.pyplot as plt 
import numpy as np 
import warnings
from turboMachine import * 

naca1212 = "../data/airfoils/naca1212"
naca5016 = "../data/airfoils/naca5016"

stage1 = stage(1, 'turbine', rotor = rotor(1), stator = stator(1))
stage1.stator.setFoil(naca1212)
stage1.rotor.setFoil(naca5016)
stage1.rotor.setAngle(10)
stage1.rotor.foilRotation()

turbine = turbomachine(1, 2, 'turbine')

# stage 1 analysis
# rotor allocation 
turbine.stage[0].rotor.setFoil(naca1212)
turbine.stage[0].rotor.setAngle(-30)
turbine.stage[0].rotor.foilRotation()
turbine.stage[0].rotor.chord = 1
# stator allocation 
turbine.stage[0].stator.setFoil(naca5016)
turbine.stage[0].stator.setAngle(10)
turbine.stage[0].stator.foilRotation()
turbine.stage[0].stator.chord = 1
# intrastage allocation 
turbine.stage[0].intrastageGap = 0.3
turbine.stage[0].interstageGap = 0.5

# stage 2 analysis
# rotor allocation
turbine.stage[1].rotor.setFoil(naca1212)
turbine.stage[1].rotor.setAngle(30)
turbine.stage[1].rotor.foilRotation()
turbine.stage[1].rotor.chord = 1
# stator allocation 
turbine.stage[1].stator.setFoil(naca5016)
turbine.stage[1].stator.setAngle(-10)
turbine.stage[1].stator.foilRotation()
turbine.stage[1].stator.chord = 1
# intrastage allocation 
turbine.stage[1].intrastageGap = 0.3

turbine.print()

turbine.turboPlot()

turbine.
