from geometry.bladeGenerator import *

# setting up path for NACA65 coordinates
pos = 'data/airfoils/naca65.txt' 

# HUB airfoil setup 
# getting airfoil data 
naca65_hub = geometryData(pos)

# airfoil shape 
naca65_hub.geometryFitting(Cl=2.586, chord=1, plot=False)

# airfoil 3D rotation
naca65_hub.geometryRotation(8.581, 10, plot=False)

# airfoil 3D translation 
#naca65_hub.geometryTranslation(naca65_hub.middleChord(), 0.20, plot=False)

# MIDSPAN airfoil setup 
# getting airfoil data
naca65_mid = geometryData(pos)

# airfoil shape 
naca65_mid.geometryFitting(Cl=1.031, chord=1, plot=False)

# airfoil 3D rotation
naca65_mid.geometryRotation(39.58, 0, plot=False)

# airfoil 3D translation 
naca65_mid.geometryTranslation(naca65_hub.middleChord(), (0.325-0.201)*1e+2, plot=False)

# TIP airfoil setup 
# getting airfoil data
naca65_tip = geometryData(pos)

# airfoil shape 
naca65_tip.geometryFitting(Cl=0.423, chord=1, plot=False)

# airfoil 3D rotation
naca65_tip.geometryRotation(55.285, 0, plot=False)

# airfoil 3D translation 
naca65_tip.geometryTranslation(naca65_hub.middleChord(), (0.4491-0.201)*1e+2, plot=False)

# stl generation 
airfoils = [naca65_hub, naca65_mid, naca65_tip]
STLsaving(airfoils, STLname='cad')