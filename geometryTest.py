from geometry.bladeGenerator import *

# setting up path for NACA65 coordinates
pos = 'data/airfoils/naca65.txt' 

# HUB airfoil setup 
# getting airfoil data 
naca65_hub = geometryData(pos)

# airfoil shape 
naca65_hub.geometryFitting(Cl=1.4, chord=1, plot=False)

# airfoil 3D rotation
naca65_hub.geometryRotation(0, -5, plot=False)

# MIDSPAN airfoil setup 
# getting airfoil data
naca65_mid = geometryData(pos)

# airfoil shape 
naca65_mid.geometryFitting(Cl=1.15, chord=0.8, plot=False)

# airfoil 3D rotation
naca65_mid.geometryRotation(5, 0, plot=False)

# airfoil 3D translation 
naca65_mid.geometryTranslation(naca65_hub.middleChord(), 0.35, plot=False)

# TIP airfoil setup 
# getting airfoil data
naca65_tip = geometryData(pos)

# airfoil shape 
naca65_tip.geometryFitting(Cl=1.1, chord=0.6, plot=False)

# airfoil 3D rotation
naca65_tip.geometryRotation(10, 5, plot=False)

# airfoil 3D translation 
naca65_tip.geometryTranslation(naca65_hub.middleChord(), 0.7, plot=False)

# stl generation 
airfoils = [naca65_hub, naca65_mid, naca65_tip]
STLsaving(airfoils, STLname='cad')