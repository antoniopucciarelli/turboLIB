#!/usr/bin/python3

# importing libraries
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# importing user defined packages 
sys.path.append(os.getcwd() + '/../parablade')
from parablade.blade_2D_camber_thickness import Blade2DCamberThickness

# 
# ROTOR GENERATION 
#
dim = 1000
# design variables allocation
design_variables = {}
design_variables['stagger'] = -45 * np.pi / 180
design_variables['theta_in'] = +30 * np.pi / 180
design_variables['theta_out'] = -70 * np.pi / 180
design_variables['radius_in'] = 0.05
design_variables['radius_out'] = 0.01
design_variables['dist_in'] = 0.60
design_variables['dist_out'] = 0.5
design_variables['thickness_upper_1'] = 0.15
design_variables['thickness_upper_2'] = 0.20
design_variables['thickness_upper_3'] = 0.13
design_variables['thickness_upper_4'] = 0.07
design_variables['thickness_upper_5'] = 0.03
design_variables['thickness_upper_6'] = 0.02
design_variables['thickness_lower_1'] = 0.15
design_variables['thickness_lower_2'] = 0.20
design_variables['thickness_lower_3'] = 0.13
design_variables['thickness_lower_4'] = 0.07
design_variables['thickness_lower_5'] = 0.03
design_variables['thickness_lower_6'] = 0.02

# convert standard-python scalars into singleton numpy arrays
for i in design_variables:
    design_variables[i] = np.asarray(design_variables[i])

# rotor object setup
u      = np.linspace(0.00, 1.00, dim)
rotor0 = Blade2DCamberThickness(design_variables)
rotor0.get_section_coordinates(u)
rotor0.check_analytic_curvature()

# upper and lower part of the 2D airfoil generation 
upper  = np.real(rotor0.get_upper_side_coordinates(u))
lower  = np.real(rotor0.get_lower_side_coordinates(u))
height = 0 

# rotor array declaration 
rotor0 = np.zeros(shape = (upper.shape[1]+lower.shape[1]-1, 3))
for ii in range(upper.shape[1]):
    rotor0[ii,0] = upper[0, ii]
    rotor0[ii,1] = upper[1, ii]
    rotor0[ii,2] = height

for ii in range(lower.shape[1]-1):
    rotor0[upper.shape[1]+ii, 0] = lower[0, ii+1]
    rotor0[upper.shape[1]+ii, 1] = lower[1, ii+1]
    rotor0[upper.shape[1]+ii, 2] = height

# 
# ROTOR GENERATION 
#

# design variables allocation
design_variables = {}
design_variables['stagger'] = -45 * np.pi / 180
design_variables['theta_in'] = +30 * np.pi / 180
design_variables['theta_out'] = -70 * np.pi / 180
design_variables['radius_in'] = 0.05
design_variables['radius_out'] = 0.01
design_variables['dist_in'] = 0.60
design_variables['dist_out'] = 0.5
design_variables['thickness_upper_1'] = 0.15
design_variables['thickness_upper_2'] = 0.20
design_variables['thickness_upper_3'] = 0.13
design_variables['thickness_upper_4'] = 0.07
design_variables['thickness_upper_5'] = 0.03
design_variables['thickness_upper_6'] = 0.02
design_variables['thickness_lower_1'] = 0.15
design_variables['thickness_lower_2'] = 0.20
design_variables['thickness_lower_3'] = 0.13
design_variables['thickness_lower_4'] = 0.07
design_variables['thickness_lower_5'] = 0.03
design_variables['thickness_lower_6'] = 0.02

# convert standard-python scalars into singleton numpy arrays
for i in design_variables:
    design_variables[i] = np.asarray(design_variables[i])

# rotor object setup
u       = np.linspace(0.00, 1.00, dim)
rotor1 = Blade2DCamberThickness(design_variables)
rotor1.get_section_coordinates(u)
rotor1.check_analytic_curvature()

# upper and lower part of the 2D airfoil generation 
upper  = np.real(rotor1.get_upper_side_coordinates(u))
lower  = np.real(rotor1.get_lower_side_coordinates(u))
height = 1

# translation vector 
translate = [0, 0, 0]
# rotation angle in degrees
alpha = np.deg2rad(0)
rotMatr = np.matrix([[np.cos(alpha), -np.sin(alpha)],[np.sin(alpha), np.cos(alpha)]])

# rotor array declaration 
rotor1 = np.zeros(shape = (upper.shape[1]+lower.shape[1]-1, 3))
for ii in range(upper.shape[1]):
    upper[:,ii] = np.matmul(rotMatr, upper[:,ii])
    rotor1[ii,0] = upper[0, ii] + translate[0]
    rotor1[ii,1] = upper[1, ii] + translate[1]
    rotor1[ii,2] = height       + translate[2]

for ii in range(lower.shape[1]-1):
    lower[:,ii+1] = np.matmul(rotMatr, lower[:,ii+1])
    rotor1[lower.shape[1]+ii, 0] = lower[0, ii+1] + translate[0]
    rotor1[lower.shape[1]+ii, 1] = lower[1, ii+1] + translate[1]
    rotor1[lower.shape[1]+ii, 2] = height         + translate[2]

# saving data
file = open('rotorCoords.scad', 'w')
file.write('rotor = [\n')
# rotor0 saving
for ii in range(rotor0.shape[0]):
    file.write('[{0:f}, {1:f}, {2:f}],\n'.format(rotor0[ii,0], rotor0[ii,1], rotor0[ii,2]))
# rotor1 saving 
for ii in range(rotor1.shape[0]):
    file.write('[{0:f}, {1:f}, {2:f}],\n'.format(rotor1[ii,0], rotor1[ii,1], rotor1[ii,2]))
file.write('];\n')
file.close()
    
file = open('cad.stl', 'w')
file.write('solid blade\n')
for ii in range(rotor0.shape[0]-1):
    # versor computation
    vec1 = np.array(rotor0[ii,:])
    vec2 = np.array(rotor0[ii+1,:])
    vec3 = np.array(rotor1[ii,:])
    versor = np.cross(vec3 - vec1, vec2 - vec1)
    # writing data
    versor = versor / np.linalg.norm(versor)
    file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
    file.write('\t\touter loop\n')
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
    file.write('\t\tendloop\n')
    file.write('\tendfacet\n')

    # versor computation
    vec1 = np.array(rotor1[ii,:])
    vec2 = np.array(rotor1[ii+1,:])
    vec3 = np.array(rotor0[ii+1,:])
    versor = - np.cross(vec3 - vec1, vec2 - vec1)
    # writing data
    versor = versor / np.linalg.norm(versor)
    file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
    file.write('\t\touter loop\n')
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
    file.write('\t\tendloop\n')
    file.write('\tendfacet\n')

for ii in range(rotor0.shape[0]-1):
    # versor computation 
    vec1 = np.array([0.3,0,0])
    vec2 = np.array(rotor0[ii,:])
    vec3 = np.array(rotor0[ii+1,:])
    # versor computation
    versor = np.cross(vec3 - vec1, vec2 - vec1)
    # writing data
    versor = versor / np.linalg.norm(versor)
    file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
    file.write('\t\touter loop\n')
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
    file.write('\t\tendloop\n')
    file.write('\tendfacet\n')

for ii in range(rotor1.shape[0]-1):
    # versor computation 
    vec1 = np.array([0.3,0,1])
    vec2 = np.array(rotor1[ii,:])
    vec3 = np.array(rotor1[ii+1,:])
    # computing versor 
    versor = - np.cross(vec3 - vec1, vec2 - vec1)
    # writing data
    versor = versor / np.linalg.norm(versor)  
    file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
    file.write('\t\touter loop\n')
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
    file.write('\t\tendloop\n')
    file.write('\tendfacet\n')

file.write('endsolid blade\n')
file.close()

