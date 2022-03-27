#!/usr/bin/python3

""" Demonstration of the 2D blade parametrization based on camberline and thickness """

#----------------------------------------------------------------------------------------------------------------------#
# Importing general packages
#----------------------------------------------------------------------------------------------------------------------#
import sys
import os
import numpy as np
import matplotlib.pyplot as plt


#----------------------------------------------------------------------------------------------------------------------#
# Importing user-defined packages
#----------------------------------------------------------------------------------------------------------------------#
sys.path.append(os.getcwd() + '/../parablade')
from parablade.blade_2D_camber_thickness import Blade2DCamberThickness
import parablade

#----------------------------------------------------------------------------------------------------------------------#
# Define the section design variables
#----------------------------------------------------------------------------------------------------------------------#
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

# Convert standard-python scalars into singleton numpy arrays
for i in design_variables:
    design_variables[i] = np.asarray(design_variables[i])


#----------------------------------------------------------------------------------------------------------------------#
# Create the 2D section
#----------------------------------------------------------------------------------------------------------------------#
u = np.linspace(0.00, 1.00, 1000)
my_blade = Blade2DCamberThickness(design_variables)
my_blade.get_section_coordinates(u)
my_blade.check_analytic_curvature()

# openscad file generation
upper  = np.real(my_blade.get_upper_side_coordinates(u))
lower  = np.real(my_blade.get_lower_side_coordinates(u))
height = 0 

file = open('rotorCoords.scad', 'w')
file.write('rotor1 = [\n')
for ii in range(upper.shape[1]):
    #file.write('[{0:f}, {1:f}, {2:f}],\n'.format(upper[0,ii], upper[1,ii], height))
    file.write('[{0:f}, {1:f}],\n'.format(upper[0,ii], upper[1,ii]))
for ii in range(lower.shape[1]): 
    #file.write('[{0:f}, {1:f}, {2:f}],\n'.format(lower[0,ii], lower[1,ii], height))
    file.write('[{0:f}, {1:f}],\n'.format(lower[0,ii], lower[1,ii]))
file.write('];\n')
file.close()

file = open('statorCoords.scad', 'w')
file.write('stator1 = [\n')
for ii in range(upper.shape[1]):
    #file.write('[{0:f}, {1:f}, {2:f}],\n'.format(upper[0,ii], upper[1,ii], height))
    file.write('[{0:f}, {1:f}],\n'.format(upper[0,ii], upper[1,ii]))
for ii in range(lower.shape[1]): 
    #file.write('[{0:f}, {1:f}, {2:f}],\n'.format(lower[0,ii], lower[1,ii], height))
    file.write('[{0:f}, {1:f}],\n'.format(lower[0,ii], lower[1,ii]))
file.write('];\n')
file.close()


#----------------------------------------------------------------------------------------------------------------------#
# Plot the 2D section
#----------------------------------------------------------------------------------------------------------------------#
# Plot a single blade
#my_blade.plot_blade_section(upper_side='yes', upper_side_control_points='yes',
#                            lower_side='yes', lower_side_control_points='yes',
#                            camberline='yes', camberline_control_points='yes', camberline_sample_points='yes',
#                            leading_edge_radius='yes', trailing_edge_radius='yes')

# Plot a cascade of blades
#my_blade.plot_blade_cascade()

# Plot the thickness distribution
#my_blade.plot_thickness_distribution()

# Plot the curvature distribution
#my_blade.plot_curvature_distribution()


# Show the figure
#plt.show()
