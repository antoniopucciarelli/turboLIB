#----------------------------------------------------------------------------------------------------------------------#
# Create the 2D section -- STATOR
#----------------------------------------------------------------------------------------------------------------------#
# variable allocation
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

# stator object setup
u = np.linspace(0.00, 1.00, 1000)
stator10 = Blade2DCamberThickness(design_variables)
stator10.get_section_coordinates(u)
stator10.check_analytic_curvature()

# openscad file generation
upper  = np.real(rotor10.get_upper_side_coordinates(u))
lower  = np.real(rotor10.get_lower_side_coordinates(u))
height = 0 

# saving data
file = open('statorCoords.scad', 'w')
file.write('stator0 = [\n')
for ii in range(upper.shape[1]):
    file.write('[{0:f}, {1:f}],\n'.format(upper[0,ii], upper[1,ii]))
for ii in range(lower.shape[1]): 
    file.write('[{0:f}, {1:f}],\n'.format(lower[0,ii], lower[1,ii]))
file.write('];\n')
file.close()

#----------------------------------------------------------------------------------------------------------------------#
# Create the 2D section -- STATOR
#----------------------------------------------------------------------------------------------------------------------#
# variables allocation
design_variables = {}
design_variables['stagger'] = -45 * np.pi / 180
design_variables['theta_in'] = +30 * np.pi / 180
design_variables['theta_out'] = -70 * np.pi / 180
design_variables['radius_in'] = 0.05
design_variables['radius_out'] = 0.01
design_variables['dist_in'] = 0.60
design_variables['dist_out'] = 0.5
design_variables['thickness_upper_1'] = 0.15
design_variables['thickness_upper_2'] = 0.10
design_variables['thickness_upper_3'] = 0.13
design_variables['thickness_upper_4'] = 0.07
design_variables['thickness_upper_5'] = 0.03
design_variables['thickness_upper_6'] = 0.02
design_variables['thickness_lower_1'] = 0.10
design_variables['thickness_lower_2'] = 0.15
design_variables['thickness_lower_3'] = 0.1
design_variables['thickness_lower_4'] = 0.07
design_variables['thickness_lower_5'] = 0.03
design_variables['thickness_lower_6'] = 0.02

# Convert standard-python scalars into singleton numpy arrays
for i in design_variables:
    design_variables[i] = np.asarray(design_variables[i])

# rotor object allocation 
u = np.linspace(0.00, 1.00, 1000)
stator11 = Blade2DCamberThickness(design_variables)
stator11.get_section_coordinates(u)
stator11.check_analytic_curvature()

# openscad file generation
upper  = np.real(rotor11.get_upper_side_coordinates(u))
lower  = np.real(rotor11.get_lower_side_coordinates(u))
height = 0 

file = open('statorCoords.scad', 'a')
file.write('\nstator1 = [\n')
for ii in range(upper.shape[1]):
    file.write('[{0:f}, {1:f}],\n'.format(upper[0,ii], upper[1,ii]))
for ii in range(lower.shape[1]): 
    file.write('[{0:f}, {1:f}],\n'.format(lower[0,ii], lower[1,ii]))
file.write('];\n')
file.close()
