# TURBOMACHINERY -- ENGINEERING DIMENSION AND COEFFICIENT FUNCTIONS LIBRARY
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY: compressor
#   CONTENT: engineering coefficient functions for the computation of turbomachinery performace
#  

# importin libraries
import numpy as np 
import matplotlib.pyplot as plt 

def efficiency(filePath,plot=False):
    '''
    This function describes the adimensional parameters for the turbomachinery design.
        !!! data are taken from BASKHARONE book -- fig 5.8 !!! 
    '''

    # importing data 
    radialEfficiency = np.float32(np.loadtxt(filePath + 'data/perfCoeff/radialEfficiency.txt'))
    axialEfficiency  = np.float32(np.loadtxt(filePath + 'data/perfCoeff/axialEfficiency.txt'))

    # plotting efficiency
    if plot:
        plt.figure(figsize=(8,8))
        plt.title('Efficiency from Baskharone')
        plt.semilogx(radialEfficiency[:,0],radialEfficiency[:,1],linewidth=3,label='radial compressor')
        plt.semilogx(axialEfficiency[:,0],axialEfficiency[:,1],linewidth=3,label='axial compressor')
        plt.xlabel('specific speed')
        plt.ylabel('efficiency')
        plt.xticks([0.1,1,10])
        plt.legend()  
        plt.grid(which='both')
        plt.show()