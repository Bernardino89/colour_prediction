#!/usr/bin/python


import numpy as np
import sys 
import math
import scipy.integrate as scint


    
filename=sys.argv[1]

#ine the parametres to find in a TD calculation and generate list of its

SCF = [] 			# SCF groud state
VEE = []			# Vertical Excitation Energy
L   = []			# wavelenght
OS  = []			# Oscillator str


#nerate lists

with open (filename, 'r') as data:
    
     for line in data:
         if 'SCF Done' in line:
              SCF.append(float(line.split()[4]))
#         if 'Excited State   1:' in line:
         if 'Excited State' in line:
             fragments = line.split()        		# split the line
             VEE.append(float(fragments[4]))
             L.append(float(fragments[6]))
             OS.append(float(fragments[8][2:]))

	 #print ('vertical excitation energie',VEE,'\n')
	 #print ('L',L, '\n')
	 #print ('oscillator',OS,'\n')


lmax = 1000			# wavelight maximum
lmin = 0.1			# wavelight minimum
FWHM = 2066.40322               # full width at half maximum 0.6 eV valore report
#FWHM = 2254.25805	        # full width at half maximum 0.55 eV
#FWHM = 2479.68386	        # full width at half maximum 0.5 eV 
#FWHM = 4132.80643	        # full width at half maximum 0.3 eV 
#FWHM = 3542.40551	        # full width at half maximum 0.35 eV 
#FWHM = 3099.60483	        # full width at half maximum 0.4 eV 
#FWHM = 6452
#FWHM = 6199.20965	        # full width at half maximum 0.2 eV 
#FWHM = 12398.41929	        # full width at half maximum 0.1 eV 
#FWHM = 8265.61287	        # full width at half maximum 0.15 eV 
#FWHM = 10000000.0
C = 2.99792458E+08		# m/s
h = 6.626068E-34		# j*s
CeVJ = 1.60217646E-19		# conversion from eV to J
cte = h*C/1.0e-09/CeVJ		# 1239.84171 nm is the wavelight of 1 eV photon
constant = 130629740.0		# constant in Gaussian pakege
step = 0.01
conv=10000000.0/FWHM		# conversion Gaussian pakege
recipFWHM=1.0/FWHM		# inverse of FWHM




# define the gaussian function and generate foreach EVV a gaussian function 


recipVEE = []
    
for g in VEE:
    nm= g/cte
    recipVEE.append(nm)

#print ('inverse VEE', recipVEE)


recipOS = []
for p in OS:
    con = p/conv
    #print("con:", con)
    recipOS.append(con)

#print ('inverse OS', recipOS)


# --- generate a file name comulative_spectra.txt 

WL = []
AB = []
WLEN = []


SQR = np.sqrt(4*np.log(2)/math.pi)
M   =  0.9
R = 1 - M
#a = - 4700
a = - 3000

GR = []
sigmoide = []


x_values = np.linspace(0.1,1000, 500)



ABG=[]
for wl in x_values:
   
    wlen = 1./wl
    absorbgaussian = 0.
    for i in range(len(VEE)):
   
        WE  = wlen - recipVEE[i]  
        sigm = (2  * recipFWHM)/(1 + np.exp( a * WE))        
        gr1 = 4 * np.log(2)* (WE/sigm)**2
              
        GF = np.exp(- gr1)                     
        absorbgaussian += (recipOS [i]) * GF 
    ABG.append(absorbgaussian)
    WL.append(wl)  

NNG = scint.simps(ABG,x = WL,dx = 0.01)
#NNG = 1
normg = []
for g in ABG:
       gau = R * 217500000 * g/NNG
       normg.append(gau)
        


ABL=[]
for wl in x_values:
   
    wlen = 1./wl
    absorbLF = 0.
    for i in range(len(VEE)):
   
        WE  = wlen - recipVEE[i]  
        sigm = (2  * recipFWHM)/(1 + np.exp( a * WE))        
        gr = (WE/sigm)**2
              
        LF = 1/(1 + 4 * gr)                     
        absorbLF += (recipOS [i]) * LF
         
    ABL.append(absorbLF)



NNL = scint.simps(ABL,x= WL, dx= 0.01)
#NNL = 1

norml = []
for l in ABL:
       lor = M * 217500000 * l/NNL
       norml.append(lor)

shape =[]
for i,l in zip(normg,norml):
       PV = i + l
       shape.append(PV)
 

#for j,i in zip(WL,shape):
#    
#    print ("{:.4f} {:4f}".format(j, i))


N_N = []
max_value = None
for num in shape:
    max_value = max(shape)
    new_norm = num/ max_value
    N_N.append(new_norm)

for j,i in zip(WL,N_N):
    
    print ("{:.8f} {:.8f}".format(j, i))
        
        









for wl in x_values:
   
    wlen = 1./wl

    for i in range(len(VEE)):
   
        WE  = wlen - recipVEE[i]  
        sigm = (2  * recipFWHM)/(1 + np.exp( a * WE))        
        gr1 = 4 * np.log(2)* (WE/sigm)**2
        gr = (WE/sigm)**2
        #gr = (WE/sigm)**2
        LF = 1/(1 + 4 * gr)
              
        GF = np.exp(- gr1)                     
       
        absorb0 = ( R * GF + M * LF )
        
        WL.append(wl)        
        AB.append(absorb0)
        
# Normalization factor (NNa)
        
#NN = scint.simps(AB,x = WL,dx = 0.01)

#F = []  
#for h in AB:
#    
#    f = (constant * recipOS[i]) * h /NN
#    F.append(f)
#
#N_N = []
#max_value = None
#for num in F:
#    max_value = max(F)
#    new_norm = num/ max_value
#    N_N.append(new_norm)

#for j,i in zip(WL,N_N):
#    
#    print ("{:.4f} {:4f}".format(j, i))
        
        

