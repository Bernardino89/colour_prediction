#!/usr/bin/python

"""

@author: bernardino.tirri@chimieparistech.psl.eu
"""


##################################################################################################
# 												 
#												 
# The program generates the spectra from Gaussian (.log) file. For each excitation energy the   #
# program generate a gauusan curve. At the end, gaussian_all_tran sum all gaussian in order to find the #
# spectra. gaussian_all_tran.py generates a file, spectra.txt, that contain the epsilon and the vawelight. #
#												 
#												 
#												 
# command line to use the program:								 
# 												 
#  				python.3.5 gaussian_all_tran file_name.log > test			 
# 												 
#												
##################################################################################################



import numpy as np
import os,sys,string
import pylab as pl
import matplotlib.pyplot as plt
import scipy as sy

filename=sys.argv[1]

# define the range of the spactra


#print ('==========================================================================================')

#print ('for UV-VIS spectra wavelenght minimum',)
#print ('for UV-VIS spectra wavelenght maximum',)

#print ('==========================================================================================')

# define the parametres to find in a TD calculation and generate list of its

SCF = [] 			# SCF groud state
VEE = []			# Vertical Excitation Energy
L   = []			# wavelenght
OS  = []			# Oscillator str

# generate lists

with open (filename, 'r') as data:

	for line in data:

		if 'SCF Done' in line:
			SCF.append(float(line.split()[4]))
#			print  ('SCF',SCF)

		if 'Excited State' in line:
			fragments = line.split()        		# split the line
			VEE.append(float(fragments[4]))
			L.append(float(fragments[6]))
			OS.append(float(fragments[8][2:]))

#	print ('vertical excitation energie',VEE,'\n')
#	print ('L',L, '\n')
#	print ('oscillator',OS,'\n')

# cycle to assicure that the lists have the same lenght 

#if len(VEE)==len(OS):

#	print ('VEE and OS lists have the same lenght')
#else:
#	print ('there is something wrong the lists do not have the same lenght')


#print ('I m calculating UV/Visible spectrum')


# define the costatnt to solve the gaussian equation


lmax = 1000			# wavelight maximum
lmin = 200			# wavelight minimum
#FWHM = 2066.40322               # full width at half maximum 0.6 eV valore report
#FWHM = 2479.68386	        # full width at half maximum 0.5 eV 
#FWHM = 4132.80643	        # full width at half maximum 0.3 eV 
#FWHM = 6199.20965	        # full width at half maximum 0.2 eV 
FWHM = 3099.60483	        # full width at half maximum 0.4 eV 
C = 2.99792458E+08		# m/s
h = 6.626068E-34		# j*s
CeVJ = 1.60217646E-19		# conversion from eV to J
cte = h*C/1.0e-09/CeVJ		# 1239.84171 nm is the wavelight of 1 eV photon
constant = 130629740.0		# constant in Gaussian pakege
step = 0.1			# step about wavelight
conv=10000000.0/FWHM		# conversion Gaussian pakege
recipFWHM=1.0/FWHM		# inverse of FWHM



# define the gaussian function and generate foreach EVV a gaussian function 


recipVEE = []
for g in VEE:
	nm = g/cte
	recipVEE.append(nm)

#print ('inverse VEE', recipVEE)


recipOS = []
for p in OS:
	con = p/conv
	recipOS.append(con)

#print ('inverse OS', recipOS)

# generate a file name comulative_spectra.txt 

WL = []
AB = []
with open("spectra.txt","w") as cmsp:

	for i in range(int((lmax-lmin)/step)+1):

		wl=lmin+i*step
		wlen=1./wl 

		absorb=0.
		for j in range(len(VEE)):
			#print ("recipOS [j]", recipOS [j])

			absorb += (recipOS [j]) * (217500000 * 1/recipFWHM) * np.exp(-4 * np.log(2)* (((wlen-recipVEE [j] )/recipFWHM )**2))

		#cmsp.write("%f %e \n" % (wl, absorb))

		WL.append(wl)
		AB. append(absorb)


		cmsp.write("%f %e \n" % (wl, absorb))

	#AB. append(absorb)
#print ('wavelight range', WL)
#print ('absorbance', AB)


N_N = []
max_value = None
for num in AB:
    max_value = max(AB)
    new_norm = num/ max_value
    N_N.append(new_norm)

for j,i in zip(WL,N_N):
    
    print ("{:.4f} {:4f}".format(j, i))



