#!/usr/bin/env python3

import sys, getopt
import click
from colorama import Fore,Style
import os.path
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import re 



def import_bar(bar,nsim):
	ddG = []
	with open(bar,'r') as file:
		i = 0
		for line in file:
			if 'Final results in kJ/mol:' in line:
				i = 1
			elif 'total' in line:
				i=2
				break
			elif i==1:
				ddG.append(re.sub("\s\s+", " ", line).split(' ')) # remove whitespaces and split to array
		ddG=[item for sublist in ddG for item in sublist] # flatten array
		ddG=ddG[1:-1] # remove star/end "\n"
		ddG_y=np.asarray(ddG[5::8],dtype='float64') # import every 8th element
		ddG_x=np.asarray(ddG[1::8],dtype='float64')
#		ddG_y_interp=np.interp(ddG_x,ddG_x,ddG_y) # 1D interpolation
		ddG_y_interp = interp1d(ddG_x,ddG_y,kind='cubic')
		ddG_x_interp = np.linspace(0, nsim-1, num=25*nsim, endpoint=True)
		print(len(ddG_x),len(ddG_y),len(ddG_y_interp(ddG_x)),len(ddG_x_interp))
		plt.plot(ddG_x,ddG_y,'o')
		plt.plot(ddG_x_interp, ddG_y_interp(ddG_x_interp), '-x')
		plt.legend(['data', 'cubic interpolation'], loc='best')
		plt.xlabel('N') 
		plt.ylabel('ddG')
		plt.savefig('ddG_interpolation.pdf')
		# nesto nije uredu, visak elemenata u equi_lambdas
		#mean_ddG = np.mean(ddG_y_interp(ddG_x_interp))
		#std_ddG = np.std(ddG_y_interp(ddG_x_interp))
		equi_ddG = np.abs(ddG_y_interp(ddG_x_interp)[0]-ddG_y_interp(ddG_x_interp)[-1])/nsim
		#equi_ddGs = [k*equi_ddG for k in range(nsim)]
		equi_ddGs = []
		equi_lambdas = []
		for k in ddG_x_interp:
			min_ddG = ddG_y_interp(ddG_x_interp)[0]
			if ddG_y_interp(k)>=min_ddG:
				equi_ddGs.append(ddG_y_interp(k))
				equi_lambdas.append("{0:.4f}".format(k))
				min_ddG=min_ddG+equi_ddG
		equi_lambdas[-1] = '1.0000'
		print(equi_ddGs,np.std(equi_ddGs))
		print(equi_lambdas)

		print(equi_ddG)
		#print(equi_ddGs)
		#print(ddG_y_interp(ddG_x_interp))
		#print(ddG_y_interp(ddG_x))
		#print(ddG_y_interp(ddG_x_interp))



def create_lambdas_equiG(nsim,start,end,bar_array):
	# Creates lambda values equidistant with respect to DDG
	
	return None

#nsim=int(input("Input number of sims: \n >>>"))
nsim=20
import_bar('./bar_results.txt',nsim)