#!/usr/bin/env python3

import sys, getopt
import click
from colorama import Fore,Style
import os.path
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import re 


"""

1) DONE import bar_results file and extract ddG for nsims from bar_results into an array 
2) import Production MDP file and read lambda ranges for each lambda type into Ntype-dim array
		- exclude inactive lambda types
3) specify total nsim and divide nsim into nsubsim for each lambda type [in array of 2)]
*) check if lambdas overlap
4) create equidistant lambdas for each lambda type  (asuming no overlap between lambdas)
		- interpolate ddG array and find stationary points then create [ddG_x,ddG_y] arrays on intervals inbetween
		- compute optimal lambdas for each interval 
			- DONE for one interval
5) forward new lambdas to gen_barinp.py
6) clean up code

*) if lambdas overlap implement multivariate interpolation in 4)
	- requires additional simulation (impossible to deconvolve ddG)
	- not worth implementing for now

"""

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
		ddG = [item for sublist in ddG for item in sublist] # flatten array
		ddG = ddG[1:-1] # remove start/end "\n"
		#ddG_x = np.asarray(ddG[1::8],dtype='float64') # import every 8th element --> gives sim N
		ddG_y = np.asarray(ddG[5::8],dtype='float64')
		ddG_x = np.linspace(0,1,num=len(ddG_y))

		ddG_x_interp = np.linspace(0, 1, num=25*nsim, endpoint=True)
		ddG_y_interp_func = interp1d(ddG_x,ddG_y,kind='cubic')
		ddG_y_interp = ddG_y_interp_func(ddG_x_interp)

		ddG_y_first = ddG_y_interp[0]
		ddG_y_last = ddG_y_interp[-1]
		ddG_y_min = np.amin(ddG_y_interp_func(ddG_x))
		ddG_y_max = np.amax(ddG_y_interp_func(ddG_x))
		#equi_ddG = np.abs(ddG_y_last-ddG_y_first)/nsim
		equi_ddG = np.abs(ddG_y_max-ddG_y_min)/nsim

		ddG_y2 = np.asarray([ddG_y_min + i*equi_ddG for i in range(nsim-1)])
		ddG_x2_func = interp1d(ddG_y_interp,ddG_x_interp,kind='linear')
		ddG_x2 = ddG_x2_func(ddG_y2)

		# append missing lambda up to 1
		ddG_x2 = np.append(ddG_x2,0.0)
		ddG_y2 = np.append(ddG_y2,ddG_y[0])
		#ddG_x2 = np.insert(ddG_x2,0,1.0)
		#ddG_y2 = np.insert(ddG_y2,0,ddG_y[-1])
		#print(np.insert(ddG_y2,1,ddG_y[-1]))
		#print(ddG_x)
		#print(ddG_y)
		#print(ddG_y2)
		print(ddG_x2[::-1]) # print lambdas in reverse order
		#print("_lambdas = "," ".join(list(map(str,ddG_x2[::-1].tolist()))))
		print("len(ddG_x) = ",len(ddG_x),"\nlen(ddG_y) = ",len(ddG_y),"\nlen(ddG_y_interp) = ",len(ddG_y_interp),"\nlen(ddG_x2) = ",len(ddG_x2),"\nlen(ddG_y2) = ",len(ddG_y2))
		

		# Plot results
		plt.rc('text', usetex=True)
		plt.plot(ddG_x,ddG_y,'*',markersize=5,fillstyle='none')
		plt.plot(ddG_x2,ddG_y2,'o',markersize=5,fillstyle='none')
		plt.plot(ddG_x_interp, ddG_y_interp, '-')
		plt.legend(['input lambdas', 'cub/lin interpolated lambdas','cubic interpolation func'], loc='best')
		plt.xlabel(r'$\lambda$') 
		plt.ylabel(r"$\Delta \Delta G / \mathrm{kJ mol^{-1}}$")
		plt.savefig('ddG_interpolation.pdf')

		return ddG_x2



def create_lambdas_equiG(nsim,start,end,bar_array):
	# Creates lambda values equidistant with respect to DDG
	
	return None

if __name__ == '__main__':
	#nsim=int(input("Input number of sims: \n >>>"))
	nsim=40
	import_bar('./bar_results.txt',nsim)