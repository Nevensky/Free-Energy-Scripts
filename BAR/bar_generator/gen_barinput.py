#!/usr/bin/env python3

# skripta za generiranje Free Energy BAR inputa

import sys, getopt
import click
import os
import numpy as np
from scipy.interpolate import interp1d
#from numpy import interp as interp1d
import matplotlib.pyplot as plt
import re 
from termcolor import colored

import ddGintervals
import importbar as bar

unil='\u03BB'

root = os.getcwd()

@click.command()

@click.option('-nsim',default=20,help='Total number of simulations.')
@click.option('-vdw',default=False,help='Mutate van der Waals interactions.')
@click.option('-coul',default=False,help='Mutate Coulomb interactions.')
@click.option('-mass',default=False,help='Mutate mass of species.')
@click.option('-bonded',default=False,help='Mutate bonded interactions.')
@click.option('-restraint',default=False,help='Mutate restrained interactions.')
@click.option('-temp',default=False,help='Simulated temptering.')

@click.option('-root',default=root,help='Total number of simulations.')
@click.option('-mdp',default=root+'/MDP',help='MDP folder.')
@click.option('-min1',default=root+'/MDP/EM/em_steep.mdp',help='MDP of steepest decent minimization.')
@click.option('-min2',default=root+'/MDP/EM/em_l-bfgs.mdp',help='MDP of L-BFGS minimization.')
@click.option('-nvt',default=root+'/MDP/NVT/nvt.mdp',help='MDP of NVT equilibration.')
@click.option('-npt',default=root+'./MDP/NPT/npt.mdp',help='MDP of NPT equilibration.')
@click.option('-prod',default=root+'/MDP/Production_MD/md.mdp',help='MDP of the production run.')

@click.option('-ncores',default=4,help='Number of cores per simulation.')
@click.option('-pin',default='auto',help='Pin to cores.')

def inputs(nsim,vdw,coul,mass,bonded,restraint,temp,root,mdp,min1,min2,nvt,npt,prod,ncores,pin):
	"""Free Energy BAR input generator. \n

Developed by Neven Golenic | neven.golenic@gmail.com

\b
Number of simulations [ N ] default:  20 
\u03BB coupling interval [N_min,N_max] default: 0,0

If BAR results from a previous simulation are imported, lambdas of different type must not overlap!

\b
Recommended inital folder structure:
.
└── analysis
|
└── MDP
|   |
|   ├─ EM
|   |  ├── em_steep.mdp
|   |  └── em_l-bfgs.mdp
|   ├─ NVT
|   |   └── nvt.mdp 
|   ├─ NPT
|   |   └── npt.mdp 
|   └─ Production_MD
|       └── md.mdp
|
└── STRUCT
      |
      ├─ system.top
      └─ system.gro
	"""
	check_analysis(root)
	check_mdp(mdp,min1,min2,nvt,npt,prod) # ovo je dobro treba reaktivirati
	coupling(nsim,vdw,coul,mass,bonded,restraint,temp)

	return None

def check_analysis(root):
	if os.path.isdir(root+"/analysis"):
		pass
	else:
		os.makedirs(root+"/analysis")

def check_mdp(mdp,min1,min2,nvt,npt,prod):
	""" 
	Check if MDP file exists or if all 
	.mdp files were specified manually 
	"""
	if os.path.exists(min1) and os.path.exists(min2) and os.path.exists(nvt) and os.path.exists(npt) and os.path.exists(prod):
		pass	
	elif os.path.isdir(mdp) :
		if os.path.exists(mdp+'/EM/em_steep.mdp') and os.path.exists(mdp+'/EM/em_l-bfgs.mdp') and os.path.exists(mdp+'/NVT/nvt.mdp')==True and os.path.exists(mdp+'/NPT/npt.mdp') and os.path.exists(mdp+'/Production_MD/md.mdp'):
			pass
		else: 
			print(colored("""ERROR: MDP folder corrupt. Expected folder structure:\n
				└── analysis
				|
				└── MDP
				|   |
				|   ├─ EM
				|   |  ├── em_steep.mdp
				|   |  └── em_l-bfgs.mdp
				|   ├─ NVT
				|   |   └── nvt.mdp 
				|   ├─ NPT
				|   |   └── npt.mdp 
				|   └─ Production_MD
				|       └── md.mdp
				|
				└── STRUCT
				      |
				      ├─ system.top
				      └─ system.gro
				""","red",attrs=["bold"]))
			exit()
	else:
		print(colored("ERROR: MDP folder not found, please specify correct MDP file paths manually.","red",attrs=["bold"]))
		auto_initialize = input("Intialize MDP files automatically? yes/no :")
		if auto_initialize=="yes" or auto_initialize=="y":
			print("Currently not implemented. You must intialize MDP files manually.")
			exit()
		else:
			exit()

# def import_mdp_old(min1,min2,nvt,npt,prod,vdw_lambdas):
# #   pass # jer je ova funkcija losa
# 	# Imports mdp files and replaces them with new ones with correct lambdas
# 	with open(min1) as file:
# 		filedata = open('MDP/em_steep.mdp','w')
# 		for line in file:
# 			#filedata = file.read(line.replace('vdw_lambdas',vdw_lambdas))
# 			if "vdw_lambdas" in line:
# 				filedata.write(vdw_lambdas)
# 			else:
# 				filedata.write(line)
# 				filedata.close()
# 	#    filedata = filedata.replace('vdw_lambdas', vdw_lambdas)
# 	# Write the file out again
# #    with open('file.txt', 'w') as file:
# #        file.write(filedata)


def write_mdp():
	return None

def gen_mdp():
	return None

def import_mdp(mdp_file):
	""" 
	Import MDP file and read lambda values 
	for each lambda type into globar vars.
	"""
	with open(mdp_file,'r+') as f:
		k = f.readlines()      
		for line in k:
			if "coul_lambdas" in line:
				global coul_lambdas_mdp
				coul_lambdas_mdp = line.split(' ')
			#         print(line)
			if "vdw_lambdas" in line:
				global vdw_lambdas_mdp
				coul_lambdas_mdp = line.split(' ')
				print(line)         
				if "mass_lambdas" in line:
					global mass_lambdas_mdp
					mass_lambdas_mdp = line.split(' ')
			#         print(line)
			if "bonded_lambdas" in line:
				global bond_lambdas_mdp
				bond_lambdas_mdp = line.split(' ')
			#         print(line)  
			if "restraint_lambdas" in line:
				global rest_lambdas_mdp
				rest_lambdas_mdp = line.split(' ')
			#         print(line)
			if "temperature_lambdas" in line:
				global temp_lambdas_mdp
				temp_lambdas_mdp = line.split(' ')                         

def coupling(nsim,vdw,coul,mass,bonded,restraint,temp):
	if vdw!=False:
		vdw=vdw.split(',')
		vdwstr='Coupling van der Waals interactions from '+colored(unil+'({0})',"magenta",attrs=["bold"])+' to '+colored(unil+'({1})',"magenta",attrs=["bold"])
		print(vdwstr.format(vdw[0],vdw[1]))
		vdw_lambdas = 'vdw_lambdas = '+" ".join(create_lambdas(nsim,int(vdw[0]),int(vdw[1])))
		println(vdw_lambdas)
	#		import_mdp(nvt)
	if coul!=False:
		coul=coul.split(',')
		coulstr='Coupling Coulomb interactions from '+colored(unil+'({0})',"magenta",attrs=["bold"])+' to '+colored(unil+'({0})',"magenta",attrs=["bold"])
		print(coulstr.format(coul[0],coul[1]))
		coul_lambdas='coul_lambdas = '+" ".join(create_lambdas(nsim,int(coul[0]),int(coul[1])))
		println(coul_lambdas)
	if mass!=False:
		mass=mass.split(',')
		massstr='Coupling mass lambdas from '+colored(unil+'({0})',"magenta",attrs=["bold"])+' to '+colored(unil+'({0})',"magenta",attrs=["bold"])
		print(massstr.format(mass[0],mass[1]))
		mass_lambdas='mass_lambdas = '+" ".join(create_lambdas(nsim,int(mass[0]),int(mass[1])))
		println(mass_lambdas)


def create_lambdas(nsim,start,end):
	""" 
	Creates lambda values partiationed as a linspace. 
	"""
	if start==end:
		print(colored("ERROR: Interval [0,0] specified. \nPlease do not specify intervals for lambdas which won't change during the simulation.","red",attrs=["bold"]))
		exit()
	lambdas = []
	lambda_0 = 1.0/(end-start)
	for i in range(nsim+1):
		lambda_i=lambda_0*(i-start)
		if i>=start and i<=end:
			lambdas.append("{0:.4f}".format(lambda_i))
		else:
			lambdas.append("0.0000")
	return lambdas

def str2array(string):
	"""
	Ignore text and convert numerical
	values in a given string to floats.
	"""
	z=[]
	for k in string:
		try: 
			z.append(float(k))
		except ValueError:
			pass
	return z
#print(str2array(temp_lambdas_mdp))

def println(x):
	"""	
	Prints x in newline sparated by _______ symbols.
	"""
	print(15*'_','\n',x,'\n',15*'_')
	return None

if __name__ == '__main__':
	inputs()