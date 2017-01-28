#!/usr/bin/env python3

# skripta za generiranje Free Energy BAR inputa

import sys, getopt
import click
import os
import numpy as np
#from scipy.interpolate import interp1d
#from numpy import interp as interp1d
import matplotlib.pyplot as plt
import re 
#from termcolor import colored

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


tree = """.
└── analysis
|		|
|		└─ bar
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
└── SYSTEM
		|
		├─ system.top
		└─ system.gro"""

def inputs(nsim,vdw,coul,mass,bonded,restraint,temp,root,mdp,min1,min2,nvt,npt,prod,ncores,pin):
	"""Free Energy BAR input generator. \n
Developed by Neven Golenic | neven.golenic@gmail.com
\b
Number of simulations [ N ] default:  20 
\u03BB coupling interval [N_min,N_max] default: 0,0
If BAR results from a previous simulation are imported, lambda vectors of different types must not overlap!
\b
Recommended inital folder structure:
	"""+tree

	check_dirs(root)
	mdp_types= [min1,min2,nvt,npt,prod]
	check_mdps(mdp,mdp_types)
	cdict = coupling(nsim,vdw,coul,mass,bonded,restraint,temp)
	cdict = fill_inactive_lambdas(nsim,cdict)
	print(cdict)

	return None

def check_dirs(root):
	"""
	Check if given directories exist and create them if they don't.
	"""
	dir_array = [root+"/analysis",root+"/analysis/bar",root+"/jobs"]
	for item in dir_array
		try:
			os.makedirs(item)
		except FileExistsError:
			pass
	return None

def check_mdps(mdp,mdp_types):
	""" 
	Check if MDP file exists or if all 
	.mdp files were specified manually 
	"""
	if os.path.isdir(mdp):
		pass
	else:
		print(color("ERROR: MDP folder corrupt. Expected folder structure:\n"+tree,"red"))
		
	for item in mdp_types:
		if os.path.exists(item):
			pass
		else:
			print(color("ERROR: MDP folder corrupt. Expected folder structure:\n"+tree,"red"))
			exit()

#		auto_initialize = input("Intialize MDP files automatically? yes/no :")
#		if auto_initialize=="yes" or auto_initialize=="y":
#			print("Currently not implemented. You must intialize MDP files manually.")
#			exit()

def gen_nsim_mdps(nsim,mdp_types):
	""" 
	! ne valja zamijeniti
	Generates mdp files of each type for 
	each simulation up to total nsim.
	"""
	mdp_types= [min1,min2,nvt,npt,prod]
	for mdp_type in mdp_types:
		with open(mdp_type,'r') as f1:
			mdp_i = f1.read()
			for i in range(nsim+1):
				mdp_path_i = mdp_type.replace('.mdp','')+str(i)+'.mdp'
				with open(mdp_path_i,'a') as f2: 
					for line in mdp_i:
						if 'init_lambda_state' in line:
							f2.write('init_lambda_state      = '=str(i)='\n')
						else:
							f2.write(line)

def import_mdp(mdp_file,dict):
	""" 
	Import a MDP file and read lambda values 
	for each lambda type if coupling not specified.
	"""
	with open(mdp_file,'r+') as f:
		k = f.readlines()      
		for line in k:
			for key in dict:
				if (dict[key] == False) and (key in line):
					dict[key] = key+" = "+" ".join(str2array(line))
	return dict

def coupling(nsim,vdw,coul,mass,bonded,restraint,temp):
	"""
	Defines a coupling dictionary which contains specified lambda intervals.
	Intervals are replaced by lambda vectors and written back to the dictionary. 
	"""
	coupling_dict = {'vdw_lambdas':vdw,'coul_lambdas':coul,'mass_lambdas':mass,'bonded_lambdas':bonded,'restraint_lambdas':restraint,'temperature_lambdas':temp}	
	for key in coupling_dict:
		if coupling_dict[key] != False
			lint = coupling_dict[key].split(',') #split coupling interval
			lstring = 'Coupling {type} from '+color(unil+'({0})',"magenta")+' to '+color(unil+'({1})',"magenta")
			print(lstring.format(lint[0],lint[1],type=key))
			coupling_dict[key] = key+" = "+" ".join(create_lambdas(nsim,int(lint[0]),int(lint[1])))
			println(coupling_dict[key])
	return coupling_dict
	
def fill_inactive_lambdas(nsim,dict):
	""" 
	Fill lambdas values with zeros up to nsim for all inactive lambda types. 
	"""
	for key in dict:
		if dict[key]==False:
			dict[key] = key+' = '+nsim*' 0.0000'
	return dict

def create_lambdas(nsim,start,end,bar):
	""" 
	Creates lambda values partiationed as a linspace if bar_results file not provided.
	Otherwise create equidistant lambdas with respect to ddG(lambda). 
	"""
	if start==end:
		print(color("ERROR: Interval [0,0] specified. \nPlease do not specify intervals for inactive lambda types.","red"))
		exit()	

	lambdas = []	
	if bar==False:
		lambda_0 = 1.0/(end-start)
		for i in range(nsim+1):
			lambda_i=lambda_0*(i-start)
			if i>=start and i<=end:
				lambdas.append("{0:.4f}".format(lambda_i))
			else:
				lambdas.append("0.0000")
	else:
		# MISSING CALL TO barmain()
		exit()
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

def color(x,c="yellow"):
	"""
	Colors and bolds string for print() output.
	"""
	x = colored(x,c,attrs=["bold"])
	return x

def println(x):
	"""	
	Prints x in newline sparated by _______ symbols.
	"""
	print(15*'_','\n',x,'\n',15*'_')
	return None

if __name__ == '__main__':
	inputs()