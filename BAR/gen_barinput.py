#!/usr/bin/env python3
 
# skripta za generiranje Free Energy BAR inputa
 
import sys, getopt
import click
from colorama import Fore,Style
import os.path
import numpy as np
from scipy.interpolate import interp1d
#from numpy import interp as interp1d
import matplotlib.pyplot as plt
import re 

import ddGintervals
import importbar as bar
 
unil='\u03BB'
 
@click.command()

@click.option('-nsim',default=20,help='Total number of simulations.')
@click.option('-vdw',default=False,help='Mutate van der Waals interactions.')
@click.option('-coul',default=False,help='Mutate Coulomb interactions.')
@click.option('-mass',default=False,help='Mutate mass of species.')
@click.option('-bonded',default=False,help='Mutate bonded interactions.')
@click.option('-restraint',default=False,help='Mutate restrained interactions.')
@click.option('-temp',default=False,help='Simulated temptering.')

@click.option('-mdp',default='./MDP',help='MDP folder.')
@click.option('-min1',default='./MDP/EM/em_steep.mdp',help='MDP of steepest decent minimization.')
@click.option('-min2',default='./MDP/EM/em_l-bfgs.mdp',help='MDP of L-BFGS minimization.')
@click.option('-nvt',default='./MDP/NVT/nvt.mdp',help='MDP of NVT equilibration.')
@click.option('-npt',default='./MDP/NPT/npt.mdp',help='MDP of NPT equilibration.')
@click.option('-prod',default='./MDP/Production_MD/md.mdp',help='MDP of the production run.')
@click.option('-ncores',default=4,help='Number of cores per simulation.')
@click.option('-pin',default='auto',help='Pin to cores.')
 
def inputs(nsim,vdw,coul,mass,bonded,restraint,temp,mdp,min1,min2,nvt,npt,prod,ncores,pin):
	'''Free Energy BAR input generator. \n

Developed by Neven Golenic | neven.golenic@gmail.com

\b
Number of simulations [ N ] default: [ 20 ]
\u03BB parameters [ N_min,N_max ]

If BAR results from a previous simulation are imported, lambdas of different type must not overlap!

\b
Recommended inital folder structure:
.
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
'''
	print('')
#    bonded=bonded.split(',')
#    restraint=restraint.split(',')
#    temp=temp.split(',')
	check_mdp(mdp,min1,min2,nvt,npt,prod) # ovo je dobro treba reaktivirati

	if vdw!=False:
		vdw=vdw.split(',')
		vdwstr='Coupling van der Waals interactions from '+Fore.MAGENTA+unil+'({0})'+Fore.RESET+' to '+Fore.MAGENTA+unil+'({1})'+Fore.RESET
		print(vdwstr.format(vdw[0],vdw[1]))
		vdw_lambdas = 'vdw_lambdas = '+" ".join(create_lambdas(nsim,int(vdw[0]),int(vdw[1])))+'\n'
		print(10*'_','\n',vdw_lambdas,10*'_','\n')
		import_mdp(min1,min2,nvt,npt,prod,vdw_lambdas)
	if coul!=False:
		coul=coul.split(',')
		coulstr='Coupling Coulomb interactions from '+Fore.MAGENTA+unil+'({0})'+Fore.RESET+' to '+Fore.MAGENTA+unil+'({1})'+Fore.RESET
		print(coulstr.format(coul[0],coul[1]))
		coul_lambdas='coul_lambdas = '+" ".join(create_lambdas(nsim,int(coul[0]),int(coul[1])))
		print(10*'_','\n',coul_lambdas,'\n',10*'_','\n')
	if mass!=False:
		mass=mass.split(',')
		massstr='Coupling mass lambdas from '+Fore.MAGENTA+unil+'({0})'+Fore.RESET+' to '+Fore.MAGENTA+unil+'({1})'+Fore.RESET
		print(massstr.format(mass[0],mass[1]))
		mass_lambdas='mass_lambdas = '+" ".join(create_lambdas(nsim,int(mass[0]),int(mass[1])))
		print(1*'_','\n',mass_lambdas,'\n',10*'_','\n')
	
#	return None
 
# click.echo(vdw)

# def import_mdp(min1,min2,nvt,npt,prod):
#     if os.path.isdir("./MDP"):
#         with open(min1) as file:
#             min1_lines=[]
#                 for line in file:
#                     min1_lines.append(line.strip().split(','))
#                     if line =="":
#     return None

def check_mdp(mdp,min1,min2,nvt,npt,prod):
	# Check if MDP file exists or if all .mdp files were specified manually.
	if os.path.exists(min1) and os.path.exists(min2) and os.path.exists(nvt) and os.path.exists(npt) and os.path.exists(prod):
		pass	
	elif os.path.isdir(mdp) :
		if os.path.exists(mdp+'/EM/em_steep.mdp') and os.path.exists(mdp+'/EM/em_l-bfgs.mdp') and os.path.exists(mdp+'/NVT/nvt.mdp')==True and os.path.exists(mdp+'/NPT/npt.mdp') and os.path.exists(mdp+'/Production_MD/md.mdp'):
			pass
		else: 
			print(Fore.RED,'ERROR: MDP folder corrupt. Expected folder structure:',Fore.YELLOW,"""\n.
└── MDP
    |
    ├─ EM
    |  ├── em_steep.mdp
    |  └── em_l-bfgs.mdp
    ├─ NVT
    |   └── nvt.mdp 
    ├─ NPT
    |   └── npt.mdp 
    └─ Production_MD
        └── md.mdp""",Fore.RESET)
			exit()
		#print("MDP exists!")
	else:
		print(Fore.RED,"ERROR: MDP folder not found, please specify paths to ALL .mdp files manually.",Fore.RESET)
		exit()

def import_mdp_old(min1,min2,nvt,npt,prod,vdw_lambdas):
#   pass # jer je ova funkcija losa
	# Imports mdp files and replaces them with new ones with correct lambdas
	with open(min1) as file:
		filedata = open('MDP/em_steep.mdp','w')
		for line in file:
			#filedata = file.read(line.replace('vdw_lambdas',vdw_lambdas))
			if "vdw_lambdas" in line:
				filedata.write(vdw_lambdas)
			else:
				filedata.write(line)
		filedata.close()
	#    filedata = filedata.replace('vdw_lambdas', vdw_lambdas)
	# Write the file out again
#    with open('file.txt', 'w') as file:
#        file.write(filedata)


def write_mdp():
	return None

def gen_mdp():
	return None
#ncpus = 1
#nsim = 20

def import_mdp(mdp_file):
   global coul_lambdas_mdp
   global vdw_lambdas_mdp
   global mass_lambdas_mdp
   global rest_lambdas_mdp
   global bond_lambdas_mdp
   global temp_lambdas_mdp
   with open(mdp_file,'r+') as f:
      k = f.readlines()      
   for line in k:
      if "coul_lambdas" in line:
         coul_lambdas_mdp = line.split(' ')
#         print(line)
      if "vdw_lambdas" in line:
         coul_lambdas_mdp = line.split(' ')
#         print(line)         
      if "mass_lambdas" in line:
         mass_lambdas_mdp = line.split(' ')
#         print(line)
      if "bonded_lambdas" in line:
         bond_lambdas_mdp = line.split(' ')
#         print(line)  
      if "restraint_lambdas" in line:
         rest_lambdas_mdp = line.split(' ')
#         print(line)
      if "temperature_lambdas" in line:
         temp_lambdas_mdp = line.split(' ')                         

def str2array(string):
   '''
   Ignore text and convert numerical
   values in a given string to floats.
   '''
   z=[]
   for k in string:
      try: 
         z.append(float(k))
      except ValueError:
         pass
   return z

#print(str2array(temp_lambdas_mdp))

         
 
def create_lambdas(nsim,start,end):
	""" Creates lambda values equidistant with respect to lambdas """
	lambdas = []
	lambda_0 = 1.0/(end-start)
	for i in range(nsim+1):
		lambda_i=lambda_0*(i-start)
		if i>=start and i<=end:
			lambdas.append("{0:.4f}".format(lambda_i))
		else:
			lambdas.append("0.0000")
	return lambdas
 
 
if __name__ == '__main__':
	inputs()