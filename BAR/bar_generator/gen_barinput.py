#!/usr/bin/env python3

# skripta za generiranje Free Energy BAR inputa

import sys, getopt
import click
import os
import numpy as np
import re 
from termcolor import colored

# import ddGintervals
# import importbar as bar
import bar_main as barmain

root = os.getcwd()
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

@click.command(context_settings=dict(help_option_names=['-h', '-help']))

@click.option('-nsim',default=20,help='Total number of simulations.')
@click.option('-vdw',default=False,help='Mutate van der Waals interactions.')
@click.option('-coul',default=False,help='Mutate Coulomb interactions.')
@click.option('-mass',default=False,help='Mutate mass of species.')
@click.option('-bonded',default=False,help='Mutate bonded interactions.')
@click.option('-restraint',default=False,help='Mutate restrained interactions.')
@click.option('-temp',default=False,help='Simulated temptering.')

@click.option('-root',default=root,help='Simulation root folder. See help for folder structure.')
@click.option('-mdp',default=root+'/MDP',help='MDP folder.')
@click.option('-min1',default=root+'/MDP/EM/em_steep.mdp',help='MDP file of steepest decent minimization.')
@click.option('-min2',default=root+'/MDP/EM/em_l-bfgs.mdp',help='MDP file of L-BFGS minimization.')
@click.option('-nvt',default=root+'/MDP/NVT/nvt.mdp',help='MDP file of NVT equilibration.')
@click.option('-npt',default=root+'/MDP/NPT/npt.mdp',help='MDP file of NPT equilibration.')
@click.option('-prod',default=root+'/MDP/Production_MD/md.mdp',help='MDP file of the Production run.')

@click.option('-bar_file',default=root+'/analysis/bar_results.txt',help='BAR results file from a previous simulation.')

@click.option('-ncores',default=8,help='Number of cores per simulation.')
@click.option('-pin',default='auto',help='Pin to cores.')


def inputs(nsim,vdw,coul,mass,bonded,restraint,temp,root,mdp,min1,min2,nvt,npt,prod,bar_file,ncores,pin):
	"""Free Energy BAR input generator. \n
Developed by Neven Golenic | neven.golenic@gmail.com

\b
Number of simulations [ N ] default:  20 
λ coupling interval [N_min,N_max] default: 0,0
If BAR results from a previous simulation are imported, lambda vectors of different types must not overlap!

\b
Recommended inital folder structure:
└── analysis
|	|
|	└─ bar
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
	└─ system.gro
	"""

	check_dirs(root)
	mdp_types= [min1,min2,nvt,npt,prod]
	check_mdps(mdp,mdp_types)
	cdict = coupling(nsim,vdw,coul,mass,bonded,restraint,temp,bar_file)
	cdict = fill_inactive_lambdas(nsim,cdict)
#	print(cdict)

	gen_mdps(min2,nsim,cdict)
	gen_mdps(nvt,nsim,cdict)
	gen_mdps(npt,nsim,cdict)
	gen_mdps(prod,nsim,cdict)
	gen_jobs(nsim,root,ncores,pin)

	run_at_once = False # Currently not implemented
	gen_runs(nsim,root,run_at_once=run_at_once)
	if run_at_once:
		gen_bar_results(nsim,root)

	return None

############################
# Check folder integrity   #
############################

def check_dirs(root):
	"""
	Check if given directories exist and create them if they don't.
	"""
	dir_array = [root+"/analysis",root+"/analysis/bar",root+"/JOBS",root+"SYSTEM"]
	for item in dir_array:
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


####################################################
# Coupling ddictionary and lambda vector creation. #
####################################################

def coupling(nsim,vdw,coul,mass,bonded,restraint,temp,bar_file):
	"""
	Defines a coupling dictionary which contains specified lambda intervals.
	Intervals are replaced by lambda vectors and written back to the dictionary. 
	"""
	coupling_dict = {'vdw_lambdas':vdw,'coul_lambdas':coul,'mass_lambdas':mass,'bonded_lambdas':bonded,'restraint_lambdas':restraint,'temperature_lambdas':temp}	
	for key in coupling_dict:
		if coupling_dict[key] != False:
			lint = coupling_dict[key].split(',') #split coupling interval
			lstring = 'Coupling {type} from '+color('λ({0})',"magenta")+' to '+color('λ({1})',"magenta")
			print(lstring.format(lint[0],lint[1],type=key))
			coupling_dict[key] = key+" = "+" ".join(create_lambdas(nsim,int(lint[0]),int(lint[1]),bar_file=bar_file))
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

def create_lambdas(nsim,start,end,bar_file):
	""" 
	Creates lambda values partiationed as a linspace if bar_results file not provided.
	Otherwise create equidistant lambdas with respect to ddG(lambda). 
	"""
	if start==end:
		print(color("ERROR: Interval [0,0] specified. \nPlease do not specify intervals for inactive lambda types.","red"))
		exit()	

	lambdas = []	
	if os.path.isfile(bar_file)==False:
		lambda_0 = 1.0/(end-start)
		for i in range(nsim+1):
			lambda_i=lambda_0*(i-start)
			if i>=start and i<=end:
				lambdas.append("{0:.4f}".format(lambda_i))
			elif i<end:
				lambdas.append("0.0000")
			elif i>end:
				lambdas.append("1.0000")
	else:
		nsim_per_interval = int(end-start)
		lambda_matrix = barmain.inputs(nsim_per_interval,bar_file)
		lambda_string = [str(k) for k in lambda_matrix]
		for i in range(start):
			lambdas.append("0.0000")		
		for k in lambda_matrix:
			lambdas.append("{0:.4f}".format(k))
		for i in range(nsim-end):
			lambdas.append("1.0000")
	return lambdas

########################################
# MDP and executable file generation  #
########################################

def gen_nsim_mdps(nsim,mdp_types):
	""" 
	!!! ne valja zamijeniti s gen_mdps()
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
							f2.write('init_lambda_state      = '+str(i)+'\n')
						else:
							f2.write(line)

def gen_mdps(file,nsim,cdict):
	"""
	Generates MDP files of each type in /MDP folder,
	sets init_lambda_state to current lambda up to nsim
	and updates the lambda vector for each lambda type.
	!!! Should take array of lambda types and lambda type names
	"""
#	print('Executing gen_mdps()')
	file_dir = os.path.dirname(file)
	file_name = os.path.basename(file).split('.')[0]
	mdp_file = []
	with open(file,'r') as f:
		mdp_file = f.readlines()
	for i in range(nsim):
		new_file_path = str(file_dir+'/'+file_name+'_'+str(i)+'.mdp')
		with open(new_file_path,'w') as f2:
#			print("Written to file: "+new_file_path)
			wrt1 = False
			wrt2 = [False for k in range(6)]
			z1 = 0
			for line in mdp_file:
				if "init_lambda_state" in line and wrt1==False:
					f2.write('init_lambda_state         = '+str(i)+'\n')
					wrt1=True
				else:
					wrt3 = False
					for key in cdict:
						if key in line and wrt2[z1]==False:
							f2.write(cdict[key]+"\n")
							wrt2[z1]=True
							z1 += 1
						# explicitly stated conditions (very bad!), needs rewritting
						elif wrt3==False and "vdw_lambdas" not in line and "coul_lambdas" not in line and "mass_lambdas" not in line and "bonded_lambdas" not in line and "restraint_lambdas" not in line and "temperature_lambdas" not in line:
							f2.write(line)
							wrt3 = True
	return None

def gen_jobs(nsim,root,ncores,pin):
	"""
	Generates job_{lambda i}.sh files for each lambda in the root dir.
	"""
	jobsh ='''#!/bin/bash
# Set some environment variables 
FREE_ENERGY={root}
echo "Free energy home directory set to $FREE_ENERGY"
MDP=$FREE_ENERGY/MDP
echo ".mdp files are stored in $MDP"
SYSTEM=$FREE_ENERGY/SYSTEM
# unset OpenMP threads enviroment var to avoid clashes with Gromacs parallelization
unset OMP_NUM_THREADS
# A new directory will be created for each value of lambda and
# at each step in the workflow for maximum organization.
mkdir lambda_{lambda_i}
cd lambda_{lambda_i}
#################################
# ENERGY MINIMIZATION 1: STEEP  #
#################################
echo "Starting minimization for lambda = {lambda_i}..." 
mkdir EM_1 
cd EM_1
# Iterative calls to grompp and gmx mdrun to run the simulations
gmx grompp -f $MDP/EM/em_steep_{lambda_i}.mdp -c $SYSTEM/system.gro -p $SYSTEM/system.top -o min_{lambda_i}.tpr
gmx mdrun -nt {ncpus} -deffnm min_{lambda_i}
sleep 5
#################################
# ENERGY MINIMIZATION 2: L-BFGS #
#################################
cd ../
mkdir EM_2
cd EM_2
# We use -maxwarn 1 here because grompp incorrectly complains about use of a plain cutoff; this is a minor issue
# that will be fixed in a future version of Gromacs
gmx grompp -f $MDP/EM/em_l-bfgs_{lambda_i}.mdp -c ../EM_1/min_{lambda_i}.gro -p $SYSTEM/system.top -o min_{lambda_i}.tpr -maxwarn 1
# Run L-BFGS in serial (cannot be run in parallel)
gmx mdrun -nt 1 -deffnm min_{lambda_i}
echo "Minimization complete."
sleep 5
#####################
# NVT EQUILIBRATION #
#####################
echo "Starting constant volume equilibration..."
cd ../
mkdir NVT
cd NVT
gmx grompp -f $MDP/NVT/nvt_{lambda_i}.mdp -c ../EM_2/min_{lambda_i}.gro -p $SYSTEM/system.top -o nvt_{lambda_i}.tpr 
gmx mdrun -nt {ncpus} -pin {pin} -deffnm nvt_{lambda_i}
echo "Constant volume equilibration complete."
sleep 5
#####################
# NPT EQUILIBRATION #
#####################
echo "Starting constant pressure equilibration..."
cd ../
mkdir NPT
cd NPT
gmx grompp -f $MDP/NPT/npt_{lambda_i}.mdp -c ../NVT/nvt_{lambda_i}.gro -p $SYSTEM/system.top -t ../NVT/nvt_{lambda_i}.cpt -o npt{lambda_i}.tpr 
gmx mdrun -nt {ncpus} -pin {pin} -deffnm npt_{lambda_i}
echo "Constant pressure equilibration complete."
sleep 5
#################
# PRODUCTION MD #
#################
echo "Starting production MD simulation..."
cd ../
mkdir Production_MD
cd Production_MD
gmx grompp -f $MDP/Production_MD/md_{lambda_i}.mdp -c ../NPT/npt_{lambda_i}.gro -p $SYSTEM/system.top -t ../NPT/npt_{lambda_i}.cpt -o md_{lambda_i}.tpr 
gmx mdrun -nt {ncpus} -pin {pin} -deffnm md_{lambda_i}
echo "Production MD complete."
# End
echo "Ending. Job completed for lambda = {lambda_i}"
	'''
	for i in range(nsim):
		job=jobsh.format(root=root,ncpus=ncores,pin=pin,lambda_i=str(i))
		with open(root+'/job_'+str(i)+'.sh','w') as f:
			f.write(job)
	return None

def gen_runs(nsim,root,run_at_once=False):
	""" Generate exec_runs.sh script to run simulations for all lambdas. """
	runall='rm -rf lambda_* ;'
	for i in range(nsim):
		runall += '\n {root}/job_{lambda_i}.sh ;'.format(root=root,lambda_i=str(i))
	with open(root+'/exec_runs.sh','w') as f:
		f.write(runall)
	if run_at_once:
		os.system('nohup '+root+'/exec_runs.sh &')
	return None

def gen_bar_results(nsim,root):
	"""
	Generates dhdl xvgs from each lambda production run
	and copies them to analysis/bar then runs gmx bar on them.
	"""
	bar_dir = root+'/analysis/bar'
	mdxvgs = ''
	for i in range(nsim):
		prod_dir = root+'lambda_'+str(i)+'/Production_MD'
		os.system('cp '+prod_dir+'/md_{}.xvg '.format(str(i))+bar_dir)
		mdxvgs += bar_dir+'/md_{}.xvg '.format(str(i))
	os.system('gmx bar -o -oi -oh -f '+mdxvgs+' > '+root+'/analysis/bar_results.txt')
	return None

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


#######################
# Formating functions #
#######################

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