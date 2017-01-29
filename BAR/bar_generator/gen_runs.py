#!/usr/bin/env python3

import click
import os
#from termcolor import colored

root = os.getcwd()

@click.command()

@click.option('-nsim',default=20,help='Total number of simulations.')


@click.option('-root',default=root,help='Simulation root folder. See help for folder structure.')
@click.option('-mdp',default=root+'/MDP',help='MDP folder.')
@click.option('-min1',default=root+'/MDP/EM/em_steep.mdp',help='MDP of steepest decent minimization.')
@click.option('-min2',default=root+'/MDP/EM/em_l-bfgs.mdp',help='MDP of L-BFGS minimization.')
@click.option('-nvt',default=root+'/MDP/NVT/nvt.mdp',help='MDP of NVT equilibration.')
@click.option('-npt',default=root+'/MDP/NPT/npt.mdp',help='MDP of NPT equilibration.')
@click.option('-prod',default=root+'/MDP/Production_MD/md.mdp',help='MDP of the production run.')

@click.option('-ncores',default=8,help='Number of cores per simulation.')
@click.option('-pin',default='auto',help='Pin to cores.')

def inputs(nsim,root,mdp,min1,min2,nvt,npt,prod,ncores,pin):
	gen_mdps(min1,nsim)
	gen_mdps(min2,nsim)
	gen_mdps(nvt,nsim)
	gen_mdps(npt,nsim)
	gen_mdps(prod,nsim)
	gen_jobs(nsim,root,ncores,pin)
	gen_runs(nsim,root)
	gen_bar_results(nsim,root)
	return None
	
def gen_mdps(file,nsim):
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
			wrt1,wrt2 = False,False
			for line in mdp_file:
				if "init_lambda_state" in line and wrt1==False:
					f2.write('init_lambda_state         = '+str(i)+'\n')
					wrt1=True
				elif "vdw_lambdas" in line and wrt2==False:
				'how to add for loop of elifs with lambda type array'
					f2.write('vdw_lambdas ='+" ".join(vdw_lambdas))
					wrt2=False
				else:
					f2.write(line)
	return None

def gen_jobs(nsim,root,ncores,pin):
	"""
	Generates job_{lambda i}.sh files for each lambda in the root dir.
	"""
	jobsh ='''#!/bin/bash
# Set some environment variables 
FREE_ENERGY={root}
echo "Free energy home directory set to $FREE_ENERGY"
MDP=$FREE_ENERGY/../MDP
echo ".mdp files are stored in $MDP"
SYSTEM=$FREE_ENERGY/SYSTEM
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

if __name__ == '__main__':
	inputs()