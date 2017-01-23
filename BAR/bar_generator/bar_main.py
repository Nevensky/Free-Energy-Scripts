#!/usr/bin/env python3

import numpy as np
import click
import os
from termcolor import colored
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib

import ddGintervals as barint
import importbar as bar

matplotlib.rcParams['text.usetex'] = True

@click.command()
@click.option('-nsim',default=20,help='Total number of simulations.')
@click.option('-bar_file',default=os.getcwd()+'/bar_results.txt',help='BAR results file from a previous simulation.')
#@click.option('-dir',default=os.getcwd(),help='Simulation root directory.')

def inputs(nsim,bar_file):
	# interpolate inital function over the whole interval [0,1]	
	ddG_x,ddG_y = bar.import_bar(bar_file,nsim)
#	print(ddG_x,ddG_y)
	ddG_x_interp,ddG_y_interp = bar.interpolate_ddG(ddG_x,ddG_y,nsim)

	ddG_x2,ddG_y2 = bar.create_lambdas_equiG(ddG_x,ddG_y,nsim)
	#bar.plot_interpolation(ddG_y,ddG_x,ddG_x2,ddG_y2,ddG_x_interp,ddG_y_interp)
	print("Equidistant lambdas with respect to ddG on interval [0,1]:\n",colored(ddG_x2,"yellow",attrs=["bold"]))
	
	# partition function into intervals
	intervals,ddGfunction = barint.find_intervals(ddG_x_interp,ddG_y_interp)
	weights = barint.interval_weights(ddGfunction)
	lambdas_per_interval = barint.lambdas_per_interval(nsim,weights)
	print("Weights: ",weights)
	print("Number of lambdas per interval & sum: ",lambdas_per_interval)
	barint.plot_intervals(intervals,ddGfunction,close=False,save=False)
	
	plt.plot(ddG_x[0],ddG_y[0],'o',markersize=5)
	lambda_matrix = partition_lambdas(ddGfunction,intervals,lambdas_per_interval)
	print(70*"~",u"\n ########  Final lambdas equidistant with respect to \u0394\u0394G(\u03BB)  ########","\n",70*"~","\n",colored(lambda_matrix,"green",attrs=["bold"]),"\n",70*"~")

def partition_lambdas(ddGfunction,intervals,lambdas_per_interval):
	""" Outputs equidistant lambdas with respect to ddG. 
	!! Currently also implements plot_ddG(x,y) functionality."""
	plt.rc('legend',**{'fontsize':6})
	lambdas = []
	plot_labels = [r"$\Delta\Delta G{}([{:.2f}, {:.2f}])$".format("_\\mathrm{interp}",k[0],k[1]) for k in intervals]
	plot_labels.append(r"$\Delta\Delta G(\lambda=0)$")
	for idx,item in enumerate(ddGfunction):
		lints_i = lambdas_per_interval[0][idx]
		if lints_i != 0:
			x,y=bar.create_lambdas_equiG(item[0][::-1],item[1][::-1],lints_i,already_interp=True)
			lambdas.append(x)
			print("Equidistant lambdas with respect to ddG on interval: ",colored("[{:.2f}, {:.2f}]".format(intervals[idx][0],intervals[idx][1]),"yellow",attrs=["bold"]),"\n",colored(x,"green"))
			# plt.plot(item[0][::-1],item[1][::-1],'-')
			plt.plot(x,y,'o',markersize=5)
			plot_labels.append(r"$\Delta\Delta G(\lambda \vert \lambda \subset [{:.2f}, {:.2f}])$".format(intervals[idx][0],intervals[idx][1]))
		else:
			print("Skipping interval:",colored("[{:.2f}, {:.2f}]".format(intervals[idx][0],intervals[idx][1]),"red",attrs=["bold"]))

	# fix first and last value of lambda, λ(0)=0 λ(1)=1	
	lambdas = np.concatenate(lambdas)
	lambdas[0] = 0
	lambdas[-1] = 1
	
	#plt.show()
	plt.legend(plot_labels,loc='best')
	plt.xlabel(r'$\lambda$')
	plt.ylabel(r"$\Delta \Delta G / \mathrm{kJ mol^{-1}}$")
	plt.savefig("analysis/ddG_plot.pdf")
	plt.close()
	return lambdas

def plot_ddG(x,y):
	""" Plots inputs and equidistant lambdas 
	(with respect to ddG) on partitioned intervals. """
	pass

if __name__ == '__main__':
	inputs()