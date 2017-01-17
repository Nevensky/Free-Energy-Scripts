#!/usr/bin/env python3

import numpy as np
from colorama import Fore,Style
import click
import os
from termcolor import colored
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import ddGintervals as barint
import importbar as bar

#print(colored(outut,"yellow",attrs=["bold"]))

@click.command()

@click.option('-nsim',default=20,help='Total number of simulations.')
@click.option('-bar_file',default=os.getcwd()+'/bar_results.txt',help='BAR results file from previous simulation.')
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
#	barint.plot_intervals(intervals,ddGfunction)

	for idx,item in enumerate(ddGfunction):
		nsim_i = lambdas_per_interval[0][idx]
		if nsim_i != 0:
			x,y=bar.create_lambdas_equiG(item[0][::-1],item[1][::-1],nsim_i,already_interp=True)
			print("Equidistant lambdas with respect to ddG on interval ",colored(str(intervals[idx]),"yellow",attrs=["bold"]),":\n",colored(x,"green",attrs=["bold"]))
			plt.plot(item[0][::-1],item[1][::-1],'-',markersize=5,fillstyle='none')
			plt.plot(x,y,'o',markersize=5)
			#plt.plot(ddG_x_interp, ddG_y_interp, '-')
	plt.show()
	plt.close()


if __name__ == '__main__':
	inputs()