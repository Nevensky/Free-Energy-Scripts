#!/usr/bin/env python3

import numpy as np
from colorama import Fore,Style
import click
import os
from termcolor import colored
from scipy.interpolate import interp1d
import matplotlib.pyplot as testplt

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

	# print("TEST TEST TEST TEST"*3)
	# test_x,test_y=ddGfunction[1][0],ddGfunction[1][1]
	# #print(test_x,"\n",test_y)
	# test_x = test_x[::-1]
	# test_y = test_y[::-1]
	# test_x_interp = np.linspace(0,1,num=100,endpoint=True)
	# test_x_interp = test_x
	# test_y_interp_func = interp1d(test_x,test_y, bounds_error=False,kind='cubic')
	# test_y_interp = test_y_interp_func(test_x_interp)
	# print(test_x_interp,"\n",test_y_interp)
	# testplt.plot(test_x_interp,test_y_interp)
	# testplt.show()
	# testplt.close()
	# print("TEST TEST TEST TEST"*3)

	for idx,item in enumerate(ddGfunction):
		nsim_i = lambdas_per_interval[0][idx]
#		print(colored("NSIM:","red",attrs=["bold"]),nsim_i)
#		print(item[0][::-1],item[1][::-1])
		if nsim_i != 0:
			zx,zy = item[0],item[1]
			zx,zy = zx[::-1],zy[::-1]
			x,y=bar.create_lambdas_equiG(zx,zy,nsim_i,already_interp=True)
			print("Equidistant lambdas with respect to ddG on interval ",colored(str(intervals[idx]),"yellow",attrs=["bold"]),":\n",colored(ddG_x2,"green",attrs=["bold"]))


if __name__ == '__main__':
	inputs()