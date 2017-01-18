#!/usr/bin/env python3
# skripta za određivanje dG intervala funkcija koje ne rastu/padaju monotono
# trenutno ne prisljava da se prvom intervalu dodijeli ijedna lambda -> vidjeti reperkusije

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy import interpolate

import importbar as bar

plt.style.use('ggplot')


def find_intervals(x,y):
	""" Find and divide intervals where ddG is a monotonously ascending or descening function. 
	The input should correspond to x = ddG_x_interp and y =ddG_y_interp."""
	
	# DEBUG
	# x = np.linspace(0,20,num=500)
	intervals = []
	dividedfunc = []
	start_i=0
	dy = np.gradient(y)
	k=0
	for i in range(len(dy)-1):
		if np.sign(dy[i]) != np.sign(dy[i+1]):
			k+=1
			print("monotonic ddG(x) on interval: [",x[start_i],x[i],"]")
			intervals.append([x[start_i],x[i]])
			dividedfunc.append([x[start_i:i],y[start_i:i]])
			#plt.plot(x[start_i:i],y[start_i:i])
			start_i = i	
	# If ddG is monotnous over the whole range return one interval
	if k==0:
		intervals.append([0.00,1.00])
		dividedfunc.append([x,y])
#		plt.plot(x,y)
	# fix last interval
	if start_i != len(dy)-1:
		print("monotonic ddG(x) on interval: [",x[start_i],x[-1],"]")
		intervals.append([x[start_i],x[-1]])
		dividedfunc.append([x[start_i:],y[start_i:]])
		#plt.plot(x[start_i:],y[start_i:])
		return intervals,dividedfunc


def interval_weights(dividedfunc):
	""" Calculate how many lambdas should span each interval
	based on the maximal change in ddG. """
	weights = []
	for item in dividedfunc:
		ddGfunction_i = item[1]
		weights.append(np.abs(np.amax(ddGfunction_i)-np.amin(ddGfunction_i)))
	weights= weights/np.sum(weights)
	return weights

def lambdas_per_interval(nsim,weights):
	"""Assign proper number of lambdas per interval	based on calculated wieghts. """
	num_lambdas = []
	rounding_errors = []
	for item in weights:
		num_lambdas.append(int(np.floor(item*nsim)))
		rounding_errors.append(np.abs(item*nsim-np.floor(item*nsim)))
	sum_lambdas = np.sum(num_lambdas)
	if sum_lambdas<nsim:
		while True:
			print("Sum of weighted lambdas: ",sum_lambdas,"\n\tRounding errors: ",rounding_errors,"\n\tMax:",np.amax(rounding_errors),"Max_idx: ",np.argmax(rounding_errors))
			if sum_lambdas<nsim:
				if np.abs(sum_lambdas-nsim)>2 or num_lambdas[-1]!=0 or num_lambdas[0]!=0:
					max_idx = np.argmax(rounding_errors)
					num_lambdas[max_idx] += 1
					rounding_errors[max_idx] = 0
					sum_lambdas = np.sum(num_lambdas)
				elif num_lambdas[0]==0 and np.abs(sum_lambdas-nsim)==2:
					print("Fixed first interval instead of weighted interval with corresponding maximum variance.")
					num_lambdas[0] += 1
					sum_lambdas = np.sum(num_lambdas)
				elif num_lambdas[-1]==0 and np.abs(sum_lambdas-nsim)==1:
					print("Fixed last interval instead of weighted interval with corresponding maximum variance.")
					num_lambdas[-1] += 1
					sum_lambdas = np.sum(num_lambdas)
			else:
				break
	elif sum_lambdas>nsim:
		print("EROOR: Could not assign weighted number of lambdas per each interval. \n Nsim > N(weighted lambdas)")
		raise SystemExit(0) # prekida program
	return num_lambdas,sum_lambdas

def plot_intervals(intervals,dividedfunc,close=True,save=True):
	""" Plot ddG interpolated function colored according to intervals 
	and separated by vertical lines. """
	plt.rc('text', usetex=True)
	axes = plt.gca()
	# fig, ax = plt.subplots(figsize=(10, 8))
	for pl in dividedfunc:
		plt.plot(pl[0],pl[1])
	
	ymin,ymax = axes.get_ylim()
	for pl in dividedfunc:
		plt.vlines(pl[0][0],ymax=pl[1][0],ymin=ymin,linestyle='dotted',color='gray',linewidth=1)

	plot_labels = ["[{:.2f}, {:.2f}]".format(k[0],k[1]) for k in intervals]
	plt.xlabel(r'$\lambda$') 
	plt.ylabel(r"$\Delta \Delta G / \mathrm{kJ mol^{-1}}$")
	plt.legend(plot_labels, loc='best')
#	plt.show()
	if save:
		plt.savefig("ddG_intervals.pdf")
	if close:
		plt.close()

if __name__ == '__main__':
	# DEBUG
	# x = np.linspace(0,1,100)
	# y = np.sin(12*x)-2*x
	nsim = 30
	ddG_x,ddG_y = bar.import_bar('./bar_results.txt',nsim)
	x,y = bar.interpolate_ddG(ddG_x,ddG_y,nsim)

	intervals,ddGfunction = find_intervals(x,y)
	weights = interval_weights(ddGfunction)
	print("Weights: ",weights)
	print("Number of lambdas per interval & sum: ",lambdas_per_interval(nsim,weights))
	plot_intervals(intervals,ddGfunction)

	ddG_x2,ddG_y2 = bar.create_lambdas_equiG(ddG_x,ddG_y,nsim)
	print("Equidistant lambdas with respect to ddG:\n",ddG_x2)


