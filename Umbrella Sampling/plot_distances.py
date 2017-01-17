#!/usr/bin/env python3
#skripta za plotanje udaljenosti i generiranje optimalnog razmaka konfiguracija
# generira pocetni email report u analysis/email_report.txt

import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
import numpy as np
import click
import os #to get workign dir
import socket #get computer name

plt.style.use('ggplot')

#max_dist=1.1 # nm
#min_sep=0.2 # nm

@click.command()

@click.option('-max_dist',default=1.0,help='Maximum distance [nm] of coord solvent from starting position.')
@click.option('-min_dist',default=False,help='Minimum distance [nm] of coord solvent from starting position.')
@click.option('-min_sep',default=0.2,help='Minimum sepearation distance [nm] between subsequent windows.')
@click.option('-omit',default=0,help='Omit all confs up to defined number.')
@click.option('-umbrella_dir',default=os.getcwd(),help='Path to simulation root dir.')
@click.option('-nskip',default=1,help='How many frames were skipped in trjconv.')

def inputs(max_dist,min_dist,min_sep,umbrella_dir,nskip,omit):
	'''Umbrella sampling mid-step. Plot distances between coord solv and ligand for all generated configurations and select meaningful ones based on 
	maximum distance and window sepearation distance. \n

Developed by Neven Golenic | neven.golenic@gmail.com'''
	
	path = umbrella_dir+"/CONF"
#	filenames = [item for item in next(os.walk(path))[2] if ".gro" in item]
#	num_files = len(filenames)-1 # ide od 0 - 5000 .... jedna viska za plotanje

	# max_dist=float(max_dist)
	# min_dist=float(min_dist)
	# min_sep=float(min_sep)
	# omit=int(omit)
	data = np.genfromtxt(umbrella_dir+'/summary_distances_confs.dat', delimiter = '\t',skip_header=1)

	order = data[omit:, 1].argsort()
	sorted_data = np.take(data, order, 0)

	desired_confs=["0 \n"]
	desired_confs_vis=["./CONF/conf0.gro "]
#	desired_confs_value[sorted_data[0,:]]
	print(sorted_data[:,1])
	
	track_k=0
	if min_dist==False:
		z=sorted_data[0,1]
	else:
		for v in range(len(sorted_data[:,0])):
			if sorted_data[v,1]>=min_dist and track_k==0:
				z=sorted_data[v,1]
				track_k=v
				print(track_k)
			else:
				pass

	desired_confs_values=[]
	for k in range(len(sorted_data[:,0])-1):
		if k>=track_k:
			if np.isnan(sorted_data[k+1,1])!=True:
				if np.abs(sorted_data[k+1,1]-z)>min_sep and sorted_data[k+1,1]<max_dist:
					#print(k+1,np.abs(sorted_data[k+1,1]-z))
					print(sorted_data[k+1,:])
					#print(np.abs(sorted_data[k+1,1]),"    ",k+1)
					z=sorted_data[k+1,1]
					desired_confs.append(str(k+1)+"\n")
					desired_confs_vis.append("./CONF/conf"+str(k+1)+".gro ")
					desired_confs_values.append(np.ndarray.tolist(sorted_data[k+1,:]))
	#desired_confs_values=[str(k) for k in desired_confs_values]
	print("maximum distance from initial positon = ",max_dist,"[nm]")
	print("minimum separation between windows = ",min_sep,"[nm]")
	print("number of windows = ",len(desired_confs))
	print("list of window confs = ",[int(l) for l in desired_confs])
	print("### Seleted confs ###\n",6*"t[ps] \t d[nm] \t ","\n",fmtcols(desired_confs_values,6))

#################################################################################
############################ EMAIL REPORT #######################################
	top_info="simulation title = \t"
	with open("system.top","r") as f3:
		topology = f3.readlines()
		for i,line in enumerate(topology):
			#print(i,line)
			if "[ system ]" in line:
				#print(i,line)
				top_info+=topology[i+1]


	with open("analysis/email_report.txt","w") as f3:
		f3.write(top_info)
		f3.write("directory = \t"+umbrella_dir+" \n")
		f3.write("computer name = \t"+socket.gethostname()+" \n")
		f3.write("Umbrella sampling:"+" \n")
		f3.write("maximum distance from initial positon = \t"+str(max_dist)+" [nm] \n")
		f3.write("minimum separation between windows = \t"+str(min_sep)+" [nm] \n")
		f3.write("number of windows = \t"+str(len(desired_confs))+"\n")
		f3.write("list of window confs = \t"+str([int(l) for l in desired_confs])+"\n")
#################################################################################

	# for k in range(len(data[:,0])-1):
	# 	z=data[0,1]
	# 	print(np.abs(data[k+1,1]-z)>0.2)
	# 	if np.abs(data[k+1,1]-z)>0.2:
	# 		z=data[k+1,1]
			#print(k,data[k,1])

	with open("confs.txt","w") as f2:
		f2.writelines(desired_confs)
		# *10 jer sum dist ispisuje u ps umjesto u broju conf...treba popraviti
		#f2.writelines([str(int(k)*10)+"\n" for k in desired_confs])

	with open("vmd_visconfs.sh","w") as f2:
		f2.write("vmd ")
		f2.writelines(desired_confs_vis)

# scale t[ps] |--> t[conf]
#	data_conf_ratio = data[-1,0]/num_files
#	print(num_files)
#	print(data[-1,0])
#	print(data_conf_ratio)


#Pr. data[:,0]*10 jer x[ps] |--> x[conf] simulacija traje 500ps ali ispisujemo 5000 frameova
	plt.plot(data[:,0],data[:,1],"-")
	for c in desired_confs:
		# *10 jer sum dist ispisuje u ps umjesto u broju conf...treba popraviti
		plt.axvline(float(c), linestyle='-',linewidth="0.05")
	plt.text(4, 3.8, "max_dist"+str(max_dist)+"nm"+"    min_sep"+str(min_sep))
	plt.xlabel("conf [number]")
	plt.ylabel("distance [nm]")
	plt.title("Distance between the ligand and coordinated solvent")
	plt.savefig("analysis/pull_distances.pdf")


def auto_omit():
	""" Fit sigmoid to data and determine up to which conf to omit all data. """
	pass

def fmtcols(mylist, cols):
	"""	Format columns. """
	lines = ("".join(str(mylist[i:i+cols])) for i in range(0,len(mylist),cols))
	return '\n'.join(lines)

if __name__ == '__main__':
	inputs()

