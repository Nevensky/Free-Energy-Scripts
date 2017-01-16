#!/usr/bin/env python3
import os
from subprocess import call
import numpy as np


def parallel(x):
	if x<ncpus:
		return "&"
	else:
		return "\nwait"


def gen_distances(path,ncpus):
	""" Generates and executes gmx distance input files for all confs in /CONF folder. """
	filenames = [item for item in next(os.walk(path))[2] if ".gro" in item and "#" not in item]
	num_files = len(filenames)
	with open(path+"/run_gendistances.sh","w") as f:
		f.write("#!/usr/bin/env bash\n")
		# prallelization counter: "&" * ncpus then "wait"
		x=0
		for k in range(num_files):
			if x<ncpus:
				x+=1
			else:
				x=0
			commandstring="gmx distance -s "+umbrella_dir+"/Production_PULL/md_pull.tpr -f "+path+"/conf{0}.gro -n "+umbrella_dir+"/index.ndx -oall "+path+"/dist{1}.xvg -select \'com of group \"DHK\" plus com of group \"COORDSOLVENT\"\' {2}\n"
			f.write(commandstring.format(str(k),str(k),parallel(x)))

	call(["chmod","+x",path+"/run_gendistances.sh"])
	print("Running gmx distance")
	call([path+"/run_gendistances.sh",">/dev/null"])
#	os.system(path+"/run_gendistances.sh") #deprecated but works with one argument


def write_distances(path):
	""" Concatenates .xvg files generate by gmx distance and writes them to a file."""
	filenames = [item for item in next(os.walk(path))[2] if ".xvg" in item and "#" not in item]
	num_files = len(filenames)
	distances = []
	for k in range(num_files):
		with open(path+"/dist{0}.xvg".format(str(k))) as f:
			dist_xvgs=f.readlines()
		for i,line in enumerate(dist_xvgs):
			if "@TYPE xy" in line:
				distances.append(dist_xvgs[i+1].split())

	distances=np.asarray(distances,dtype=np.float64)
	np.savetxt("summary_distances.dat", distances, fmt='%1.3f', delimiter='\t', newline='\n',header='t[ps]  d[nm]')

	# ruzno, prepraviti nekako
	distances_inconfs = list(enumerate(distances[:,1]))
	distances_inconfs = np.asarray([[i,j] for i,j in distances_inconfs])
	np.savetxt("summary_distances_confs.dat", distances_inconfs, fmt=['%i','%1.3f'], delimiter='\t', newline='\n',header='conf  d[nm]')

if __name__ == '__main__':
	ncpus = 8
	umbrella_dir = os.getcwd()
	path = umbrella_dir+"/CONF"
	gen_distances(path,ncpus)
	write_distances(path)
