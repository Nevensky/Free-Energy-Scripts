#!/usr/bin/env python
import click


@click.command()

@click.option('-conf',default='./confs.txt',help='Array of configurations. Defaults to import numbers in ./confs.txt')
 
def inputs(conf):
	'''Grompp and mdrun script generator for Umbrella sampling. \n

Developed by Neven Golenic | neven.golenic@gmail.com'''
	if ".txt" in conf:
		with open("./confs.txt","r") as file:
			confs=file.readlines()
			confs=[int(c) for c in confs]
			print(confs)
	else:
		confs=conf
		confs=confs.split(" ")

	out="mkdir UMBRELLA_NPT\n"+"cd UMBRELLA_NPT\n"
	out2="mkdir UMBRELLA_PROD\n"+"cd UMBRELLA_PROD\n"

	#create grompps
	for k in confs:
		s="gmx grompp -f $MDP/npt_umbrella.mdp -c $CONF/conf{0}.gro -p $UMBRELLA_HOME/system.top -n $UMBRELLA_HOME/index.ndx -o umbrella_npt{1}.tpr -maxwarn 1\n"
		s=s.format(str(k),str(k))
		out+=s
		s2="gmx grompp -f $MDP/md_umbrella.mdp -c $UMBRELLA_HOME/UMBRELLA_NPT/umbrella_npt{0}.gro -t $UMBRELLA_HOME/UMBRELLA_NPT/umbrella_npt{1}.cpt -p $UMBRELLA_HOME/system.top -n $UMBRELLA_HOME/index.ndx -o umbrella_prod{2}.tpr\n"
		s2=s2.format(str(k),str(k),str(k))
		out2+=s2

	out+="\n \n \n"
	out2+="\n \n \n"

	#create mdruns
	for k in confs:
		s3="gmx mdrun -nt $NCPUS -pin auto -nb gpu -deffnm umbrella_npt{0}\n"
		s3=s3.format(str(k))
		out+=s3
		s4="gmx mdrun -nt $NCPUS -pin auto -nb gpu -deffnm umbrella_prod{0} -pf pullf-umbrella_prod{1}.xvg -px pullx-umbrella_prod{2}.xvg\n"
		s4=s4.format(str(k),str(k),str(k))
		out2+=s4

	out3,out4="",""
	#create list of tpr files and pullf files
	for k in confs:
		out3+="umbrella_prod{0}.tpr\n".format(str(k))
		out4+="pullf-umbrella_prod{0}.xvg\n".format(str(k))

	#out
	with open("exec_umbrella_npt.sh","w") as f:
		f.writelines(out)

	#out2
	with open("exec_umbrella_prod.sh","w") as f:
		f.writelines(out2)

	#out3
	with open("tpr-files.dat","w") as f:
		f.writelines(out3)

	#out4
	with open("pullf-files.dat","w") as f:
		f.writelines(out4)
 
if __name__ == '__main__':
	inputs()
