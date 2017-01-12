#!/usr/bin/env python

import os
from subprocess import call


dir=os.getcwd()

k=[]
with open(dir+"/confs.txt","r") as f:
	k = f.readlines()

z=["vmd "]
for i in k:
	z.append("".join([str(dir),"/conf",i.rstrip(),".gro "]))

#print("".join(z))

with open(dir+"/vmd_visconfs.sh","w") as f:
	f.writelines(z)


call(["chmod","+x","vmd_visconfs.sh"])
