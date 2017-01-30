#!/usr/bin/env python3

import os
from subprocess import call

print("Installing dependencies...")
call(["conda","install","numpy"])
call(["conda","install","pip"])
call(["pip","install","termcolor"])
call(["conda","install","scipy"])
call(["conda","install","matplotlib"])
call(["conda","install","click"])

call(["chmod","+x",os.getcwd()+"/bar_main.py"])
call(["rm","-rf",os.getcwd()+"/bin/genbar"])
call(["ln","-s",os.getcwd()+"/bar_main.py",os.getcwd()+"/bin/genbar"])

home_dir = os.getenv("HOME")
with open(home_dir+"/.bashrc","a") as f:
  f.write("\n#BAR automation script\n"+"export PATH=$PATH:"+os.getcwd()+"/bin/genbar")

print("Command genbar added to .bashrc\n Instalation successful.")