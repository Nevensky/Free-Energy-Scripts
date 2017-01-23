#!/usr/bin/env python3

import os
from subprocess import call

call(["chmod","+x",os.getcwd()+"/bar_main.py"])
call(["rm","-rf",os.getcwd()+"/bin/genbar"])
call(["ln","-s",os.getcwd()+"/bar_main.py",os.getcwd()+"/bin/genbar"])

home_dir = os.getenv("HOME")
with open(home_dir+"/.bashrc","a") as f:
  f.write("\n#BAR automation script\n"+"export PATH=$PATH:"+os.getcwd()+"/bin/genbar")