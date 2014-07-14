#!/usr/bin/python3

import subprocess
import os
import sys

os.system("git add --all .")
subprocess.call(["git", "commit", "-m", sys.argv[1]])
os.system("git push origin master")
