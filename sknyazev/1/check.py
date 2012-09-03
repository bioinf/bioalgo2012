#! /usr/bin/env python

import os
import commands
import sys

results = commands.getoutput("ls result/").split('\n')
os.chdir("/storage/labnas/bioalgo")

for result in results:
	print result + ":"
	print commands.getoutput("./check ../students/sknyazev/bioalg/1/result/" + result)
	print ""
