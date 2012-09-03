#! /usr/bin/env python
from featureExtractor import *
import sys
import os

if len(sys.argv) < 2:
	print ("Usage: " + sys.argv[0] + "trainingsetDir")

print ("Vector file generating ...")
for subdir, dirs, files in os.walk(sys.argv[1]):
	for file in files:
##		print(file)
		generateTrainVectorsFile(sys.argv[1] + file)
