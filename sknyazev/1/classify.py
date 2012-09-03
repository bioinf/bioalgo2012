#! /usr/bin/env python
from multiprocessing import Pool
import commands
import sys
import os

def f(arg):
	print commands.getoutput(arg)

tray = []
resultPath = "result/"
launchCommand = ["./classifier.py \"vectors\" \"testData/0\" \"result/\""]
modes = ["svm", "neighbors", "tree", "naive"]

for subdir, dirs, files in os.walk("vectors/"):
	for file in files:
		for mode in modes:
			pathName = resultPath + file + "_" + mode + "/"
			commands.getoutput("mkdir " + pathName)
			tray.append("./classifier.py vectors/" + file +" testData/0 " + pathName + " " + mode)

p = Pool(16)
p.map(f, tray)
