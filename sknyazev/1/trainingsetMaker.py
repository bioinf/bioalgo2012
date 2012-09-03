#! /usr/bin/env python

import sys
import os
import math
import random

def printUsage():
	print "Usage: " + sys.argv[0] + " datasetPath pointsCount accuracy"

def getDatasetCount(fileName):
	datasetCount = 0
	datasetFile = open(fileName, 'r')
	for line in datasetFile:
		if '>' in line:
			datasetCount += 1
	datasetFile.close()
	return datasetCount

def getTrainingsetCountCoef(x):
	return math.log(x) ** 5

def getTrainingsetSize(classes, trainingsetThreshold):
	result = {}
	for i in classes:
		result[i] = getTrainingsetCountCoef(classes[i])
	maxCoef = max(result.viewvalues())
	for i in result:
		result[i] = int(float(result[i])/maxCoef * trainingsetThreshold)
	return result

def getMeaningfullClassesSize(datasetDir, fractionThreshold):
	total = 0
	datasetCount = {}
	for subdir, dirs, files in os.walk(datasetDir):
		for file in files:
			count = getDatasetCount(datasetDir + file)
			datasetCount[int(file)] = count
			total += count
	i = 1
	covarage = 0
	classes = 0
## meaningfullClasses contain members count
	meaningfullClasses = {}
	for cl in datasetCount:
		portion = float(datasetCount[cl]) * 100 / total
		if portion >= fractionThreshold:
			covarage += portion
			classes += 1
			meaningfullClasses[cl] = datasetCount[cl]
##		print(str(i) + " " + str(x) + " " + str(portion) + "% training set:" + \
##		str(getTrainingsetCountCoef(x)))
		i += 1
##	print("total " + str(total))
	print("Threshold: " + str(fractionThreshold) + "%, coverage: " + \
	str(covarage) + ", classes: " + str(classes))
	return meaningfullClasses

def getPoint(file):
	point = []
	string = ""
	while('>' not in string):
		string = file.readline()
	point.append(string)
	point.append(file.readline())
	return point

def getTrainingset(path, meaningfullClasses, trainingSetCount):
	classesPointsNumbers = {}
	for cl in meaningfullClasses:
		dataSet = random.sample(xrange(1,meaningfullClasses[cl]-1,1),trainingSetCount[cl])
		dataSet.sort()
		classesPointsNumbers[cl] = dataSet
	classesPoints = {}
	for cl in classesPointsNumbers:
		print str(cl)
		stringNumber = 0
		points = []
		file = open(path + str(cl), 'r')
		for i in classesPointsNumbers[cl]:
			while(stringNumber - 1 != i):
				getPoint(file)
				stringNumber += 1
			points.append(getPoint(file))
#			print (str(cl) + " : " + str(i))
			stringNumber += 1
		classesPoints[cl] = points
	return classesPoints
		

if len(sys.argv) != 4:
	printUsage();
	exit()
meaningfullClasses = getMeaningfullClassesSize(sys.argv[1], float(sys.argv[3]))
trainingsetSize = getTrainingsetSize(meaningfullClasses, int(sys.argv[2]))

print("Meaningfull classes:")
print(meaningfullClasses)
print("Training set size:")
print(trainingsetSize)

print("Training set generating: ...")
trainingset = getTrainingset(sys.argv[1], meaningfullClasses, trainingsetSize)
print("output file writing: ...")

trainingsetFilePath = "trainingset/trainingset_" + sys.argv[2] + "_" + sys.argv[3]
trainingsetFile = open(trainingsetFilePath, 'w')

for cl in trainingset:
	for p in trainingset[cl]:
		trainingsetFile.write(p[0] + str(cl) + "\n" + p[1])

trainingsetFile.close()
