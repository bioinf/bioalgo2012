#! /usr/bin/env python

import sys
import re

def getPoint(file, mode):
	point = []
	linesCount = 0
	if mode == "test":
		linesCount = 2
	elif mode == "train":
		linesCount = 3
	else:
		return []
	for i in range(0,linesCount):
		line = file.readline()
		if not line:
			point = []
			break
		point.append(line)
	return point

def getPoints(file, mode):
	points = []
	trainingsetFile = open(file, 'r')
	vectors = []
	while True:
		point = getPoint(trainingsetFile, mode)
		if point == []:
			break
		points.append(point)
	trainingsetFile.close()
	return points

def getChromosomeNumber(string):
	res = string.split(" ")
	res = re.split('chr', res[1])
	if 'X' in res[1]:
		res[1] = str(23)
	if 'Y' in res[1]:
		res[1] = str(24)
	return int(res[1])

acgtFilter = ['A','C','G','T','N']

def getAcgtCount(string):
	vector = []
	for i in range(0,len(acgtFilter)):
#		print acgtFilter[i] + " : " + str(string.count(acgtFilter[i]))
		vector.append(string.count(acgtFilter[i]))
	return vector

def getAcgtFract(string):
	result = []
	vector = getAcgtCount(string)
	meaningfullStrLen = len(string) - vector[4]
	for i in vector:
		if meaningfullStrLen != 0:
			result.append(float(i)/meaningfullStrLen)
		else:
			result.append(0)
	result[4] = float(vector[4]) / len(string)
	return result

def getGcFract(string):
	vector = getAcgtCount(string)
	return float(vector[1] + vector[2])/(len(string) - vector[4])

nucleotides = ['A','C','G','T']

def formKmerTemplate(k):
	templates = []
	if k == 1:
		for nucl in nucleotides:
			templates.append(nucl)
	else:
		for nucl in nucleotides:
			newTemplates = []
			kmers = formKmerTemplate(k-1)
			for kmer in kmers:
				newTemplates.append(nucl + kmer)
			templates += newTemplates
	return templates

def getKmerCount(string,k):
	templates = []
	templates = formKmerTemplate(k)
	kmers = []
	for template in templates:
		kmers.append(string.count(template))
#	print kmers
	return kmers

def getKmerCountNorm(string,k):
	kmers = getKmerCount(string,k)
	result = []
	for x in kmers:
		result.append(float(x)/len(string))
	return result

def getFeatures(descriptor, sequence):
	vector = [len(sequence)] \
	+ [getChromosomeNumber(descriptor)] \
	+ getAcgtFract(sequence) \
	+ getKmerCountNorm(sequence, 2)
##	+ [getGcFract(sequence)] \
#	print vector
#	vector = getKmerCount(sequence, 2)
	return vector

def generateTrainVectorsFile(file):
	points = getPoints(file, "train")
	vectors = []
	for point in points:
		vectors.append([int(point[1])] + getFeatures(point[0], point[2]))
	vectorsFile = open("vectors/" + file.split('/')[-1], 'w')
	for vector in vectors:
		vectorsFile.write(str(vector[0]))
		for i in range(1,len(vector)):
			vectorsFile.write(" " + str(vector[i]))
		vectorsFile.write('\n')
	vectorsFile.close()

