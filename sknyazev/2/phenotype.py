#! /usr/bin/env python

import csv
import sys
import re

def getRelatives(file, n):
	relatives = []
	f = open(file,'r')
	for line in f:
		if n == 0:
			break
		n -= 1
		t = line.split('\t')
		relatives.append(t[0])
	f.close()
	return relatives

def findMax(matrix, column):
	vector = [element[column] for element in matrix]
	return max((set(vector)).difference(set(["-"])), key=vector.count)

def getPhenotype(relatives, path):
	
	names = []
	for relative in relatives:
		name = relative.split("_")
		names.append(re.findall( r'\d+', name[1]))

	allPhenotypes = csv.reader(open(path,"rb"), delimiter=";")
	phenotypes = []
	for line in allPhenotypes:
		ID = re.findall(r"^\d+$", line[0])
		if ID and ID in names:
			phenotypes.append(line)	
	phenotype = []
	for i in range(3,14):
		phenotype.append(findMax(phenotypes,i))
	return phenotype

relations = getRelatives(sys.argv[1], 10)
f = open("myPhenotype",'w')
phenotype = getPhenotype(relations,"/storage/labnas/students/slebedev/opensnp/raw/phenotypes_201207241103.csv")
for i in phenotype:
	f.write(str(i) + "\t")
f.write("\n")
f.close()
