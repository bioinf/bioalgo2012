#! /usr/bin/env python

import sys
import commands
from sklearn import svm
from sklearn import cross_validation
from sklearn import neighbors
from sklearn import tree
from sklearn.naive_bayes import GaussianNB
from featureExtractor import *

if(len(sys.argv) != 5):
	print("Usage: " + sys.argv[0] + " trainingsetVectorsFile testFile resultDir mode")
	exit()

vectorsFile = open(sys.argv[1], 'r')
X = []
Y = []
for line in vectorsFile:
	vector = line.split(" ")
	Y.append(int(vector[0]))
	xvector = []
	for i in range(1,len(vector)):
		xvector.append(float(vector[i]))
	X.append(xvector)
X_train, X_test, y_train, y_test = cross_validation.train_test_split(X, Y, test_size = 0.3, random_state = 0)

if (sys.argv[4] == "svm"):
	clf = svm.SVC()
elif (sys.argv[4] == "neighbors"):
	clf = neighbors.KNeighborsClassifier(1)
elif (sys.argv[4] == "tree"):
	clf = tree.DecisionTreeClassifier()
elif (sys.argv[4] == "naive"):
	clf = GaussianNB()
else:
	exit()
print "Learning: ..."
print "Test points extracting: ..."
#clf.fit(X,Y)
clf.fit(X_train, y_train)
points = getPoints(sys.argv[2], "test")
if sys.argv[3][-1] != '/':
	sys.argv[3].append('/')
print "Predicting: ..."
commands.getoutput("rm " + sys.argv[3] + "*")
currenFileName = "chr1"
file = open(sys.argv[3] + currenFileName, "a")
print currenFileName
for point in points:
	vector = getFeatures(point[0], point[1])
	descriptor = point[0].split(" ")
	if currenFileName != descriptor[1]:
		file.close()
		currenFileName = descriptor[1]
		print currenFileName
		file = open(sys.argv[3] + currenFileName, "a")
	file.write(str(int(descriptor[2])) + " " + str(int(descriptor[3])) + " " + str(int(clf.predict(vector)[0])) + "\n") 
file.close()
