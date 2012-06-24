import os
import re
import sys

if len(sys.argv) != 3:
    print "Usage: <phenotype csv file> <user list>"
    exit(0)

user = []
userListFile = open(sys.argv[2])
for line in userListFile:    
    num =  re.findall('\d+', line)
    if len(num) > 0:
        user.append(num[0])
userListFile.close()

phenotype = []
phenotypeFile = open(sys.argv[1])
header = phenotypeFile.readline().strip().split(';')
for line in phenotypeFile:
    line = line.strip().split(';')
    if line[0] in user:
        phenotype.append(line)
phenotypeFile.close()
for i in range(1, len(header)):
    print header[i], ':'
    trait = dict()
    for record in phenotype:
        nowTrait = record[i]
        if nowTrait not in trait:
            trait[nowTrait] = 0
        trait[nowTrait] += 1
    trait = [(trait[key], key) for key in trait]
    trait.sort(reverse=True)
    for t in trait:
        print t
