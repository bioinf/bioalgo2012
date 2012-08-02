import os
import sys
import csv

reader = open(sys.argv[1])
user_ids = set()
for line in reader:
    user_ids.add(int(line[4:]))       

pheno_reader = csv.reader(open(sys.argv[2]), delimiter=";")

nearest = []
#skipping header
names = pheno_reader.next()
for i in  pheno_reader:
    print i[0]
    if int(i[0]) in user_ids:
        print "Found ", i[0]
        nearest.append(i)


l = len(nearest[0])


for i in range(1, l):
    stats = dict()
    for rec in nearest:
        if rec[i] in stats:
            stats[rec[i]] += 1
        else:
            stats[rec[i]] = 1
    print "\n"
    print names[i]
    for k,v in stats.iteritems():
        print k, v

print user_ids
