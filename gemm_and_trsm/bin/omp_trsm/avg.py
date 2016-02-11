#!/usr/bin/python
import numpy as np
from locale import atof

size = list(range(100,3100,100))
timeList = []
time = []

for i in range(1,11):
	filename = "time_" + str(i) + ".txt"
	pf = open(filename, "r")
	for line in pf:
		content = line.split()
		timeList.append(atof(content[1]))
	pf.close()
	time.append(timeList)

timeArray = np.array(time)
timeAvg = np.sum(timeArray, axis=0) / len(size)

pf = open("avgTime.txt", "a")
for i in size:
	pf.write("%s %s\n" % (i, timeAvg[size.index(i)]))
pf.close()

