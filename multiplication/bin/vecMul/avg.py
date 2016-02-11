#!/usr/bin/python
import numpy as np
from locale import atof

size = list(range(10,1010,10))
timeList = []
flopsList = []
time = []
flops = []

for i in range(1,11):
	filename = "time_" + str(i) + ".txt"
	pf = open(filename, "r")
	for line in pf:
		content = line.split()
		timeList.append(atof(content[1]))
		flopsList.append(atof(content[2]))
	pf.close()
	time.append(timeList)
	flops.append(flopsList)

timeArray = np.array(time)
flopsArray = np.array(flops)
timeAvg = np.sum(timeArray, axis=0)
flopsAvg = np.sum(flopsArray, axis=0)

pf = open("avgVecMul.txt", "a")
for i in size:
	if timeAvg[size.index(i)] > 1e-16 and flopsAvg[size.index(i)] != 'inf':	
		pf.write("%s %s %s\n" % (i, timeAvg[size.index(i)], flopsAvg[size.index(i)]))
pf.close()



