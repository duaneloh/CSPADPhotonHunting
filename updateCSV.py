#!/usr/bin/env python

import time
import csv
import collections as C

fn = "LCLS-cxi73013-Maia.csv"

class Found(Exception): pass

cf = csv.DictReader(open(fn,'rU'), delimiter=',')
dataDict = C.OrderedDict({})
priorityDict = C.OrderedDict({})
headers = C.OrderedDict([(a, None) for a in cf.fieldnames])
for line in cf:
	priorityDict.setdefault(line['Priority'], [])
	priorityDict[line['Priority']].append(line)
	dataDict.setdefault(line['RunNum'], C.OrderedDict(line))

sortedKeys = priorityDict.keys()
sortedKeys.sort()
try:
	for k in reversed(sortedKeys):
		for l in priorityDict[k]:
			if l['TestField'] == '0':
				dataDict[l['RunNum']]['TestField'] = '1'
				raise Found
except Found:
	pass

listWriter = csv.DictWriter(open(fn, "wb"), delimiter=',', fieldnames=headers, quotechar='|', quoting=csv.QUOTE_MINIMAL)
listWriter.writeheader()
for k, v in dataDict.iteritems():
	listWriter.writerow(v)
