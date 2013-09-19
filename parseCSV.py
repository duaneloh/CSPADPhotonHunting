#!/usr/bin/env python
import csv
import collections as C

fn = "LCLS-cxi73013-Maia.csv"
cf = csv.DictReader(open(fn,'rU'), delimiter=',')
dataDict = C.OrderedDict({})
headers = C.OrderedDict([(a, None) for a in cf.fieldnames])

for line in cf:
	dataDict.setdefault(line['RunNum'], C.OrderedDict(line))

listWriter = csv.DictWriter(open(fn, "wb"), delimiter=',', fieldnames=headers, quotechar='|', quoting=csv.QUOTE_MINIMAL)

listWriter.writeheader()

for k, v in dataDict.iteritems():
	listWriter.writerow(v)