#!/usr/bin/env python

#script to make mapfile

import os
import sys

number = range(51)
number = number[1:51]

lineageDef = {}
for key in [43, 44, 45]:
    lineageDef[key] = 'atricapillus'
for key in [4, 6, 8, 13, 14, 18, 19, 20, 23, 24, 25, 34, 35, 36, 46, 50]:
    lineageDef[key] = 'central'
for key in [42]:
    lineageDef[key] = 'lateralis'
for key in [7, 10, 11, 12, 22, 27, 28, 33, 47, 48, 49]:
    lineageDef[key] = 'northern'
for key in [9, 31, 32]:
    lineageDef[key] = 'variegatus'
for key in [1, 2, 3, 5, 15, 16, 17, 21, 26, 29, 30, 37, 38, 39, 40, 41]:
    lineageDef[key] = 'southern'

outfile = 'mapMitofile.txt'
out = open(outfile, 'w')

for num in number:
	
	element0 = 'index' + str(num) + '\t'
	element1 = lineageDef[num] + '\n'
	out.write(element0)
	out.write(element1)

out.close()
