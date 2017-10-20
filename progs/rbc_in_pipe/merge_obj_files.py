#!/usr/bin/python

import sys
import os
import glob
from collections import defaultdict

if len(sys.argv) < 2:
	print 'Usage: merge_obj_files folder'
folder = str(sys.argv[1])
if folder[len(folder) - 1] != '/':
	folder = folder + '/'

print 'Merging files in ' + folder

# List files
files = glob.glob(folder + 'objsP*.txt')

# Iteration => file map
md = defaultdict(list)

for file in files:
	# Read current iteration
	with open(file, 'r') as f:
		first_line = f.readline()
		result = first_line.replace(':', ' ').split()
		if result[0] != 'iteration':
			print 'Expected first line of ' + file + ' to be of the form iteration:%d'
		md[int(result[1])].append(file)

# Read files and merge them
for iteration,files in md.iteritems():
	fname = folder + 'objects_' + str(iteration) + '.txt'
	with open(fname, 'w') as fout:
		print 'Writing to file ' + fname 
		first_file = True
		for file in files:
			print ' Reading ' + file
			with open(file, 'r') as fin:
				lines = fin.readlines()
				lines = [l for l in lines if first_file or not 'iteration:' in l]
				fout.writelines(lines)
			first_file = False
			print ' Removing ' + file
			#os.remove(file)
