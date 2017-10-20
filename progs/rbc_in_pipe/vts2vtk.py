import sys
from paraview.simple import *

#if not servermanager.ActiveConnection:connection = servermanager.Connect()

if len(sys.argv) < 2:
	print "Useage\n\tpvpython vts2vtk.py FILE.vtk [more FILES.vtk...]\n"

for arg in sys.argv[1:]:
	reader = OpenDataFile((arg + ".vts"))
	writer = servermanager.writers.DataSetWriter(Input=reader,FileName=(arg + ".vtk"))
	writer.UpdatePipeline()
