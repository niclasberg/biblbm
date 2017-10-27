import numpy as np
import glob
import matplotlib.pyplot as plt

folders = ['tweezers_G6e-06_Kb4e-18_dx1e-06', 'tweezers_G6e-06_Kb4e-18_dx5e-07', 'tweezers_G6e-06_Kb4e-18_dx2.5e-07']
legends = ['$\Delta x = 1\ \mu m$', '$\Delta x = 0.5\ \mu m$', '$\Delta x = 0.25\ \mu m$']
styles = ['r:', 'g.-', 'k--']

# Read experimental data
expData = np.genfromtxt('mills_optical_tweezer_final.txt', delimiter=',')
expForce = expData[:, 0]
expDx = expData[:, 3]
expDxErr = expData[:, 4] - expData[:, 3]
expDyz = expData[:, 1]
expDyzErr = expData[:, 1] - expData[:, 2]

cmpPlot, ax = plt.subplots(figsize=(6,6))
plt.errorbar(expForce, expDx, yerr=expDxErr, fmt='sb', label='Experiments')
plt.errorbar(expForce, expDyz, yerr=expDyzErr,fmt='sb')

# Read simulation data
for folder, legend, linestyle in zip(folders, legends, styles):
	files = glob.glob(folder +'/deformation*.txt')

	data = None
	for i, f in enumerate(files):
		tmp = np.genfromtxt(f, delimiter=',')
		if tmp.size > 0:
			if not data:
				data = tmp
			else:
				data = np.append(data, tmp, axis = 0)

	data = np.array(sorted(data, key = lambda x: x[0]))

	s = expDx[0] / (1e6*data[0, 2])
	it = data[:, 0]
	force = data[:, 1]
	Dx = data[:, 2]
	Dyz = data[:, 3]

	h = plt.figure()
	plt.plot(it, 1e6*Dx, label='Dx')
	plt.plot(it, 1e6*Dyz, label='Dyz')
	plt.xlabel('Iteration [-]')
	plt.ylabel('Deformation [$\mu$m]')
	plt.savefig(folder+'/deformationConvergence.png', bbox_inches='tight')
	plt.close(h)

	# Find where the force changes 
	inds = np.where(np.diff(force) > 0)[0]
	inds = np.concatenate(([0], inds))

	simDx = 1e6*Dx[inds]
	simDyz = 1e6*Dyz[inds]
	simForce = 1e12*force[inds]
	simForce[0] = 0.
	
	plt.figure(cmpPlot.number)
	plt.plot(simForce, simDx, linestyle, label=legend)
	plt.plot(simForce, simDyz, linestyle)
plt.legend(loc='upper left')
#xmin,xmax = ax.get_xlim()
#ymin,ymax = ax.get_ylim()
#asp = abs((xmax-xmin)/(ymax-ymin))
#ax.set_aspect(asp)

plt.ylim([0, 20])
plt.xlim([0, 200])
plt.xlabel('Force [pN]', fontsize=16)
plt.ylabel('Axis length [$\mu$m]', fontsize=16)
plt.savefig('deformationComparison.pdf', bbox_inches='tight')
