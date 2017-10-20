import numpy as np
import glob
import matplotlib.pyplot as plt

folders = ['tweezers_G6e-06_Kb4e-18_dx1e-06', 'tweezers_G6e-06_Kb4e-18_dx5e-07', 'tweezers_G6e-06_Kb4e-18_dx3e-07', 'tweezers_G6e-06_Kb4e-18_dx2.5e-07']

# Read experimental data
expData = np.genfromtxt('mills_optical_tweezer_final.txt', delimiter=',')
expForce = expData[:, 0]
expDx = expData[:, 3]
expDxErr = expData[:, 4] - expData[:, 3]
expDyz = expData[:, 1]
expDyzErr = expData[:, 1] - expData[:, 2]

cmpPlot = plt.figure()
plt.errorbar(expForce, expDx, yerr=expDxErr, fmt='sb')
plt.errorbar(expForce, expDyz, yerr=expDyzErr,fmt='sr')

# Read simulation data
for folder in folders:
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
	plt.plot(simForce, simDx, '.-', label='Dx')
	plt.plot(simForce, simDyz, '.-', label='Dyz')
plt.ylim([0, 18])
plt.xlim([0, 200])
plt.xlabel('Force [pN]')
plt.ylabel('Deformation [$\mu$m]')
plt.savefig('deformationComparison.png', bbox_inches='tight')
