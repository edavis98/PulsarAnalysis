# -----------------------------------------------------------------------------
# File: relativistic_trajectory.py
# Name: Evan Davis
# Description: Creates 2D visualizations of relativistic trajectory
# and gamma evolution over time of a particle in varying magnetic/electic
# field configurations
# Usage: External file allprt01_125133.dat is needed to run this program.
# One may use a os.chdir command to change directory to where the .dat file is
# located if the program does not detect the .dat file
# -----------------------------------------------------------------------------

# %% Load libraries

import os
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# %% Read in data

statinfo = os.stat('allprt01_125133.dat')  # info from data file
sfile = statinfo.st_size  # size of data file
nattr = 61  # number of attributes per particles

# determine number of particles by dividing size of file by each element
# with 8 bytes per element
numprt = math.floor(sfile / (nattr*8))

# create array of data (stored as one long row listing all 61 attributes per
# particle in order)
with open('allprt01_125133.dat', 'rb') as ff:
    data = np.fromfile(ff, dtype=np.float64, count=nattr*numprt)

# reshape data into array  with each row being a particle with 61
# matching attribute columns
datars = np.reshape(data, (numprt, nattr))

# define constants (charge, mass, speed of light, time span)
q = -4.8E-10
m = 9.1E-28
c = 3E10
init = (0, 0, 0, 10, 10, 0)
t = np.linspace(0, 0.0001, 5000)

# %% ODE to be used for gamma evolution over time


# define ODE
def deri(yvec, t, qmc, c, Ex, Ey, Ez, Bx, By, Bz):
    gamma = np.sqrt(1 + yvec[3]**2 + yvec[4]**2 + yvec[5]**2)
    return (yvec[3] * c / gamma, yvec[4] * c / gamma, yvec[5] * c / gamma,
            qmc * (Ex + (yvec[4] * Bz - yvec[5] * By) / gamma),
            qmc * (Ey + (yvec[5] * Bx - yvec[3] * Bz) / gamma),
            qmc * (Ez + (yvec[3] * By - yvec[4] * Bx) / gamma))


# %% Create plot showing gamma evolution and trajectory over time with Ey=0.0

# initialize electric and magnetic field strengths
Bx = 1
By = 0
Bz = 0
Ex = 0
Ey = 0.0*Bx
Ez = 0
qmc = q/(m*c)

# integrate ODE to get fucntion
yarr = odeint(deri, init, t, args=(qmc, c, Ex, Ey, Ez, Bx, By, Bz))

# plot relativistic trajectory
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(yarr[:, 0], yarr[:, 1], yarr[:, 2], marker='.')
plt.title('Relativistic Trajectory of Particle with Bx = 1.0')
plt.savefig('rel_traj_Bx_1.png', dip=300)
plt.show()

# plot gamma evolution
fig = plt.figure()
gtlist = []
time = range(len(yarr))

# calculate gamma for plot
for i in time:
    gammat = np.sqrt(yarr[i, 3]**2 + yarr[i, 4]**2 + yarr[i, 5]**2)
    gtlist.append(gammat)

plt.plot(time, gtlist)
plt.title('Gamma Evolution With Time for Bx = 1.0')
plt.xlabel('Time (seconds)')
plt.ylabel('Gamma Value')
plt.savefig('gamma_t_Bx_1.png', dip=300)
plt.show()

# %% Create plot showing gamma evolution and trajectory over time with Ey=0.8

# initialize electric and magnetic field strengths
Bx = 1
By = 0
Bz = 0
Ex = 0
Ey = 0.8*Bx
Ez = 0
qmc = q/(m*c)

# integrate ODE to get fucntion
yarr = odeint(deri, init, t, args=(qmc, c, Ex, Ey, Ez, Bx, By, Bz))

# plot relativistic trajectory
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(yarr[:, 0], yarr[:, 1], yarr[:, 2], marker='.')
plt.title('Relativistic Trajectory of Particle with Bx = 1.0, Ey = 0.8')
plt.savefig('rel_traj_Bx_1_Ey_0.8.png', dip=300)
plt.show()

# plot gamma evolution
fig = plt.figure()
gtlist = []
time = range(len(yarr))

# calculate gamma for plot
for i in time:
    gammat = np.sqrt(yarr[i, 3]**2 + yarr[i, 4]**2 + yarr[i, 5]**2)
    gtlist.append(gammat)

plt.plot(time, gtlist)
plt.title('Gamma Evolution With Time for Bx = 1.0, Ey = 0.8')
plt.xlabel('Time (seconds)')
plt.ylabel('Gamma Value')
plt.savefig('gamma_t_Bx_1_Ey_0.8.png', dip=300)
plt.show()

# %% Create plot showing gamma evolution and trajectory over time with Ey=1.2

# initialize electric and magnetic field strengths
Bx = 1
By = 0
Bz = 0
Ex = 0
Ey = 1.2*Bx
Ez = 0
qmc = q/(m*c)

# integrate ODE to get fucntion
yarr = odeint(deri, init, t, args=(qmc, c, Ex, Ey, Ez, Bx, By, Bz))

# plot relativistic trajectory
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(yarr[:, 0], yarr[:, 1], yarr[:, 2], marker='.')
plt.title('Relativistic Trajectory of Particle with Bx = 1.0, Ey = 1.2')
plt.savefig('rel_traj_Bx_1_Ey_1.2.png', dip=300)
plt.show()

# plot gamma evolution
fig = plt.figure()
gtlist = []
time = range(len(yarr))

# calculate gamma for plot
for i in time:
    gammat = np.sqrt(yarr[i, 3]**2 + yarr[i, 4]**2 + yarr[i, 5]**2)
    gtlist.append(gammat)

plt.plot(time, gtlist)
plt.title('Gamma Evolution With Time for Bx = 1.0, Ey = 1.2')
plt.xlabel('Time (seconds)')
plt.ylabel('Gamma Value')
plt.savefig('gamma_t_Bx_1_Ey_1.2.png', dip=300)
plt.show()

# %% Create plot showing gamma evolution and trajectory over time with
# Ey=0.8, Ex=0.1

# initialize electric and magnetic field strengths
Bx = 1
By = 0
Bz = 0
Ex = 0.1
Ey = 0.8*Bx
Ez = 0
qmc = q/(m*c)

# integrate ODE to get fucntion
yarr = odeint(deri, init, t, args=(qmc, c, Ex, Ey, Ez, Bx, By, Bz))

# plot relativistic trajectory
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(yarr[:, 0], yarr[:, 1], yarr[:, 2], marker='.')
plt.title('Relativistic Trajectory of Particle with Bx = 1.0, Ey = 0.8,'
          'Ex = 0.1')
plt.savefig('rel_traj_Bx_1_Ey_0.8_Ex_0.1.png', dip=300)
plt.show()

# plot gamma evolution
fig = plt.figure()
gtlist = []
time = range(len(yarr))

# calculate gamma for plot
for i in time:
    gammat = np.sqrt(yarr[i, 3]**2 + yarr[i, 4]**2 + yarr[i, 5]**2)
    gtlist.append(gammat)

plt.plot(time, gtlist)
plt.title('Gamma Evolution With Time for Bx = 1.0, Ey = 0.8, Ex = 0.1')
plt.xlabel('Time (seconds)')
plt.ylabel('Gamma Value')
plt.savefig('gamma_t_Bx_1_Ey_0.8_Ex_0.1.png', dip=300)
plt.show()

# %% Create plot showing gamma evolution and trajectory over time with
# Ey=0.8, Ex=0.01

# initialize electric and magnetic field strengths
Bx = 1
By = 0
Bz = 0
Ex = 0.01
Ey = 0.8*Bx
Ez = 0
qmc = q/(m*c)

# integrate ODE to get fucntion
yarr = odeint(deri, init, t, args=(qmc, c, Ex, Ey, Ez, Bx, By, Bz))

# plot relativistic trajectory
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(yarr[:, 0], yarr[:, 1], yarr[:, 2], marker='.')
plt.title('Relativistic Trajectory of Particle with Bx = 1.0, Ey = 0.8,'
          'Ex = 0.01')
plt.savefig('rel_traj_Bx_1_Ey_0.8_Ex_0.01.png', dip=300)
plt.show()

# plot gamma evolution
fig = plt.figure()
gtlist = []
time = range(len(yarr))

# calculate gamma for plot
for i in time:
    gammat = np.sqrt(yarr[i, 3]**2 + yarr[i, 4]**2 + yarr[i, 5]**2)
    gtlist.append(gammat)

plt.plot(time, gtlist)
plt.title('Gamma Evolution With Time for Bx = 1.0, Ey = 0.8, Ex = 0.01')
plt.xlabel('Time (seconds)')
plt.ylabel('Gamma Value')
plt.savefig('gamma_t_Bx_1_Ey_0.8_Ex_0.01.png', dip=300)
plt.show()

# %% Create plot showing gamma evolution and trajectory over time with
# Ey=0.8, Ex=1

# initialize electric and magnetic field strengths
Bx = 1
By = 0
Bz = 0
Ex = 1
Ey = 0.8*Bx
Ez = 0
qmc = q/(m*c)

# integrate ODE to get fucntion
yarr = odeint(deri, init, t, args=(qmc, c, Ex, Ey, Ez, Bx, By, Bz))

# plot relativistic trajectory
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(yarr[:, 0], yarr[:, 1], yarr[:, 2], marker='.')
plt.title('Relativistic Trajectory of Particle with Bx = 1.0, Ey = 0.8,'
          'Ex = 1')
plt.savefig('rel_traj_Bx_1_By_0.8_Bx_1.png', dip=300)
plt.show()

# plot gamma evolution
fig = plt.figure()
gtlist = []
time = range(len(yarr))

# calculate gamma for plot
for i in time:
    gammat = np.sqrt(yarr[i, 3]**2 + yarr[i, 4]**2 + yarr[i, 5]**2)
    gtlist.append(gammat)

plt.plot(time, gtlist)
plt.title('Gamma Evolution With Time for Bx = 1.0, Ey = 0.8, Ex = 1')
plt.xlabel('Time (seconds)')
plt.ylabel('Gamma Value')
plt.savefig('gamma_t_Bx_1_Ey_0.8_Ex_1.png', dip=300)
plt.show()
