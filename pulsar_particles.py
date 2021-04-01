# -----------------------------------------------------------------------------
# File: pulsar_particles.py
# Name: Evan Davis
# Description: Creates 2D and 3D visualizations/animations of a pulsar and the
# particles in its magnetosphere and their characteristics
# Usage: External file allprt01_125133.dat is needed to run this program.
# One may use a os.chdir command to change directory to where the .dat file is
# located if the program does not detect the .dat file
# -----------------------------------------------------------------------------

# %% Load libraries

import os
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import interpolate

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

# %% Create distribution of gamma values of particles

gvals = np.zeros(numprt)  # create array for gamma values for all particles
nbins = 500  # number of bins
bins = np.zeros(nbins)  # create array of bins

# calculate log of gamma values for each particle (x, y, z components are
# attributes 3, 4, 5)
for i in np.arange(numprt):
    gvals[i] = np.log10(np.sqrt(1 + (datars[i, 3])**2 + (datars[i, 4])**2 +
                        (datars[i, 5])**2))

# add infinitessimal to max and min of gamma values to avoid encountering 0
# when calculating bin placement
gminimum = gvals.min()-1.4599E-6
gmaximum = gvals.max()+1.4599E-6
gbinlength = (gmaximum-gminimum)/nbins  # calculate bin length

# place gamma values for each particle into bins
for i in range(numprt):
    kbin = np.int((gvals[i] - gminimum) / gbinlength)
    bins[kbin] = bins[kbin] + 1

# calculate middle of the bins for plotting visuals
gmidbins = np.linspace(gminimum + gbinlength / 2, gmaximum - gbinlength / 2,
                       nbins)


# use lists for positive/negative particles since we don't know how large to
# make an array for them
listp = []
listn = []
nbins = 500
binsp = np.zeros(nbins)
binsn = np.zeros(nbins)

# calculate log of gamma values, assign to list for positive/negative
# where attribute 6 represents charge
for i in np.arange(numprt):
    clg = np.log10(np.sqrt(1 + (datars[i, 3])**2 + (datars[i, 4])**2 +
                           (datars[i, 5])**2))
    if datars[i, 6] == 1.0:
        listp.append(clg)
    elif datars[i, 6] == -1.0:
        listn.append(clg)

# convert lists into arrays and find size
listp = np.array(listp)
listn = np.array(listn)
listpnum = np.size(listp)
listnnum = np.size(listn)

# add infintessimals and calculate binlength
minp = listp.min() - 1.4599E-6
maxp = listp.max() + 1.4599E-6
minn = listn.min() - 1.4599E-6
maxn = listn.max() + 1.4599E-6
pbinlength = (maxp-minp) / nbins
nbinlength = (maxn-minn) / nbins

# place gamma values for positive particles in bins
for i in range(listpnum):
    kpbin = np.int((listp[i] - minp) / pbinlength)
    binsp[kpbin] = binsp[kpbin] + 1

# place gamma values for negative particles in bins
for i in range(listnnum):
    knbin = np.int((listn[i] - minn) / nbinlength)
    binsn[knbin] = binsn[knbin] + 1

# calculate middle of the bins for positive/negative for visuals
midbinsp = np.linspace(minp + pbinlength / 2, maxp - pbinlength / 2, nbins)
midbinsn = np.linspace(minn + nbinlength / 2, maxn - nbinlength / 2, nbins)

# plot gamma distribution for negative, positive, and total particles
plt.plot(midbinsp, np.log10(binsp), label='Positron Distribution')
plt.plot(midbinsn, np.log10(binsn), label='Electron Distribution')
plt.plot(gmidbins, np.log10(bins), label='Total Distribution')
plt.legend(loc='lower left')
plt.xlabel('Gamma Value')
plt.ylabel('Frequency of Gamma Values')
plt.title('Distribution of Gamma Values')
plt.savefig('pos_neg_gammadist.png', dip=300)

# %% Create 3D plots/animation showing pulsar and the surrounding particles

# initialize lists to place x, y, z coords for positive/negative particles
xplist = []
yplist = []
zplist = []
xnlist = []
ynlist = []
znlist = []

# calculate gamma values and distance from pulsar center for each particle
for i in np.arange(numprt):
    r = np.sqrt(datars[i, 0]**2 + datars[i, 1]**2 + datars[i, 2]**2)
    g = np.sqrt(1 + (datars[i, 3])**2 + (datars[i, 4])**2 + (datars[i, 5])**2)

# only add high energy (high gamma) particles that are far away enough
    if r <= 2 and g > 10**2.5:
        if datars[i, 6] == 1.0:
            xplist.append(datars[i, 0])
            yplist.append(datars[i, 1])
            zplist.append(datars[i, 2])
        elif datars[i, 6] == -1.0:
            xnlist.append(datars[i, 0])
            ynlist.append(datars[i, 1])
            znlist.append(datars[i, 2])

# plot 3D figure of high energy positive particle positions
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(xplist, yplist, zplist, c='r', marker='.')
plt.title('Positions of Positrons')
plt.xlabel('Radial X Distance')
plt.ylabel('Radial Y Distance')
ax.set_zlabel('Radial Z Distance')
plt.savefig('pos_highenergy_position.png', dip=300)

# plot 3D figure of high energy negative particle positions
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(xnlist, ynlist, znlist, c='b', marker='.')
plt.title('Positions of Electrons')
plt.xlabel('Radial X Distance')
plt.ylabel('Radial Y Distance')
ax.set_zlabel('Radial Z Distance')
plt.savefig('neg_highenergy_position.png', dip=300)

# plot 3D figure of all high energy particle positions
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(xplist, yplist, zplist, c='r', marker='.', label='positrons')
ax.scatter(xnlist, ynlist, znlist, c='b', marker='.', label='electrons')
plt.title('Positions of Positrons and Electrons')
plt.xlabel('Radial X Distance')
plt.ylabel('Radial Y Distance')
plt.legend()
ax.set_zlabel('Radial Z Distance')
plt.savefig('pos_neg_highenergy_position.png', dip=300)

# create solid sphere to place into 3D plots to represent the pulsar
theta = np.linspace(0, 2*np.pi, 200)
phi = np.linspace(0, np.pi, 200)
xsphere = 0.25*np.outer(np.cos(theta), np.sin(phi))
ysphere = 0.25*np.outer(np.sin(theta), np.sin(phi))
zsphere = 0.25*np.outer(np.ones(np.size(theta)), np.cos(phi))

# plot sphere in the middle of the middle of the last plot
ax.plot_surface(xsphere, ysphere, zsphere, color='k')
plt.savefig('pos_neg_highenergy_position_sphere.png', dip=300)


def init():
    ax.plot_surface(xsphere, ysphere, zsphere, color='k')
    return fig


def animate(i):
    ax.view_init(elev=0, azim=i)
    return fig


# create animation showing a rotation around the 3D plot
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=300,
                               interval=50, blit=False)
anim.save('highenergy_particle_animation.mp4', fps=30,
          extra_args=['-vcodec', 'libx264'])

# %% Create phase plot of pulsar

# define limits of phi and xi, determine bin length for both
p2 = 2*np.pi
phibins = 100
xibins = 100
phimax = p2
ximax = np.pi
phibinlength = (phimax-0)/phibins
xibinlength = (ximax-0)/xibins
skymap = np.zeros((xibins, phibins))

# calculation of modulus of velocity, distance, simplified gamma
for i in np.arange(numprt):
    v = [datars[i, 3], datars[i, 4], datars[i, 5]]
    modv = np.sqrt(datars[i, 3]**2 + datars[i, 4]**2 + datars[i, 5]**2)
    r = [datars[i, 0], datars[i, 1], datars[i, 2]]
    modr = np.sqrt(datars[i, 0]**2 + datars[i, 1]**2 + datars[i, 2]**2)
    gmsim = np.sqrt(1 + (datars[i, 3])**2 + (datars[i, 4])**2 +
                    (datars[i, 5])**2)
    gamma = datars[i, 31]
    Rc = datars[i, 8]

# calculate energy of particles not too close nor too far with high energy
# and place them into a phase plot (xi is vertical angle, phi is horizontal)
    if modr >= 0.4 and modr <= 2.0 and gmsim > 30 and Rc > 0:
        emsc = gamma**4 / Rc**2
        xi = np.arccos(datars[i, 5] / modv)
        delphi1 = np.arctan2(datars[i, 4], datars[i, 3])
        delphi2 = np.dot(r, v) / modv
        phi = (-1 * delphi1 - delphi2) % p2
        kbin = np.int((phi - 0) / phibinlength)
        jbin = np.int((xi - 0) / xibinlength)
        skymap[jbin][kbin] = skymap[jbin][kbin] + emsc

# create phase plot showing energy of particles at specific locations
# around the pulsar
fig = plt.figure()
plt.imshow(skymap, cmap='brg', origin='lower', extent=[0, phimax, 0, ximax])
plt.xlabel('Phi (Radians)')
plt.ylabel('Xi (Radians)')
plt.title('Phase Plot')
colorbar = plt.colorbar()
colorbar.ax.set_ylabel('E_msc')
plt.savefig('phi_v_xi.png', dip=300)

# %% Create animation of light curve based on xi

# interpolate equation for calculating intensity
phis = np.linspace(0 + phibinlength / 2, 2 * np.pi - phibinlength / 2, 100)
xis = np.linspace(0 + xibinlength / 2, np.pi - xibinlength / 2, 100)
fint = interpolate.interp2d(phis, xis, skymap, kind='linear')

# create range of xi values we want light curve for
fig = plt.figure(figsize=(8, 8))
frmax = 400
ax = plt.axes(xlim=(0, p2), ylim=(0, 1.1))
xirange = np.linspace(p2 / 12, p2 / 2, frmax)


# create function to plot light curve for a specific xi value
def lightcurves(i, ximax, frmax):
    ax.clear()
    xi = xirange[i]
    ylc = []

    for y in phis:
        ylc.append(fint(y, xi))

    ylc = ylc / np.max(ylc)
    plot = plt.plot(phis, ylc)
    plt.title('Light Curve for Xi = %.2f$^\\circ$' %
              (xi * 180 / np.pi))
    plt.xlabel('Phi (Radians)')
    plt.ylabel('Normalized Intensity')
    return plot


# create animation of light curve based on function and range of xi values
anim = animation.FuncAnimation(fig, lightcurves, fargs=(p2 / 2, frmax),
                               frames=frmax, interval=5, blit=False)
anim.save('light_curve_evolution_with_xi.mp4', fps=30,
          extra_args=['-vcodec', 'libx264'])
plt.show()
