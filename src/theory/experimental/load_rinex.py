import os
import matplotlib.pyplot as plt
import numpy as np

# u-blox example
datadir = '../data/wine/static'
navfile = 'base_COM7___460800_250416_143415.nav'
rovfile = 'rover1_COM13___460800_250416_143419.obs'
basefile = 'base_COM7___460800_250416_143415.obs'
os.chdir(datadir)

# set run parameters
maxepoch = 1005  # max # of epochs, used for debug, None = no limit
basepos = []  # default to not specified here

# import rtklib files
import src.__test_config as cfg
import src.rinex as rn
from src.rtkpos import rtkinit

nav = rtkinit(cfg)

# load rover obs
rov = rn.rnx_decode(cfg)
print('Reading rover obs...')
rov.decode_obsfile(nav, rovfile, maxepoch)

# load base obs
base = rn.rnx_decode(cfg)
print('Reading base obs...')
base.decode_obsfile(nav, basefile, None)
if basepos != []:
    nav.rb = basepos
elif nav.rb[0] == 0:
    nav.rb = base.pos

# load nav data from rover obs
print('Reading nav data...')
rov.decode_nav(navfile, nav)

# OBS data
# .P : pseudorange



#
print('Plotting...')
num_epochs = len(rov.obslist)
sats = []
for i in range(num_epochs):
    sats.append(rov.obslist[i].sat)
unique_sats = np.unique(np.concatenate(sats).ravel())

plt.figure()
plt.title('Satellites in view over time')
for i in range(num_epochs):
    for j in range(len(sats[i])):
        plt.plot(i, sats[i][j], color='r', marker='o', ls='none')
plt.xlabel('Epoch')
plt.ylabel('Sat ID')
plt.show()
plt.close()