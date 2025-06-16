import os

import georinex as gr
import matplotlib.pyplot as plt
import numpy as np
import pymap3d as pm


# u-blox example
datadir = '/home/guillem/Code/rtklib-py/data/wine/static'
navfile = 'base_COM7___460800_250416_143415.nav'
rovfile = 'rover1_COM13___460800_250416_143419.obs'
basefile = 'base_COM7___460800_250416_143415.obs'

# load rover obs
print('Reading rover obs...')
rov = gr.load(os.path.join(datadir, rovfile), tlim=['2025-04-16T14:34', '2025-04-16T14:36'])

df = rov.to_dataframe()
df = df.reset_index()

print(f'Position ECEF: {rov.position}')
llh = pm.ecef2geodetic(*rov.position)
print(f'Position LLH: {llh}')

# load base obs


#
print('Plotting...')

plt.figure()
ax = plt.gca()
plt.title('Satellites in view over time')
for sv in df.sv.unique():
    df_sv = df[df['sv'] == sv].sort_values(by='time')
    df.set_index('time')['C1C'].plot(ax=ax, label=sv, alpha=0.5)
plt.xlabel('Time')
plt.ylabel('C1C')

plt.legend()
plt.show()
plt.close()
