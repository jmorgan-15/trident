import h5py
import six
import requests
import pickle
import astropy
import pickle
import numpy as np
import illustris_python as il
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import LogNorm
from matplotlib import cm
from matplotlib import colorbar
from mpl_toolkits.mplot3d import Axes3D
import yt
import trident



#TO BE SAFE, MOVE THIS FILE OUTSIDE OF THE YT/TRIDENT INSTALLATION FOLDERS
#yt.enable_plugins(/path/to/plugins/file.py)  #use this line if you could not store the plugins file in the default folder (~/.config/yt/) or if it failed to load from that folder.
yt.enable_plugins()
import unyt
import random


#AFTER FULLY INSTALLING ALL CUSTOM FILES, SHOULD BE ABLE TO RUN THIS SCRIPT. 
halos=[1040]
plotpath='.'
  
for haloid in halos:  
  #filename=f'./TNG100-1/cutouts/84_layered_new/cutout_{haloid}.hdf5'
  filename=f'./cutout_{haloid}.hdf5'
  print('print statement 1 in raw_spectra_generator_test')
  ds=yt.load(filename)
  r_vir=ds.hs('Group_R_Crit200_kpc')
  #ip_array=[20, 150]
  ip_array=[10]
  for ip in ip_array:   
    print(haloid, ip) 
      
    #Down x-axis, along y
    x_ray=[unyt.unyt_array([ip, -3*r_vir.value, 0], units='kpc'), unyt.unyt_array([ip, 3*r_vir.value, 0], units='kpc')]  #we are constructing these points like the center of the halo is at [0, 0, 0]; adding the actual origin is handled within make_my_ray
    print('first ray at IP')
    ray_filename=f'{haloid}_x_alongy_'+str(format(ip, '.3f'))
    my_x_ray_stored=trident.make_my_ray(ds, start_position=x_ray[0], end_position=x_ray[1], snr=18, data_filename=f'{plotpath}/{ray_filename}_ray.h5', complete_filename=f'{plotpath}/{ray_filename}_stored.hdf5', spectral_filename=f'{plotpath}/{ray_filename}', interactive=True, halo=haloid, apply_lsf=True, store_observables=False, lines=['H', 'C', 'N', 'O', 'Mg'])
    




my_default_x_ray=yt.load(f'{ray_filename}_ray.h5')
my_instrument=trident.ray_generator.instruments['COS-G130M']
sg=trident.SpectrumGenerator(lambda_min=my_instrument[0], lambda_max=my_instrument[1], dlambda=my_instrument[2], line_database='MortCashTri.txt')
sg.make_spectrum(my_default_x_ray, lines=['H', 'C', 'N', 'O'], store_observables=True)  














