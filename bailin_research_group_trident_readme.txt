This file assumes you have passing familiarity with the content of the base Trident readme and it's contents-- specifically about installing Trident and using basic functions. DO NOT follow the base Trident installation guide, but read it so this one will be easier to understand!

1) Assuming yt is installed correctly, adding Trident should be relatively easy:

$git clone https://github.com/jmorgan-15/trident.git
$cd trident
$git checkout -b bailin_research_group_trident
$python -m pip install --upgrade pip
$python -m pip install --user -e .

2) Install ionization table. The first time you start Trident it will need an internet connection to install an ionization table. You can start this by entering python and:

  >>>import trident

This should prompt the download. If not, try running the test script contained in this folder BUT MAKE SURE TO MOVE IT OUTSIDE OF YT/TRIDENT FOLDERS FIRST

The main addition is the make_my_ray function. As a reminder, this function takes quite a few arguments:

my_ray=trident.make_my_ray(ds, start_position, end_position, snr, data_filename, spectral_filename, complete_filename, line_database, interactive)

start_position, end_position: hand the function any 2 points in 3D space here. If unitless, the function will assume you mean simulation units. If you pass it units, it will handle them correctly.

snr: signal to noise ratio, we default to 18 at the behest of Varsha/Jianghao

data_filename: the name of the initial ray file yt makes. We make our own much more complex ray files, we USED to just let these files be overwritten. However, there are intractable differences between ArepoHDF5Dataset files and YTLightRayDataset files that make storing all the data in one file problematic. Or at least, it's a problem I don't have time to fix. Now we should SAVE these files as well, as they contain spectral data that lets us post-process things and get particle-level info about spectral properties.

spectral_filename: the name for the plots of spectra and .txt files yt makes. This is only used if interactive=True. If interactive=False, the only data product will be the complete hdf5 files. Which do contain all the data in all other files.

complete_filename: this is the file we actually care about! The hdf5 file with all gas cell fields and spectral properties included

line_database: the name of our line list file; defaults to MortCashTri.txt. If you want to use a different list it needs to be in the trident/data/line_lists folder.

ineractive: if True, make_my_ray() will return the loaded dataset it created (ie, it does yt.load(complete_filename)) and also produce the spectral plots/.txt files. Otherwise it generates only the main hdf5 file and returns the filename (so if you wanted you could load it up after yourself)

store_observables: if True, the created datafile will also store particle-level spectral info for all lines and instruments by default, or only those specified by the user. 

There are other arguments, but most are either not used or have defaults that I doubt will change. If you really want to see them, make_my_ray() is in trident/ray_generator.py.





The 'store_observables=True' keyword can be passed to trident.make_my_ray() to store things like the optical depth for each line for each gas particle. This is best done for only a handful of lines, specified by using the 'lines=' keyword. Set to a list of either specific ions/wavelengths ('H I 1216') or just an ion/just an element type ('C IV' or 'C'). Each instrument stores data about the lines within its wavelength range. Recall that both line subsets and instrument subsets can be passed as arguments to the make_my_ray() function, but it does all lines and instruments by default. When accessing the data, you can see which lines were recorded like so:

  >>>my_ray.rayvals.keys()  #to give you a list of instruments used
  >>>my_ray.rayvals['COS-G130M'].keys()  #to give a list of recorded lines within a given instruments range (chose G130M here)
  >>>my_ray.rayvals['COS-G130M']['O VI 1032'].keys() #returns the list of by-particle spectral properties
  
If you already know which instrument/line ect you want, you can just use the los() function:
  >>>my_ray.los('COS-G130M', line='O VI 1032', line_spectral_prop='column_density')

Keyword arguments are necessary though. If you made a ray without storing observables, it can be found by post-processing the original ray.h5 file. These files are much smaller and only contain data about the chemical makeup, temp, density, and LOS velocity. Load it up, define an instrument, and perform the observation. This allows for post-processing one instrument and line (or a subset of instruments and lines) at a time, after you've decided which you're interested in from the full spectrum. 

  >>>my_default_x_ray=yt.load('{ray_filename}_ray.h5')
  >>>my_instrument=trident.ray_generator.instruments['COS-G130M']    #I have instrument data stored in this module, properties in order: [lambda_min, lambda_max, dlambda, LSF_filename]
  >>>sg=trident.SpectrumGenerator(lambda_min=my_instrument[0], lambda_max=my_instrument[1], dlambda=my_instrument[2], line_database='MortCashTri.txt')
  >>>sg.make_spectrum(my_default_x_ray, lines=['H', 'C', 'N', 'O'], store_observables=True) 

Sometimes, so certain are images easier to make, we want the coordinate system centered on the galaxy center, but with the same axes as the simulation. This is the one type of coordinate system not already stored. Finally-- I had to think about this a while to make sure it was right-- we can rotate our rays to fit the simulation axes easily. Below I make a ray that goes along the simulation​ y-axis with an impact param of ip on the simulation​ x-axis:
  >>>import yt
  >>>yt.enable_plugins()
  >>>import unyt
  >>>ds=yt.load(halo_file_name)
  >>>theta, phi=ds.hs('baryonic_L_orientation').value  #get rid of "dimensionless" units
  >>>R_z, R_y=np.array([[np.cos(phi), -1*np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]]), np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-1*np.sin(theta), 0, np.cos(theta)]])
  >>>x_raw_start, x_raw_end=unyt.unyt_array([ip, 1000, 0], units='kpc'), unyt.unyt_array([ip, -1000, 0], units='kpc')
  >>>x_start, x_end=x_raw_start @ (R_z @ R_y).T, x_raw_end @ (R_z @ R_y).T  #matrix multiplication

And now we can feed x_start and x_end into trident.make_my_ray! 

A final note: I mentioned stars above, and some of the scripts have code for dealing with star particles. Our newest data files contain stellar info, though I haven't used it for anything yet.
