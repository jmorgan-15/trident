"""
SpectrumGenerator class and member functions.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Trident Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#-----------------------------------------------------------------------------

from trident.light_ray import \
    LightRay
from yt.loaders import \
    load
from trident.config import \
    ion_table_filepath
from trident.line_database import \
    LineDatabase, \
    uniquify
from trident.roman import \
    from_roman
from yt.data_objects.static_output import \
    Dataset
from trident.ion_balance import \
    atomic_number

#########ADDITIONAL IMPORTS######################
#this is actually the same thing as unyt_array, but unyt used to be a part of yt before becoming its own thing. Trident still uses it like its a party of yt, but if you look into the yt code, it just says "YTArray=unyt_array". Since we are messing with Trident code, "when in Rome"....so we use YTArray
from yt.units import \
    YTArray 
from trident.spectrum_generator import \
    SpectrumGenerator 

import numpy as np
import copy
import h5py
import shutil
def make_simple_ray(dataset_file, start_position, end_position,
                    lines=None, ftype="gas", fields=None,
                    solution_filename=None, data_filename=None,
                    trajectory=None, redshift=None, field_parameters=None,
                    setup_function=None, load_kwargs=None,
                    line_database=None, ionization_table=None):
    """
    Create a yt LightRay object for a single dataset (eg CGM).  This is a
    wrapper function around yt's LightRay interface to reduce some of the
    complexity there.

    A simple ray is a straight line passing through a single dataset
    where each gas cell intersected by the line is sampled for the desired
    fields and stored.  Several additional fields are created and stored
    including ``dl`` to represent the path length in space
    for each element in the ray, ``v_los`` to represent the line of
    sight velocity along the ray, and ``redshift``, ``redshift_dopp``, and
    ``redshift_eff`` to represent the cosmological redshift, doppler redshift
    and effective redshift (combined doppler and cosmological) for each
    element of the ray.

    A simple ray is typically specified by its start and end positions in the
    dataset volume.  Because a simple ray only probes a single output, it
    lacks foreground absorbers between the observer at z=0 and the redshift
    of the dataset that one would naturally encounter.  Thus it is usually
    only appropriate for studying the circumgalactic medium rather than
    the intergalactic medium.

    This function can accept a yt dataset already loaded in memory,
    or it can load a dataset if you pass it the dataset's filename and
    optionally any load_kwargs or setup_function necessary to load/process it
    properly before generating the LightRay object.

    The :lines: keyword can be set to automatically add all fields to the
    resulting ray necessary for later use with the SpectrumGenerator class.
    If the necessary fields do not exist for your line of choice, they will
    be added to your dataset before adding them to the ray.

    **Parameters**

    :dataset_file: string or yt Dataset object
    
        Either a yt dataset or the filename of a dataset on disk.  If you are
        passing it a filename, consider usage of the ``load_kwargs`` and
        ``setup_function`` kwargs.

    :start_position, end_position: list of floats or YTArray object

        The coordinates of the starting and ending position of the desired
        ray.  If providing a raw list, coordinates are assumed to be in
        code length units, but if providing a YTArray, any units can be
        specified.

    :lines: list of strings, optional

        List of strings that determine which fields will be added to the ray
        to support line deposition to an absorption line spectrum.  List can
        include things like "C", "O VI", or "Mg II ####", where #### would be
        the integer wavelength value of the desired line.  If set to 'all',
        includes all possible ions from H to Zn. :lines: can be used
        in conjunction with :fields: as they will not override each other.
        Default: None

    :ftype: string, optional

        This is now deprecated and unnecessary.
        Default: "gas"

    :fields: list of strings, optional

        The list of which fields to store in the output LightRay.
        See :lines: keyword for additional functionality that will add fields
        necessary for creating absorption line spectra for certain line
        features.
        Default: None

    :solution_filename: string, optional

        Output filename of text file containing trajectory of LightRay
        through the dataset.
        Default: None

    :data_filename: string, optional
    
        Output filename for ray data stored as an HDF5 file.  Note that
        at present, you *must* save a ray to disk in order for it to be
        returned by this function.  If set to None, defaults to 'ray.h5'.
        Default: None

    :trajectory: list of floats, optional

        The (r, theta, phi) direction of the LightRay.  Use either end_position
        or trajectory, but not both.
        Default: None

    :redshift: float, optional

        Sets the highest cosmological redshift of the ray.  By default, it will
        use the cosmological redshift of the dataset, if set, and if not set,
        it will use a redshift of 0.
        Default: None

    :field_parameters: optional, dict
        Used to set field parameters in light rays. For example,
        if the 'bulk_velocity' field parameter is set, the relative
        velocities used to calculate peculiar velocity will be adjusted
        accordingly.
        Default: None.

    :setup_function: function, optional

        A function that will be called on the dataset as it is loaded but
        before the LightRay is generated.  Very useful for adding derived
        fields and other manipulations of the dataset prior to LightRay
        creation.
        Default: None

    :load_kwargs: dict, optional

        Dictionary of kwargs to be passed to the yt "load" function prior to
        creating the LightRay.  Very useful for many frontends like Gadget,
        Tipsy, etc. for passing in "bounding_box", "unit_base", etc.
        Default: None

    :line_database: string, optional

        For use with the :lines: keyword. If you want to limit the available
        ion fields to be added to those available in a particular subset,
        you can use a :class:`~trident.LineDatabase`.  This means when you
        set :lines:='all', it will only use those ions present in the
        corresponding LineDatabase.  If :LineDatabase: is set to None,
        and :lines:='all', it will add every ion of every element up to Zinc.
        Default: None

    :ionization_table: string, optional

        For use with the :lines: keyword.  Path to an appropriately formatted
        HDF5 table that can be used to compute the ion fraction as a function
        of density, temperature, metallicity, and redshift.  When set to None,
        it uses the table specified in ~/.trident/config
        Default: None

    **Example**

    Generate a simple ray passing from the lower left corner to the upper
    right corner through some Gizmo dataset:

    >>> import trident
    >>> import yt
    >>> ds = yt.load('path/to/dataset')
    >>> ray = trident.make_simple_ray(ds,
    ... start_position=ds.domain_left_edge, end_position=ds.domain_right_edge,
    ... lines=['H', 'O', 'Mg II'])
    """
    if load_kwargs is None:
        load_kwargs = {}
    if fields is None:
        fields = []
    if data_filename is None:
        data_filename = 'ray.h5'

    if isinstance(dataset_file, str):
        ds = load(dataset_file, **load_kwargs)
    elif isinstance(dataset_file, Dataset):
        ds = dataset_file

    lr = LightRay(ds, load_kwargs=load_kwargs)

    if ionization_table is None:
        ionization_table = ion_table_filepath

    # Include some default fields in the ray to assure it's processed correctly.

    fields = _add_default_fields(ds, fields)

    # If 'lines' kwarg is set, we need to get all the fields required to
    # create the desired absorption lines in the grid format, since grid-based
    # fields are what are directly probed by the LightRay object.

    # We first determine what fields are necessary for the desired lines, and
    # inspect the dataset to see if they already exist.  If so, we add them
    # to the field list for the ray.  If not, we have to create them.

    if lines is not None:

        ion_list = _determine_ions_from_lines(line_database, lines)
        fields = _determine_fields_from_ions(ds, ion_list, fields)

    # To assure there are no fields that are double specified or that collide
    # based on being specified as "density" as well as ("gas", "density"),
    # we will just assume that all non-tuple fields requested are ftype "gas".
    for i in range(len(fields)):
        if isinstance(fields[i], str):
            fields[i] = ('gas', fields[i])
    fields = uniquify(fields)

    return lr.make_light_ray(start_position=start_position,
                             end_position=end_position,
                             trajectory=trajectory,
                             fields=fields,
                             setup_function=setup_function,
                             solution_filename=solution_filename,
                             data_filename=data_filename,
                             field_parameters=field_parameters,
                             redshift=redshift)

def make_compound_ray(parameter_filename, simulation_type,
                      near_redshift, far_redshift,
                      lines=None, ftype='gas', fields=None,
                      solution_filename=None, data_filename=None,
                      use_minimum_datasets=True, max_box_fraction=1.0,
                      deltaz_min=0.0, minimum_coherent_box_fraction=0.0,
                      find_outputs=False, seed=None,
                      setup_function=None, load_kwargs=None,
                      line_database=None, ionization_table=None,
                      field_parameters = None):
    """
    Create a yt LightRay object for multiple consecutive datasets (eg IGM).
    This is a wrapper function around yt's LightRay interface to reduce some
    of the complexity there.

    .. note::

        The compound ray functionality has only been implemented for the
        Enzo and Gadget/Gizmo codes.  If you would like to help us implement
        this functionality for your simulation code, please contact us
        about this on the mailing list.

    A compound ray is a series of straight lines passing through multiple
    consecutive outputs from a single cosmological simulation to approximate
    a continuous line of sight to high redshift.

    Because a single continuous ray traversing a simulated volume can only
    cover a small range in redshift space (e.g. 100 Mpc only covers the
    redshift range from z=0 to z=0.023), the compound ray passes rays through
    multiple consecutive outputs from the same simulation to approximate the
    path of a single line of sight to high redshift.  By probing all of the
    foreground material out to any given redshift, the compound ray is
    appropriate for studies of the intergalactic medium and circumgalactic
    medium.

    By default, it selects a random starting location and trajectory in
    each dataset it traverses, to assure that the same cosmological structures
    are not being probed multiple times from the same direction.  In doing
    this, the ray becomes discontinuous across each dataset.

    The compound ray requires the parameter_filename of the simulation run.
    This is *not* the dataset filename from a single output, but the parameter
    file that was used to run the simulation itself.  It is in this parameter
    file that the output frequency, simulation volume, and cosmological
    parameters are described to assure full redshift coverage can be achieved
    for a compound ray.  It also requires the simulation_type of the simulation.

    Unlike the simple ray, which is specified by its start and end positions
    in the dataset volume, the compound ray requires the near_redshift and
    far_redshift to determine which datasets to use to get full coverage
    in redshift space as the ray propagates from near_redshift to far_redshift.
    
    Like the simple ray produced by :class:`~trident.make_simple_ray`,
    each gas cell intersected by the LightRay is sampled for the desired
    fields and stored.  Several additional fields are created and stored
    including ``dl`` to represent the path length in space
    for each element in the ray, ``v_los`` to represent the line of
    sight velocity along the ray, and ``redshift``, ``redshift_dopp``, and
    ``redshift_eff`` to represent the cosmological redshift, doppler redshift
    and effective redshift (combined doppler and cosmological) for each
    element of the ray.

    The :lines: keyword can be set to automatically add all fields to the
    resulting ray necessary for later use with the SpectrumGenerator class.

    **Parameters**

    :parameter_filename: string

        The simulation parameter file *not* the dataset filename

    :simulation_type: string

        The simulation type of the parameter file.  At present, this
        functionality only works with "Enzo" and "Gadget" yt frontends.

    :near_redshift, far_redshift: floats

        The near and far redshift bounds of the LightRay through the
        simulation datasets.

    :lines: list of strings, optional

        List of strings that determine which fields will be added to the ray
        to support line deposition to an absorption line spectrum.  List can
        include things like "C", "O VI", or "Mg II ####", where #### would be
        the integer wavelength value of the desired line.  If set to 'all',
        includes all possible ions from H to Zn. :lines: can be used
        in conjunction with :fields: as they will not override each other.
        Default: None

    :ftype: string, optional

        This is now deprecated and unnecessary.
        Default: "gas"

    :fields: list of strings, optional

        The list of which fields to store in the output LightRay.
        See :lines: keyword for additional functionality that will add fields
        necessary for creating absorption line spectra for certain line
        features.
        Default: None

    :solution_filename: string, optional

        Output filename of text file containing trajectory of LightRay
        through the dataset.
        Default: None

    :data_filename: string, optional

        Output filename for ray data stored as an HDF5 file.  Note that
        at present, you *must* save a ray to disk in order for it to be
        returned by this function.  If set to None, defaults to 'ray.h5'.
        Default: None

    :use_minimum_datasets: bool, optional

        Use the minimum number of datasets to make the ray continuous
        through the supplied datasets from the near_redshift to the
        far_redshift.  If false, the LightRay solution will contain as many
        datasets as possible to enable the light ray to traverse the
        desired redshift interval.
        Default: True

    :max_box_fraction: float, optional

        The maximum length a light ray segment can be in order to span the
        redshift interval from one dataset to another in units of the domain
        size.  Values larger than 1.0 will result in LightRays crossing the
        domain of a given dataset more than once, which is generally undesired.
        Zoom-in simulations can use a value equal to the length of the
        high-resolution region so as to limit ray segments to that size.  If
        the high-resolution region is not cubical, the smallest size should b
        used.
        Default: 1.0 (the size of the box)

    :deltaz_min: float, optional

        The minimum delta-redshift value between consecutive datasets used
        in the LightRay solution.
        Default: 0.0

    :minimum_coherent_box_fraction: float, optional

        When use_minimum_datasets is set to False, this parameter specifies
        the fraction of the total box width to be traversed before
        rerandomizing the ray location and trajectory.
        Default: 0.0

    :find_outputs: optional, bool

        Whether or not to search for datasets in the current
        directory. This is useful if the number of existing datasets is
        different than what would be predicted by the simulation parameter file.
        Default: False.

    :seed: int, optional

        Sets the seed for the random number generator used to determine the
        location and trajectory of the LightRay as it traverses the
        simulation datasets.  For consistent results between LightRays,
        use the same seed value.
        Default: None

    :setup_function: function, optional

        A function that will be called on the dataset as it is loaded but
        before the LightRay is generated.  Very useful for adding derived
        fields and other manipulations of the dataset prior to LightRay
        creation.
        Default: None

    :load_kwargs: dict, optional

        Dictionary of kwargs to be passed to the yt "load" function prior to
        creating the LightRay.  Very useful for many frontends like Gadget,
        Tipsy, etc. for passing in "bounding_box", "unit_base", etc.
        Default: None

    :line_database: string, optional

        For use with the :lines: keyword. If you want to limit the available
        ion fields to be added to those available in a particular subset,
        you can use a :class:`~trident.LineDatabase`.  This means when you
        set :lines:='all', it will only use those ions present in the
        corresponding LineDatabase.  If :LineDatabase: is set to None,
        and :lines:='all', it will add every ion of every element up to Zinc.
        Default: None

    :ionization_table: string, optional

        For use with the :lines: keyword.  Path to an appropriately formatted
        HDF5 table that can be used to compute the ion fraction as a function
        of density, temperature, metallicity, and redshift.  When set to None,
        it uses the table specified in ~/.trident/config
        Default: None

    :field_parameters: optional, dict
        Used to set field parameters in light rays. For example,
        if the 'bulk_velocity' field parameter is set, the relative
        velocities used to calculate peculiar velocity will be adjusted
        accordingly.
        Default: None.

    **Example**

    Generate a compound ray passing from the redshift 0 to redshift 0.05
    through a multi-output enzo simulation.

    >>> import trident
    >>> fn = 'path/to/simulation/parameter/file'
    >>> ray = trident.make_compound_ray(fn, simulation_type='Enzo',
    ... near_redshift=0.0, far_redshift=0.05, lines=['H', 'O', 'Mg II'])

    Generate a compound ray passing from the redshift 0 to redshift 0.05
    through a multi-output gadget simulation.

    >>> import trident
    >>> fn = 'path/to/simulation/parameter/file'
    >>> ray = trident.make_compound_ray(fn, simulation_type='Gadget',
    ... near_redshift=0.0, far_redshift=0.05, lines=['H', 'O', 'Mg II'])
    """
    if load_kwargs is None:
        load_kwargs = {}
    if fields is None:
        fields = []
    if data_filename is None:
        data_filename = 'ray.h5'

    lr = LightRay(parameter_filename,
                  simulation_type=simulation_type,
                  near_redshift=near_redshift,
                  far_redshift=far_redshift,
                  find_outputs=find_outputs,
                  use_minimum_datasets=use_minimum_datasets,
                  max_box_fraction=max_box_fraction,
                  deltaz_min=deltaz_min,
                  minimum_coherent_box_fraction=minimum_coherent_box_fraction,
                  load_kwargs=load_kwargs)

    if ionization_table is None:
        ionization_table = ion_table_filepath

    # We use the final dataset from the light ray solution in order to test it for
    # what fields are present, etc.  This all assumes that the fields present
    # in this output will be present in ALL outputs.  Hopefully this is true,
    # because testing each dataset is going to be slow and a pain.

    ds = load(lr.light_ray_solution[-1]['filename'])

    # Include some default fields in the ray to assure it's processed correctly.

    fields = _add_default_fields(ds, fields)

    # If 'lines' kwarg is set, we need to get all the fields required to
    # create the desired absorption lines in the grid format, since grid-based
    # fields are what are directly probed by the LightRay object.

    # We first determine what fields are necessary for the desired lines, and
    # inspect the dataset to see if they already exist.  If so, we add them
    # to the field list for the ray or add the necessary fields that can
    # generate them on the ray.

    if lines is not None:

        ion_list = _determine_ions_from_lines(line_database, lines)
        fields = _determine_fields_from_ions(ds, ion_list, fields)

    # To assure there are no fields that are double specified or that collide
    # based on being specified as "density" as well as ("gas", "density"),
    # we will just assume that all non-tuple fields requested are ftype "gas".
    for i in range(len(fields)):
        if isinstance(fields[i], str):
            fields[i] = ('gas', fields[i])
    fields = uniquify(fields)

    return lr.make_light_ray(seed=seed,
                             fields=fields,
                             setup_function=setup_function,
                             solution_filename=solution_filename,
                             data_filename=data_filename,
                             redshift=None, njobs=-1,
                             field_parameters = field_parameters)

def _determine_ions_from_lines(line_database, lines):
    """
    Figure out what ions are necessary to produce the desired lines
    """
    if line_database is not None:
        line_database = LineDatabase(line_database)
        ion_list = line_database.parse_subset_to_ions(lines)
    else:
        ion_list = []
        if lines == 'all' or lines == ['all']:
            for k,v in atomic_number.items():
                for j in range(v+1):
                    ion_list.append((k, j+1))
        else:
            for line in lines:
                linen = line.split()
                if len(linen) >= 2:
                    ion_list.append((linen[0], from_roman(linen[1])))
                elif len(linen) == 1:
                    num_states = atomic_number[linen[0]]
                    for j in range(num_states+1):
                        ion_list.append((linen[0], j+1))
                else:
                    raise RuntimeError("Cannot add a blank ion.")

    return uniquify(ion_list)

def _determine_fields_from_ions(ds, ion_list, fields):
    """
    Figure out what fields need to be added based on the ions present.

    Check if the number_density fields for these ions exist, and if so, add
    them to field list. If not, leave them off, as they'll be generated
    on the fly by SpectrumGenerator as long as we include the 'density',
    'temperature', and appropriate 'metallicity' fields.
    """
    for ion in ion_list:
        atom = ion[0].capitalize()
        ion_state = ion[1]
        nuclei_field = "%s_nuclei_mass_density" % atom
        metallicity_field = "%s_metallicity" % atom
        field = "%s_p%d_number_density" % (atom, ion_state-1)

        # Check to see if the ion field exists.  If so, add
        # it to the ray.  If not, then append the density and the appropriate
        # metal field so one can create the ion field on the fly on the
        # ray itself.
        if ("gas", field) not in ds.derived_field_list:
            fields.append(('gas', 'density'))
            if ('gas', metallicity_field) in ds.derived_field_list:
                fields.append(('gas', metallicity_field))
            elif ('gas', nuclei_field) in ds.derived_field_list:
                fields.append(('gas', nuclei_field))
            elif atom != 'H':
                fields.append(('gas', 'metallicity'))
            else:
                # Don't need metallicity field if we're just looking
                # at hydrogen
                pass
        else:
            fields.append(("gas", field))

    return fields

def _add_default_fields(ds, fields):
    """
    Add some default fields to rays to assure they can be processed correctly.
    """
    if ("gas", "temperature") in ds.derived_field_list:
        fields.append(("gas", 'temperature'))

    # H_nuclei_density should be added if possible to assure that the _log_nH
    # field, which is used as "density" in the ion_balance interpolation to
    # produce ion fields, is calculated as accurately as possible.
    if ('gas', 'H_nuclei_density') in ds.derived_field_list:
        fields.append(('gas', 'H_nuclei_density'))

    return fields
    
    
##########################################################################################################################################################
##################################################################ADDITIONS###############################################################################
#ALWAYS SPECIFY INSTRUMENT PROPERTIES IN SCRIPTS AND DO NOT USE THE INSRTUMENT NAMES! THOSE NAMES ARE LOADED FROM ANOTHER LOCATION (SpectrumGenerator) WHERE THEY HAVE A DIFFERENT DEFINITION! Also, consider changing the definition in SpectrumGenerator to match our preferred properties
cos_grisms=['COS-G130M', 'COS-G160M', 'COS-G185M', 'COS-G225M', 'COS-G285M']
grism_info=[[899, 1469, 0.00997, 'avg_COS_G130M.txt'], [1342, 1798, 0.01223, 'avg_COS_G160M.txt'], [1670, 2127, 0.037, None], [2070, 2527, 0.033, None], [2480, 3229, 0.04, None]]
instruments={cos_grisms[i]:grism_info[i] for i in range(len(cos_grisms))}

#cos_grisms=['COS-G130M', 'COS-G160M']
#grism_info=[[1300, 1400, 0.01, 'avg_COS_G130M.txt'], [1405, 1777, 0.012, 'avg_COS_G160M.txt']]
#grism_info=[[1300, 1400, 0.01, 'avg_COS_G130M.txt'], [1342, 1798, 0.01223, 'avg_COS_G160M.txt']]
#instruments={cos_grisms[i]:grism_info[i] for i in range(len(cos_grisms))}
#sgs=[trident.SpectrumGenerator(lambda_min=grism_info[i][0], lambda_max=grism_info[i][1], dlambda=grism_info[i][2], line_database='lines.txt') for i in range(3)]

#These fields are only stored under 'gas', not 'PartType0'
only_gas_fields=(('gas', 'entropy'), )
#These fields are generated during the initial ray creation (NOT spectra creation). They are then used for spectra post-processing. They are searched for under the 'gas' keyword-- tried to put them under PartType0 and it led to cascading errors with units. Instead, just make a seperate 'gas' group in file and store there
necessary_ray_fields={'velocity_los':'code_velocity', 'v_los':'dimensionless', 'redshift':'dimensionless',  'redshift_dopp':'dimensionless', 'redshift_eff':'dimensionless'}
ray_fields_to_skip=('ParticleIDs', 'density', 'mass', 'metallicity', 'temperature', 'velocity_magnitude', 'x', 'y', 'z', 'H_nuclei_density', 'H_p0_number_density')
#these are all properties that are also generated when you load up the hdf5 file as a dataset; ie, these properties are already listed under ds.derived_field_list

#These suffixes are added to unit-converted data when the file loads up (halo and ray files). So we don't have to save any of the unit-ed quantities to the ray file when making it; if our gas prop has one of these as substring, we skip
used_units=('_msun_pc3', '_flux_erg_s_cm2', '_ergcm3_s', '_solar', '_km2_s2', '_g_km', '_msun_yr', '_kpc', '_km_s', '_erg_s', '_gauss')

ds_index_holder, ray_index_holder, common_values_holder=0, 0, 0
groups=['Config', 'Header', 'halo_and_sub_properties', 'ray_properties', 'PartType5', 'PartType0']
ray_prop_names=['ray_uv', 'ip_uv', 'ip', 'ip_d_along_ray', 'ip_coordinate', 'start_position', 'end_position']
def make_my_ray(ds, start_position, end_position, instruments=instruments, snr=18,
             lines='all', ftype="PartType0", fields=None,
             solution_filename=None, data_filename=None, spectral_filename=None, complete_filename=None,
             trajectory=None, redshift=None, field_parameters=None,
             setup_function=None, load_kwargs=None,
             line_database='MortCashTri.txt', ionization_table=None, interactive=False, halo=None, lims=None, add_qso_spectrum=False, add_milky_way_foreground=False, apply_lsf=True, store_observables=False):     
  origin=YTArray(75000/2, units='code_length', registry=ds.unit_registry)
  #data_filename is for the ray.h5 object; spectral filename is for the images/data files; complete_filename is for the final hdf5 file we make here     
  #start and end positions are assumed to be input as either unyt_arrays or bare lists/arrays. If bare, assume in simulation coords. First thing we do is unscale the start and end points from kpc to code lengths

  if hasattr(start_position, 'units'):
    #start_position, end_position=start_position.to('code_length', registry=ds.unit_registry), end_position.to('code_length', registry=ds.unit_registry)  #for some reason this line breaks, but the line below does not
    start_position, end_position=start_position.to(origin.units), end_position.to(origin.units)
  else:
    start_position, end_position=YTArray(start_position, units='code_length', registry=ds.unit_registry), YTArray(end_position, units='code_length', registry=ds.unit_registry)
    #Now we find the unit vector of our ray, the impact parameter, and the unit vector from the origin to the impact point. We do all this before shifting the ray position to fit where yt thinks the actual particles are
  ray_vec=end_position-start_position
  ray_uvec=ray_vec.value/np.linalg.norm(ray_vec.value)  #unit vectors are ironically unitless
  distance_along_ray=np.sum(-start_position*ray_uvec)  #dot product tells us how far along ray the point closest to the origin is
  x=start_position+distance_along_ray*ray_uvec  #point along ray closest to origin; because origin=galactic center, can do next line easy:
  ip=YTArray(np.linalg.norm(x), units='code_length', registry=ds.unit_registry)   #this is ip in code length units, but comes out of norm() unitless
  uv_to_ip=x.value/ip.value
    
  #Now we make the actual ray
  ray=make_simple_ray(ds, start_position=start_position+origin, end_position=end_position+origin, data_filename=data_filename, lines=lines, ftype='PartType0', line_database=line_database, fields=[('PartType0', 'ParticleIDs')])
  dl=copy.deepcopy(ray.r['gas', 'dl'])
  ray_props=[ray_uvec, uv_to_ip, ip, distance_along_ray, x, start_position, end_position]
  #Everything has been calculated; now we just need to organize/store data in hdf5 file    
  #First we copy/store metadata in a new file
  with h5py.File(complete_filename, 'w') as f:
    for grp in groups:
      f.create_group(grp)
    f['Config'].attrs['VORONOI'], f['Config'].attrs['RAY']=1, 1  
    f['Config'].attrs['STORE_OBSERVABLES']=int(store_observables==True)
  #trying to move away from using attrs, but this set-up is essential to reloading the data files, as yt assumes this arrangement. Same with header below
    f['Header'].attrs['LineDatabase']=line_database
    for key, value in ds.headvals.items():
      if key=='NumPart_ThisFile':
        f['Header'].attrs[key]=np.int32([len(ray.r['gas', 'ParticleIDs']), 0., 0., 0., 0., len(ds.r['PartType5', 'ParticleIDs'])])
      else:
        f['Header'].attrs[key]=value
    for key, value in ds.hsvals.items():
      f['halo_and_sub_properties'].create_dataset(key, data=value)
  
  #Now we store data about the ray itself 
    for i in range(len(ray_props)):
      f['ray_properties'].create_dataset(ray_prop_names[i], data=np.array(ray_props[i]))
      
    #We use our instruments to create spectral data from the ray
    for inst_name, inst_props in instruments.items():   #generate/save spectra for different instruments
      sg=SpectrumGenerator(lambda_min=inst_props[0], lambda_max=inst_props[1], dlambda=inst_props[2], line_database=line_database)
      sg.make_spectrum(ray, lines=lines, store_observables=store_observables)
      if add_qso_spectrum==True:
        sg.add_qso_spectrum()
      if add_milky_way_foreground==True:
        sg.add_milky_way_foreground
      if apply_lsf==True:
        if isinstance(inst_props[3], str):
          sg.apply_lsf(filename=inst_props[3])
        else:
          sg.apply_lsf(function='gaussian', width=4)     #until we get actual LSF files
      if snr!=None:
        sg.add_gaussian_noise(snr)
      f['ray_properties'].create_group(inst_name)
      sg._write_spectrum_hdf5(f['ray_properties'][inst_name], add_to_file=True, filename=complete_filename)
      #if we're storing observables we want to put them in now
      if store_observables==True: 
        for line, properties in sg.line_observables_dict.items():
          f['ray_properties'][inst_name].create_group(line)
          for property_name, value in properties.items():
            if property_name=='EW':
              f['ray_properties'][inst_name][line].create_dataset(property_name, data=np.array(value))
              continue
            try:
              f['ray_properties'][inst_name][line].create_dataset(property_name, data=np.array(value)[ray_index])
            except IndexError as e:
              print(property_name)
              raise e
      #f['ray_properties'].create_dataset(inst_name, data=np.array([sg.lambda_field, sg.tau_field, sg.flux_field, sg.error_field]))
      if interactive==True:
        sg.save_spectrum(spectral_filename+'_'+inst_name+'.txt')
        #sg.plot_spectrum(title='Halo '+str(halo)+' IP='+str(format(ip.to('kpc'), '.3f'))+' '+inst_name, filename=spectral_filename+'_'+inst_name+'.pdf', lambda_limits=lims)
        sg.plot_spectrum(filename=spectral_filename+'_'+inst_name+'.pdf', lambda_limits=lims)

  #Now we store particle data. First, create mask for gas data in ds object
    common_values, ds_index, ray_index=np.intersect1d(ds.gas('ParticleIDs'), ray.r['gas', 'ParticleIDs'], assume_unique=True, return_indices=True)
    l_ray=ray.r['gas', 'l'][ray_index]  #this is a list of l values in the same order as ray properties will be in. We want to make a list of indeces for the ray and for the ds that go in order of increasing 'l'. We can zip these values with ds_index and ray_index to get a new order for both which goes in order of 'l' instead of randomly
    placeholder1=list(zip(l_ray, ds_index, ray_index))
    placeholder1.sort()  #sorts indeces by corresponding 'l' value
    placeholder2=list(zip(*placeholder1))
    common_values, ds_index, ray_index=np.array(placeholder2[0]), np.array(placeholder2[1]), np.array(placeholder2[2])  #These are now the right order of indeces to give particles in order along the ray
    global ds_index_holder, ray_index_holder, common_values_holder
    ds_index_holder, ray_index_holder, common_values_holder=ds_index, ray_index, common_values
    #To keep all data in right order, need to use these "masks" on both ray object and ds object. Even though all ray object particles are used, they aren't in the same order (ie, ray_index is not just [0, 1, 2, 3....]). 
    for field in ds.derived_field_list:
      if 'PartType5' in field:  #if field tuple starts with PartType5
        f['PartType5'].create_dataset(field[1], data=ds.bh(field[1]))
      elif 'PartType0' in field: 
      #elif 'PartType0' in field or 'gas' in field: #tried this line, but many fields are calculated to be the same for both gas, PartType0, so we get copies and it all breaks
        if any(unit in field[1] for unit in used_units):
          continue
        if 'count' in field:
          f['PartType0'].create_dataset(field[1], data=len(ray_index))
        else:
          print(field)
          f['PartType0'].create_dataset(field[1], data=ds.r[field][ds_index])
        #if '_84orientation_78' in field[1]:
          #f['PartType0'].create_dataset(field[1], data=ds.gas(field[1], 78)[ds_index])
        #elif '_78' in field[1]:
          #f['PartType0'].create_dataset(field[1], data=ds.gas(field[1], 78, 78)[ds_index])
        #else:
          #f['PartType0'].create_dataset(field[1], data=ds.gas(field[1])[ds_index])
    #this is annoying and a bit confusing, but in DATASETS 'entropy' is stored under 'gas', in RAYS it is stored under 'PartType0'. But it is the same values either way
 
    for field in only_gas_fields:
      #f['PartType0'].create_dataset(field[1], data=ds.gas(field[1])[ds_index]) 
      f['PartType0'].create_dataset(field[1], data=ds.r[field][ds_index])
      #f['PartType0'].create_dataset(field[1]+'_78', data=ds.gas(field[1], 78)[ds_index])
      

    #We take properties from the dataset first; some auto-generated properties in the ray may have the same name as properties I made up in the dataset. My own properties are better, so we take them. If we find the same property in the ray we just skip it
    for field in ray.derived_field_list:
      if field[1] in ray_fields_to_skip:
        continue
      if field[0]=='gas':
        #if field[1][:9]=='particle_':
        #  print(field)
        if field[1]=='dl':
          data, units=dl[ray_index], 'code_length'
          f['PartType0'].create_dataset(field[1], data=data.to(units, registry=ds.unit_registry))    
        elif field[1]=='l':
          #data, units=l_ray[ray_index], 'code_length'  #BAD this line doesn't work because the order of l_ray is already changed from what's on the disk
          data, units=ray.r['gas', 'l'][ray_index], 'code_length'
          f['PartType0'].create_dataset(field[1], data=data.to(units, registry=ds.unit_registry))    
        elif field[1] in necessary_ray_fields.keys():
          data, units=ray.r[field][ray_index], necessary_ray_fields[field[1]]
          f['PartType0'].create_dataset(field[1], data=data.to(units, registry=ds.unit_registry))  
        #elif 'gas' in field and (field[1]=='l' or 'number_density' in field[1] or 'ion_fraction'in field[1] or 'nuclei_density' in field[1]):
        elif any(substr in field[1] for substr in ('number_density', 'nuclei_density')):
          data, units=ray.r[field][ray_index], '1/cm**3'  
          f['PartType0'].create_dataset(field[1], data=data.to(units, registry=ds.unit_registry))    
        elif 'ion_fraction' in field[1]:
          data, units=ray.r[field][ray_index], 'dimensionless'
          f['PartType0'].create_dataset(field[1], data=data.to(units, registry=ds.unit_registry)) 
          
  #Finally done storing/organizing data; now just save hdf5 file. Note that we are saving several properties from ds that are normally calculated from TNG properties on disk when loaded with yt; for ray files, these properties will already be calculated and saved to disk. When you use load() on them the properties are needlessly recalculated, but nothing bad happens
    f.close()
  ray.close()
  if interactive==True:
    return load(complete_filename)
  else:
    return complete_filename
  #return complete_filename
