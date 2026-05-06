import argparse
import logging

from osgeo import osr
import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree
from spectral.io import envi

from emit_utils.daac_converter import add_variable, makeGlobalAttr
from emit_utils.file_checks import envi_header


obs_metadata = {"path_length": {"standard_name":  None,
                                "long_name": "Path length",
                                "description": "Sensor-to-ground distance",
                                "units": "m",
                                "band": 0},
                "to_sensor_azimuth": {"standard_name":  "sensor_azimuth_angle",
                                      "long_name": "To-sensor-azimuth",
                                      "description": "0 to 360 degrees clockwise from N",
                                      "units": "degrees",
                                      "band": 1},
                "to_sensor_zenith": {"standard_name":  "sensor_zenith_angle",
                                     "long_name": "To-sensor-zenith",
                                     "description": "0 to 90 degrees from zenith",
                                     "units": "degrees",
                                     "band": 2},
                "to_sun_azimuth": {"standard_name": "solar_azimuth_angle",
                                   "long_name": "To-sun-azimuth",
                                   "description": "0 to 360 degrees clockwise from N",
                                   "units": "degrees",
                                   "band": 3},
                "to_sun_zenith": {"standard_name":  "solar_zenith_angle",
                                  "long_name": "To-sun-zenith",
                                   "description":"0 to 360 degrees clockwise from N",
                                   "units": "degrees",
                                   "band": 4},
                "solar_phase": {"standard_name":  None,
                                "long_name": "Solar phase",
                                "description":"Degrees between to-sensor and to-sun vectors in principal plane",
                                "units": "degrees",
                                "band": 5},
                "slope": {"standard_name":  "ground_slope_angle",
                          "long_name": "Slope",
                          "description": "Local surface slope as derived from DEM in degrees",
                          "units": "degrees",
                          "band": 6},
                "aspect": {"standard_name":"ground_slope_direction",
                           "long_name": "Aspect",
                           "description": "Local surface aspect 0 to 360 degrees clockwise from N",
                           "units": "degrees",
                           "band": 7},
                "cosine_i": {"standard_name":  None,
                             "long_name": "Cosine i",
                           "description": "Apparent local illumination factor based on DEM slope and aspect and to sun vector, -1 to 1",
                           "units": "unitless",
                           "band": 8},
                "utc_time": {"standard_name":  "time",
                             "long_name": "UTC time",
                             "description": "Decimal hours for mid-line pixels",
                             "units": "hours",
                             "band": 9},
                "earth_sun_distance": {"standard_name":  None,
                                       "long_name": "Earth sun distance",
                                       "description": "Earth-sun distance",
                                       "units": "AU",
                                       "band": 10}}

class Gridder():
    """Nearest Neighbor Gridder class

    TODO: Will fail across antimeridian, need to update

    """
    def __init__(self):
        self.tree = None
        self.indices = None
        self.distances = None
        self.pixel_size = None
        self.input_spacing = None
        self.input_shape = None
        self.output_shape = None
        self.mask = None

    def create_tree(self, coords, input_shape, input_spacing=None):
        self.input_shape = input_shape
        self.input_spacing = input_spacing
        self.tree = cKDTree(coords, balanced_tree=False)

    def query_tree(self, ul_lon, ul_lat, lr_lon, lr_lat, pixel_size):
        self.pixel_size = pixel_size
        lines = int(np.ceil((ul_lat - lr_lat) / pixel_size))
        columns = int(np.ceil((lr_lon - ul_lon) / pixel_size))
        self.output_shape = (lines, columns)

        row_idx, col_idx = np.indices(self.output_shape)
        lon = (col_idx * pixel_size + ul_lon).flatten()
        lat = (ul_lat - row_idx * pixel_size).flatten()
        dest_points = np.column_stack([lon, lat])

        distances, flat_idx = self.tree.query(dest_points, k=1)
        self.indices = np.unravel_index(flat_idx, self.input_shape)
        self.distances = distances.reshape(self.output_shape)

        #TODO: May need to update logic to make sure no gaps inside image
        threshold = (self.input_spacing or pixel_size) * np.sqrt(2) / 2
        self.mask = ~(self.distances < threshold)

    def project_band(self, band, no_data):
        out = band[self.indices[0], self.indices[1]].reshape(self.output_shape)
        out[self.mask] = no_data
        return out

def create_CRS(nc_ds, out_lines, out_columns, pixel_size, geotransform):

    lon = nc_ds.createDimension("lon", out_columns)
    lat = nc_ds.createDimension("lat", out_lines)

    x = nc_ds.createVariable("lon", "f8", ("lon",))
    x[:] = np.arange(out_columns) * pixel_size + geotransform[0] + pixel_size / 2
    x.standard_name = "longitude"
    x.long_name = "longitude"
    x.units = "degrees_east"
    x.axis = "X"

    y = nc_ds.createVariable("lat", "f8", ("lat",))
    y[:] = geotransform[3] - np.arange(out_lines) * pixel_size - pixel_size / 2
    y.standard_name = "latitude"
    y.long_name = "latitude"
    y.units = "degrees_north"
    y.axis = "Y"

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    crs = nc_ds.createVariable("crs", "S1", ())
    crs.grid_mapping_name = "latitude_longitude"
    crs.long_name = "CRS definition"
    crs.longitude_of_prime_meridian = 0.0
    crs.semi_major_axis = srs.GetSemiMajor()
    crs.inverse_flattening = srs.GetInvFlattening()
    crs.spatial_ref = srs.ExportToWkt()
    crs.crs_wkt = srs.ExportToWkt()
    crs.GeoTransform = " ".join(str(x) for x in geotransform)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''This script \
    converts L1B G and L2A PGE outputs to L3 DAAC compatible formats, with supporting metadata''', add_help=True)

    parser.add_argument('rfl_output_filename', type=str, help="Output Reflectance netcdf filename")
    parser.add_argument('rfl_unc_output_filename', type=str, help="Output Reflectance Uncertainty netcdf filename")
    parser.add_argument('obs_output_filename', type=str, help="Output Observables netcdf filename")
    parser.add_argument('rfl_file', type=str, help="EMIT L2A reflectance ENVI file")
    parser.add_argument('rfl_unc_file', type=str, help="EMIT L2A reflectance uncertainty ENVI file")
    parser.add_argument('state_file', type=str, help="EMIT L2A reflectance state ENVI file")
    parser.add_argument('mask_file', type=str, help="EMIT L2A water/cloud mask ENVI file")
    parser.add_argument('loc_file', type=str, help="EMIT L1B location data ENVI file")
    parser.add_argument('obs_file', type=str, help="EMIT L1B observables data ENVI file")
    parser.add_argument('version', type=str, help="3 digit (with leading V) version number")
    parser.add_argument('software_delivery_version', type=str, help="The extended build number at delivery time")
    parser.add_argument('--ummg_file', type=str, help="Output UMMG filename")
    parser.add_argument('--log_file', type=str, default=None, help="Logging file to write to")
    parser.add_argument('--log_level', type=str, default="INFO", help="Logging level")
    parser.add_argument('--chunksize', type=int, nargs=3, default=None, help="Chunk size for netCDF compression as (bands, lat, lon)")
    parser.add_argument('--complevel', type=int, default=1, help="netCDF compression level (1-9)")
    parser.add_argument('--compress', action='store_true', default=False, help="Enable zlib compression")
    args = parser.parse_args()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=args.log_level, filename=args.log_file)

    logging.info(f'Creating gridder')
    loc_ds = envi.open(envi_header(args.loc_file))
    obs_ds = envi.open(envi_header(args.obs_file))

    loc = loc_ds.open_memmap(interleave='bip')
    longitude = loc[..., 0].copy()
    latitude = loc[..., 1].copy()
    in_lines, in_columns = latitude.shape

    ps = args.pixel_size
    lon_min = np.floor(longitude.min() / ps) * ps
    lon_max = np.ceil(longitude.max() / ps) * ps
    lat_min = np.floor(latitude.min() / ps) * ps
    lat_max = np.ceil(latitude.max() / ps) * ps

    input_spacing = max((longitude.max() - longitude.min()) / in_columns,
                        (latitude.max() - latitude.min()) / in_lines)

    coords = np.column_stack([longitude.ravel(), latitude.ravel()])

    grid = Gridder()
    grid.create_tree(coords, (in_lines, in_columns), input_spacing=input_spacing)
    grid.query_tree(lon_min, lat_max, lon_max, lat_min, ps)

    geotransform = (lon_min - (ps/2), ps, 0, lat_max + (ps/2), 0, -ps)

    out_lines, out_columns = grid.output_shape

    # Create reflectance netCDF
    ###########################

    logging.info(f'Creating netCDF4 file: {args.rfl_output_filename}')

    rfl_ds = envi.open(envi_header(args.rfl_file))
    mask_ds = envi.open(envi_header(args.mask_file))
    mask = mask_ds.open_memmap(interleave='bip')[..., 9].copy() != 0

    nc_ds = Dataset(args.rfl_output_filename, 'w', clobber=True, format='NETCDF4')

    makeGlobalAttr(nc_ds, args.rfl_file, args.version)

    nc_ds.title = "EMIT L3 Gridded Estimated Surface Reflectance " + args.version
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis file contains L3 estimated surface reflectances \
and geolocation data. Reflectance estimates are created using an Optimal Estimation technique - see ATBD for \
details. Reflectance values are reported as fractions (relative to 1)."

    create_CRS(nc_ds, out_lines, out_columns, args.pixel_size, geotransform)

    # Handle data pre January, where bbl was not set in ENVI header
    if 'bbl' not in rfl_ds.metadata or rfl_ds.metadata['bbl'] == '{}':
        wl = np.array(nc_ds['sensor_band_parameters']['wavelengths'])
        bbl = np.ones(len(wl))
        bbl[np.logical_and(wl > 1325, wl < 1435)] = 0
        bbl[np.logical_and(wl > 1770, wl < 1962)] = 0
    else:
        bbl = [bool(d) for d in rfl_ds.metadata['bbl']]


    bands = nc_ds.createDimension("bands", rfl_ds.nbands)
    add_variable(nc_ds, "sensor_band_parameters/wavelengths", "f4", "Wavelength Centers", "nm",
                 [float(d) for d in rfl_ds.metadata['wavelength']], {"dimensions": ("bands",)})
    add_variable(nc_ds, "sensor_band_parameters/fwhm", "f4", "Full Width at Half Max", "nm",
                 [float(d) for d in rfl_ds.metadata['fwhm']], {"dimensions": ("bands",)})
    add_variable(nc_ds, "sensor_band_parameters/good_wavelengths", "u1", "Wavelengths where reflectance is useable: 1 = good data, 0 = bad data", "unitless",
                 bbl, {"dimensions": ("bands",)})

    logging.debug('Gridding and writing refectance data')

    rfl_mmap = rfl_ds.open_memmap(interleave='bip')[...].copy()
    rfl_mmap[mask,:] = -9999
    rfl_grid = np.zeros((rfl_ds.nbands,out_lines,out_columns), dtype = np.float32)
    for band in range(rfl_ds.nbands):
        rfl_grid[band] = grid.project_band(rfl_mmap[:,:,band],-9999)

    kargs = {'zlib': args.compress,
             'complevel': args.complevel,
             'fill_value': -9999,
             'significant_digits': 5}

    if args.chunksize is not None:
        kargs['chunksizes'] = tuple(args.chunksize)

    add_variable(nc_ds,
                 "reflectance",
                 "f4",
                 "Surface hemispherical directional reflectance factor",
                 "unitless",
                 rfl_grid,
                 {"dimensions": ("bands", 'lat','lon'), **kargs},)

    nc_ds["reflectance"].grid_mapping = "crs"

    logging.debug('Gridding and writing state data')
    state_ds = envi.open(envi_header(args.state_file)).open_memmap(interleave='bip')

    #TODO: Check band names in case layer order is not constant
    aot550 = state_ds[..., 0].astype(np.float32)
    aot550[mask] =  -9999
    aot550 = grid.project_band(aot550, -9999)

    add_variable(nc_ds,
                 "state_variables/aerosol_optical_thickness",
                 "f4", "Optical thickness of atmosphere layer due to ambient aerosol particles",
                 'unitless',
                 aot550,
                 {"dimensions": ("lat", "lon"), **kargs},
                 )

    nc_ds["state_variables/aerosol_optical_thickness"].grid_mapping = "crs"

    wv = state_ds[..., 1].astype(np.float32)
    wv[mask] =  -9999
    wv = grid.project_band(wv, -9999)

    # Change chunksize for 2D case
    if args.chunksize is not None:
        kargs['chunksizes'] = tuple(args.chunksize[-2:])

    add_variable(nc_ds,
                 "state_variables/water_vapor",
                 "f4",
                 "LWE thickness of atmosphere mass content of water vapor",
                 "cm",
                 wv,
                 {"dimensions": ("lat", "lon"), **kargs},)

    nc_ds["state_variables/water_vapor"].grid_mapping = "crs"
    nc_ds.ncei_template_version = "NCEI_NetCDF_Grid_Template_v2.0"
    nc_ds.close()
    logging.debug(f'Successfully created {args.rfl_output_filename}')

    # Create uncertainty netCDF
    ###########################

    logging.info(f'Creating netCDF4 file: {args.rfl_unc_output_filename}')
    rfl_unc_ds = envi.open(envi_header(args.rfl_unc_file))
    nc_ds = Dataset(args.rfl_unc_output_filename, 'w', clobber=True, format='NETCDF4')

    makeGlobalAttr(nc_ds, args.rfl_unc_file, args.version)

    nc_ds.title = "EMIT L3 Gridded Estimated Surface Reflectance Uncertainty " + args.version
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis file contains L3 estimated surface reflectances \
and geolocation data. Reflectance uncertainty estimates are created using an Optimal Estimation technique - see ATBD for \
details. Reflectance uncertainty values are reported as fractions (relative to 1)."

    create_CRS(nc_ds, out_lines, out_columns, args.pixel_size, geotransform)

    bands = nc_ds.createDimension("bands", rfl_unc_ds.nbands)
    add_variable(nc_ds, "sensor_band_parameters/wavelengths", "f4", "Wavelength Centers", "nm",
                 [float(d) for d in rfl_ds.metadata['wavelength']], {"dimensions": ("bands",)})
    add_variable(nc_ds, "sensor_band_parameters/fwhm", "f4", "Full Width at Half Max", "nm",
                 [float(d) for d in rfl_ds.metadata['fwhm']], {"dimensions": ("bands",)})
    add_variable(nc_ds, "sensor_band_parameters/good_wavelengths", "u1", "Wavelengths where reflectance is useable: 1 = good data, 0 = bad data", "unitless",
                 bbl, {"dimensions": ("bands",)})

    logging.debug('Gridding and writing refectance uncertainty data')
    rfl_unc_mmap = rfl_unc_ds.open_memmap(interleave='bip')[...].copy()
    rfl_unc_mmap[mask,:] = -9999
    rfl_unc_grid = np.zeros((rfl_unc_ds.nbands,out_lines,out_columns), dtype = np.float32)
    for band in range(rfl_unc_ds.nbands):
        rfl_unc_grid[band] = grid.project_band(rfl_unc_mmap[:,:,band],-9999)

    #Chunksize back to 3d case
    if args.chunksize is not None:
        kargs['chunksizes'] = tuple(args.chunksize)

    add_variable(nc_ds,
                 "reflectance_uncertainty",
                 "f4",
                 "Surface hemispherical directional reflectance factor uncertainty",
                 "unitless",
                 rfl_unc_grid,
                 {"dimensions": ("bands", 'lat','lon'), **kargs},)

    nc_ds.ncei_template_version = "NCEI_NetCDF_Grid_Template_v2.0"

    nc_ds.close()
    logging.debug(f'Successfully created {args.rfl_unc_output_filename}')

    # Create observations netCDF
    ###########################
    logging.info(f'Creating netCDF4 file: {args.obs_output_filename}')

    obs_ds = envi.open(envi_header(args.obs_file))
    nc_ds = Dataset(args.obs_output_filename, 'w', clobber=True, format='NETCDF4')

    makeGlobalAttr(nc_ds, args.obs_file, args.version)

    nc_ds.title = "EMIT L3 Gridded Observation Data " + args.version
    nc_ds.summary = nc_ds.summary + \
        f"\\n\\nThis file contains L3 geometric information (path length, view and solar angles, timing) associated with \
each pixel in an acquisition."

    create_CRS(nc_ds, out_lines, out_columns, args.pixel_size, geotransform)

    obs_mmap = obs_ds.open_memmap(interleave='bip')[...].copy()
    obs_mmap[mask,:] = -9999

    for band in range(rfl_unc_ds.nbands):
        rfl_unc_grid[band] = grid.project_band(rfl_unc_mmap[:,:,band],-9999)

    logging.debug('Gridding and writing observation data')

    # Change chunksize for 2D case
    if args.chunksize is not None:
        kargs['chunksizes'] = tuple(args.chunksize[-2:])

    for band in obs_metadata.keys():

        band_gridded = grid.project_band(obs_mmap[:,:,obs_metadata[band]['band']],-9999)
        variable_name = f"observation_parameters/{band}"

        add_variable(nc_ds,
                     variable_name,
                     "f4",
                     obs_metadata[band]['long_name'],
                     obs_metadata[band]['units'],
                     band_gridded,
                     {"dimensions": ('lat', 'lon'), **kargs},)


        if obs_metadata[band]['standard_name']:
            nc_ds[variable_name].standard_name = obs_metadata[band]['standard_name']

        nc_ds[variable_name].grid_mapping = "crs"

    nc_ds.ncei_template_version = "NCEI_NetCDF_Grid_Template_v2.0"
    nc_ds.close()

    logging.debug(f'Successfully created {args.obs_output_filename}')


