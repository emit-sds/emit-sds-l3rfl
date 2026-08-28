import argparse
import json
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np
from matplotlib.path import Path
from netCDF4 import Dataset
from osgeo import osr
from PIL import Image
from scipy.spatial import cKDTree
from scipy import ndimage
from spectral.io import envi
from kerchunk.hdf import SingleHdf5ToZarr
from kerchunk.combine import MultiZarrToZarr

from emit_utils.daac_converter import add_variable, makeGlobalAttr
from emit_utils.file_checks import envi_header

osr.UseExceptions()

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
        self.coords = coords
        self.tree = cKDTree(coords, balanced_tree=False)

    def footprint(self):
        grid = self.coords.reshape(*self.input_shape, 2)
        boundary = np.concatenate([
            grid[0, :], grid[:, -1], grid[-1, ::-1], grid[::-1, 0]
        ])
        return Path(boundary)

    def query_tree(self, ul_lon, ul_lat, lr_lon, lr_lat, pixel_size):
        self.pixel_size = pixel_size
        lines = int(np.ceil((ul_lat - lr_lat) / pixel_size))
        columns = int(np.ceil((lr_lon - ul_lon) / pixel_size))
        self.output_shape = (lines, columns)

        row_idx, col_idx = np.indices(self.output_shape)
        lon = (col_idx * pixel_size + ul_lon).flatten()
        lat = (ul_lat - row_idx * pixel_size).flatten()
        dest_points = np.column_stack([lon, lat])

        distances, flat_idx = self.tree.query(dest_points, k=1, workers=-1)
        self.indices = np.unravel_index(flat_idx, self.input_shape)
        self.distances = distances.reshape(self.output_shape)

        path = self.footprint()
        inside = path.contains_points(dest_points).reshape(self.output_shape)

        half_pixel = (self.input_spacing or pixel_size) * 0.5
        self.mask = ~(inside | (self.distances <= half_pixel))

    def project_band(self, band, no_data):
        out = band[self.indices[0], self.indices[1]].reshape(self.output_shape)
        out[self.mask] = no_data
        return out

def set_statistics(var, mn=0.0, mx=1.0, mean=0.25, std=0.15, valid_pct=100.0):
    var.setncattr("STATISTICS_MINIMUM", float(mn))
    var.setncattr("STATISTICS_MAXIMUM", float(mx))
    var.setncattr("STATISTICS_MEAN", float(mean))
    var.setncattr("STATISTICS_STDDEV", float(std))
    var.setncattr("STATISTICS_VALID_PERCENT", float(valid_pct))

def create_CRS(nc_ds, out_lines, out_columns, pixel_size, geotransform):

    lon = nc_ds.createDimension("lon", out_columns)
    lat = nc_ds.createDimension("lat", out_lines)

    x = nc_ds.createVariable("lon", "f8", ("lon",))
    x[:] = np.arange(out_columns) * pixel_size + geotransform[0] + pixel_size / 2
    x.standard_name = "longitude"
    x.long_name = "Longitude (WGS-84)"
    x.units = "degrees_east"
    x.axis = "X"

    y = nc_ds.createVariable("lat", "f8", ("lat",))
    y[:] = geotransform[3] - np.arange(out_lines) * pixel_size - pixel_size / 2
    y.standard_name = "latitude"
    y.long_name = "Latitude (WGS-84)"
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



def copy_attributes(src_obj, dst_obj):
    for k, v in src_obj.attrs.items():
        if k in ["CLASS", "NAME", "REFERENCE_LIST", "DIMENSION_LIST"]:
            continue
        dst_obj.attrs[k] = v


def filter_kwargs(src_ds):
    kwargs = {}
    if src_ds.compression is not None:
        kwargs["compression"] = src_ds.compression
        if src_ds.compression_opts is not None:
            kwargs["compression_opts"] = src_ds.compression_opts
    if src_ds.shuffle:
        kwargs["shuffle"] = True
    if src_ds.fletcher32:
        kwargs["fletcher32"] = True
    if src_ds.scaleoffset is not None:
        kwargs["scaleoffset"] = src_ds.scaleoffset
    return kwargs


def finalize_layout(path):
    """Rewrite a NetCDF-4/HDF5 file so all group/variable metadata lands near byte 0.

    netCDF writes headers at current file position when each objects is created

    this function uses hdf5 to rewrite by passing through the file twice.  First pass
    creates everything without data, second pass copies data right back in, adjusting
    the order. Bit of a hack, as it relies on hitting the disk multiple times,
    but works with the libraries in hand.
    
    Important - keep libver at v108 (HDF5 1.8), for superblock verison 2.  Version 3
    is better for reads, but is incompatible with QGIS.  A sidecar will bridge the
    gap for cloud streaming.
    """
    tmp_path = path + ".finalize.tmp"

    with h5py.File(path, "r") as fs, h5py.File(tmp_path, "w", libver="v108") as fd:
        copy_attributes(fs, fd)  # global/root attrs

        groups, datasets = [], []

        def visit(name, obj):
            if isinstance(obj, h5py.Group):
                groups.append(name)
            elif isinstance(obj, h5py.Dataset):
                datasets.append(name)

        fs.visititems(visit)

        for g in groups:
            copy_attributes(fs[g], fd.create_group(g))

        # pass 1, structure and attributes
        layouts = {}
        for name in datasets:
            src_ds = fs[name]
            layout = src_ds.id.get_create_plist().get_layout()  # 1=contiguous, 2=chunked
            layouts[name] = layout
            if layout == 2:
                dst_ds = fd.create_dataset(name, shape=src_ds.shape, dtype=src_ds.dtype,
                                           chunks=src_ds.chunks, fillvalue=src_ds.fillvalue,
                                           **filter_kwargs(src_ds))
            else:
                dst_ds = fd.create_dataset(name, shape=src_ds.shape, dtype=src_ds.dtype,
                                           fillvalue=src_ds.fillvalue)
            copy_attributes(src_ds, dst_ds)

        # Recreate dimension-scale relationships (reference-typed, can't be copied as attrs).
        for name in datasets:
            if h5py.h5ds.is_scale(fs[name].id):
                h5py.h5ds.set_scale(fd[name].id, name.encode())
        for name in datasets:
            src_ds = fs[name]
            for i in range(len(src_ds.dims)):
                for j in range(len(src_ds.dims[i])):
                    scale_name = src_ds.dims[i][j].name.lstrip("/")
                    h5py.h5ds.attach_scale(fd[name].id, fd[scale_name].id, i)

        # pass 2 (small first)
        def order_key(name):
            src_ds = fs[name]
            return int(np.prod(src_ds.shape)) * src_ds.dtype.itemsize

        contiguous = [n for n in datasets if layouts[n] != 2]
        chunked = sorted((n for n in datasets if layouts[n] == 2), key=order_key)
        write_order = contiguous + chunked

        with open(path, "rb") as raw:
            for name in write_order:
                src_ds, dst_ds = fs[name], fd[name]
                if layouts[name] == 2:
                    for i in range(src_ds.id.get_num_chunks()):
                        ci = src_ds.id.get_chunk_info(i)
                        raw.seek(ci.byte_offset)
                        dst_ds.id.write_direct_chunk(ci.chunk_offset, raw.read(ci.size),
                                                     filter_mask=ci.filter_mask)
                else:
                    dst_ds[...] = src_ds[...]

    # temp file becomes actual
    os.replace(tmp_path, path)
    logging.info(f"finalized layout: {os.path.basename(path)}")


def write_combined_sidecar(paths, sidecar_path, url_basename="lp-prod-protected/EMITL2ARFL.002"):
    """
    Write an integrated kerchunk sidecar json file on the common grid,
    with template URL resolution.

    paths[0] must be the base (e.g. reflectance)
    """
    
    # Generate dictionaries - may have lingering basenames
    per_file = [SingleHdf5ToZarr(p, inline_threshold=0).translate()
                for p in paths]
    
    combined = MultiZarrToZarr(per_file, concat_dims=[],
                               identical_dims=["lat", "lon", "bands"]).translate()
    
    old_templates = combined.pop("templates", {})
    
    # map using templates
    # e.g., "emit20220818t020924-l3-obs.nc" -> "{{u}}/emit20220818t020924-l3-obs.nc"
    name_to_template = {os.path.basename(p): f"{{{{u}}}}/{os.path.basename(p)}" for p in paths}
    
    # Walk through the final references and forcefully update the URLs
    for k, v in combined.get("refs", {}).items():
        # Check if v is a chunk reference [url, offset, size]
        if isinstance(v, list) and len(v) >= 1 and isinstance(v[0], str):
            url_val = v[0]
            
            # expand Kerchunk
            for tk, tv in old_templates.items():
                url_val = url_val.replace(f"{{{{{tk}}}}}", tv)
                
            # replace literals with url
            for basename, new_url in name_to_template.items():
                if basename in url_val:
                    url_val = new_url
                    break
                    
            v[0] = url_val
            
    # inject our templates to the root of the sidecar payload
    url_basename = url_basename.strip('/')
    fid = os.path.splitext(os.path.basename(paths[0]))[0]

    combined["templates"] = {
        "u": f"s3://{url_basename}/{fid}",
        "u_https_hint": f"https://data.lpdaac.earthdatacloud.nasa.gov/{url_basename}/{fid}"
    }

    with open(sidecar_path, "w") as f:
        json.dump(combined, f, indent=2)
    logging.info(f"combined sidecar: {os.path.basename(sidecar_path)}")

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''This script \
    converts L1B G and L2A PGE outputs to L3 DAAC compatible formats, with supporting metadata''', add_help=True)

    parser.add_argument('rfl_output_filename', type=str, help="Output Reflectance netcdf filename")
    parser.add_argument('rfl_unc_output_filename', type=str, help="Output Reflectance Uncertainty netcdf filename")
    parser.add_argument('obs_output_filename', type=str, help="Output Observables netcdf filename")
    parser.add_argument('browse_output_filename', type=str, help="Output browse image filename")
    parser.add_argument('rfl_file', type=str, help="EMIT L2A reflectance ENVI file")
    parser.add_argument('rfl_unc_file', type=str, help="EMIT L2A reflectance uncertainty ENVI file")
    parser.add_argument('state_file', type=str, help="EMIT L2A reflectance state ENVI file")
    parser.add_argument('mask_file', type=str, help="EMIT L2A water/cloud mask ENVI file")
    parser.add_argument('loc_file', type=str, help="EMIT L1B location data ENVI file")
    parser.add_argument('obs_file', type=str, help="EMIT L1B observables data ENVI file")
    parser.add_argument('version', type=str, help="3 digit (with leading V) version number")
    parser.add_argument('software_delivery_version', type=str, help="The extended build number at delivery time")
    parser.add_argument('--log_file', type=str, default=None, help="Logging file to write to")
    parser.add_argument('--log_level', type=str, default="INFO", help="Logging level")
    parser.add_argument('--chunksize', type=int, nargs=3, default=None, help="Chunk size for netCDF compression as (bands, lat, lon)")
    parser.add_argument('--complevel', type=int, default=1, help="netCDF compression level (1-9)")
    parser.add_argument('--compress', action='store_true', default=False, help="Enable zlib compression")
    parser.add_argument('--pixel_size', type=float, default=0.00055, help="Pixel size for the output grid")
    parser.add_argument('--mask_band', type=int, default=5, help="Band index to apply for mask")
    parser.add_argument('--prob_band', type=int, default=4, help="Band index for specTf cloud probability")
    parser.add_argument('--max_workers', type=int, default=64, help="Maximum number of workers to used")
    parser.add_argument('--sidecar', action='store_true', default=False,
                        help="Write a single, combined kerchunk reference JSON that covers all three netCDFs")

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

    t0 = time.time()
    grid = Gridder()
    grid.create_tree(coords, (in_lines, in_columns), input_spacing=input_spacing)
    grid.query_tree(lon_min, lat_max, lon_max, lat_min, ps)
    logging.info(f"gridder: {time.time() - t0:.3f}s")

    geotransform = (lon_min - (ps/2), ps, 0, lat_max + (ps/2), 0, -ps)

    out_lines, out_columns = grid.output_shape

    # Create reflectance netCDF
    ###########################
    logging.info(f'Creating netCDF4 file: {args.rfl_output_filename}')

    rfl_ds = envi.open(envi_header(args.rfl_file))
    mask_ds = envi.open(envi_header(args.mask_file))
    mask = mask_ds.open_memmap(interleave='bip')[..., args.mask_band].copy()

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
    wl = np.array([float(d) for d in rfl_ds.metadata['wavelength']])
    add_variable(nc_ds, "sensor_band_parameters/wavelengths", "f4", "Wavelength Centers", "nm",
                 wl, {"dimensions": ("bands",)}, standard_name = "radiation_wavelength")
    add_variable(nc_ds, "sensor_band_parameters/fwhm", "f4", "Full Width at Half Max", "nm",
                 [float(d) for d in rfl_ds.metadata['fwhm']], {"dimensions": ("bands",)})
    add_variable(nc_ds, "sensor_band_parameters/good_wavelengths", "u1", "Wavelengths where reflectance is useable: 1 = good data, 0 = bad data", "unitless",
                 bbl, {"dimensions": ("bands",)})
    # It's redundant, but also add 'bands' at root, for xarray and QGIS
    add_variable(nc_ds, "bands", "f4", "Wavelength Centers", "nm", [float(d) for d in rfl_ds.metadata['wavelength']],
                 {"dimensions": ("bands",)}, standard_name = "radiation_wavelength")

    logging.debug('Gridding and writing refectance data')

    t0 = time.time()
    rfl_cube = np.array(rfl_ds.open_memmap(interleave='bip'))
    logging.info(f"read cube: {time.time() - t0:.3f}s")

    mask_grid = grid.project_band(mask, -9999)

    lbl, n = ndimage.label(mask_grid)
    eroded = ndimage.binary_erosion(mask_grid, np.ones((3, 3), bool))
    keep = np.zeros(n + 1, bool)
    keep[np.unique(lbl[eroded])] = True
    keep[0] = False

    mask_grid = keep[lbl]

    rfl_grid = np.zeros((rfl_ds.nbands, out_lines, out_columns), dtype=np.float32)
    def _proj(band):
        rfl_grid[band] = grid.project_band(rfl_cube[:,:,band], -9999)

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        list(ex.map(_proj, range(rfl_ds.nbands)))

    rfl_grid[:,mask_grid == 1] = -9999
    logging.info(f"band gridding: {time.time() - t0:.3f}s")

    t0 = time.time()
    browse_bands = [int(np.argmin(np.abs(wl-x))) for x in [660,550,440]]
    browse = rfl_grid[browse_bands]
    browse[browse == -9999] = np.nan

    lo, hi = np.nanpercentile(browse, [2, 98])
    browse = np.clip((browse - lo) / (hi - lo), 0, 1)
    browse = (browse * 255).astype(np.uint8).transpose(1, 2, 0)

    Image.fromarray(browse).save(args.browse_output_filename)
    logging.info(f"browse: {time.time() - t0:.3f}s")

    kargs = {'zlib': args.compress,
             'complevel': args.complevel,
             'fill_value': -9999,
             'least_significant_digit': 5}

    if args.chunksize is not None:
        kargs['chunksizes'] = tuple(args.chunksize)

    t0 = time.time()
    add_variable(nc_ds, "reflectance", "f4", "Hemispherical-Directional Reflectance Factor",
                 "unitless", rfl_grid, {"dimensions": ("bands", 'lat','lon'), **kargs},
                 standard_name = "surface_bidirectional_reflectance")
    logging.info(f"rfl write: {time.time() - t0:.3f}s")

    nc_ds["reflectance"].grid_mapping = "crs"

    set_statistics(nc_ds["reflectance"])

    logging.debug('Gridding and writing state data')
    t0 = time.time()
    state_obj = envi.open(envi_header(args.state_file))
    state_ds = state_obj.open_memmap(interleave='bip')

    band_names = [b.strip() for b in state_obj.metadata['band names']]
    aot_band = band_names.index('AOT550')
    h2o_band = band_names.index('H2OSTR')

    aot550 = state_ds[..., aot_band].astype(np.float32)
    aot550 = grid.project_band(aot550, -9999)
    aot550[mask_grid == 1] = -9999

    # Change chunksize for 2D case (state variables below are all 2D)
    if args.chunksize is not None:
        kargs['chunksizes'] = tuple(args.chunksize[-2:])

    add_variable(nc_ds,
                 "state_variables/aerosol_optical_thickness",
                 "f4", "Optical thickness of atmosphere layer due to ambient aerosol particles at 550 nm",
                 'unitless',
                 aot550,
                 {"dimensions": ("lat", "lon"), **kargs},
                 standard_name = "atmosphere_absorption_optical_thickness_due_to_ambient_aerosol_particles"
                 )

    nc_ds["state_variables/aerosol_optical_thickness"].grid_mapping = "crs"

    wv = state_ds[..., h2o_band].astype(np.float32)
    wv = grid.project_band(wv, -9999)
    wv[mask_grid == 1] = -9999

    add_variable(nc_ds,
                 "state_variables/water_vapor",
                 "f4",
                 "Atmospheric mass content of water vapor",
                 "g/cm^2",
                 wv,
                 {"dimensions": ("lat", "lon"), **kargs},
                 standard_name = "atmosphere_mass_content_of_water_vapor")

    nc_ds["state_variables/water_vapor"].grid_mapping = "crs"
    logging.info(f"state write: {time.time() - t0:.3f}s")

    cloud_prob = mask_ds.open_memmap(interleave='bip')[..., args.prob_band].copy()
    cloud_prob = grid.project_band(cloud_prob, -9999)

    add_variable(nc_ds,
                 "masks/cloud_probability",
                 "f4",
                 "SpecTf-Cloud Probability",
                 "unitless",
                 cloud_prob,
                 {"dimensions": ("lat", "lon"), **kargs})

    nc_ds["masks/cloud_probability"].grid_mapping = "crs"

    nc_ds.ncei_template_version = "NCEI_NetCDF_Grid_Template_v2.0"
    nc_ds.close()
    logging.debug(f'Successfully created {args.rfl_output_filename}')

    finalize_layout(args.rfl_output_filename)
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
                 wl, {"dimensions": ("bands",)}, standard_name = "radiation_wavelength")
    add_variable(nc_ds, "sensor_band_parameters/fwhm", "f4", "Full Width at Half Max", "nm",
                 [float(d) for d in rfl_ds.metadata['fwhm']], {"dimensions": ("bands",)})
    add_variable(nc_ds, "sensor_band_parameters/good_wavelengths", "u1", "Wavelengths where reflectance is useable: 1 = good data, 0 = bad data", "unitless",
                 bbl, {"dimensions": ("bands",)})
    # It's redundant, but also add 'bands' at root, for xarray and QGIS
    add_variable(nc_ds, "bands", "f4", "Wavelength Centers", "nm", [float(d) for d in rfl_ds.metadata['wavelength']],
                 {"dimensions": ("bands",)}, standard_name = "radiation_wavelength")

    logging.debug('Gridding and writing refectance uncertainty data')

    t0 = time.time()
    rfl_unc_cube = np.array(rfl_unc_ds.open_memmap(interleave='bip'))

    logging.info(f"read cube: {time.time() - t0:.3f}s")

    rfl_unc_grid = np.zeros((rfl_unc_ds.nbands, out_lines, out_columns), dtype=np.float32)

    def _proj(band):
        rfl_unc_grid[band] = grid.project_band(rfl_unc_cube[:,:,band], -9999)

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        list(ex.map(_proj, range(rfl_unc_ds.nbands)))

    rfl_unc_grid[:,mask_grid == 1] = -9999

    logging.info(f"band gridding: {time.time() - t0:.3f}s")

    #Chunksize back to 3d case
    if args.chunksize is not None:
        kargs['chunksizes'] = tuple(args.chunksize)

    t0 = time.time()
    add_variable(nc_ds, "reflectance_uncertainty", "f4", "Hemispherical-Directional Reflectance Factor Uncertainty",
                 "unitless", rfl_unc_grid, {"dimensions": ("bands", 'lat','lon'), **kargs},)
    logging.info(f"rfl unc write: {time.time() - t0:.3f}s")

    nc_ds["reflectance_uncertainty"].grid_mapping = "crs"
    set_statistics(nc_ds["reflectance_uncertainty"])

    nc_ds.ncei_template_version = "NCEI_NetCDF_Grid_Template_v2.0"

    nc_ds.close()
    logging.debug(f'Successfully created {args.rfl_unc_output_filename}')

    finalize_layout(args.rfl_unc_output_filename)
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

    logging.debug('Gridding and writing observation data')

    # Change chunksize for 2D case
    if args.chunksize is not None:
        kargs['chunksizes'] = tuple(args.chunksize[-2:])

    for band in obs_metadata.keys():

        band_gridded = grid.project_band(obs_mmap[:,:,obs_metadata[band]['band']],-9999)
        band_gridded[mask_grid == 1] = -9999
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

    finalize_layout(args.obs_output_filename)
    logging.debug(f'Successfully created {args.obs_output_filename}')

    if args.sidecar:
        write_combined_sidecar(
            [args.rfl_output_filename, args.rfl_unc_output_filename, args.obs_output_filename],
            os.path.splitext(args.rfl_output_filename)[0] + ".json")


if __name__ == "__main__":
    main()