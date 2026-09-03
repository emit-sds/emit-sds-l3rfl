# Earth Surface Mineral Dust Source Investigation (EMIT)

## Level 3 Gridded Reflectance Product User Guide

**Version:** 
**Release Date:** TBD
**JPL-D-TBD**

Jet Propulsion Laboratory
California Institute of Technology
Pasadena, California 91109

**Change Log**

| Version | Date       | Comments      |
| ------- | ---------- | ------------- |
| x.x     | YYYY-MM-DD | Initial release |

## Table of Contents

- [1 Introduction](#1-introduction)
  * [1.1 Identification](#11-identification)
  * [1.2 Overview](#12-overview)
  * [1.3 File Formats](#13-file-formats)
    + [1.3.1 Metadata Structure](#131-metadata-structure)
    + [1.3.2 Data Products](#132-data-products)
    + [1.3.3 Grid Definition](#133-grid-definition)
    + [1.3.4 Storage and Compression](#134-storage-and-compression)
  * [1.4 Product Availability](#14-product-availability)
- [2 Working with the Data](#2-working-with-the-data)
  * [2.1 Fill Values and Screening](#21-fill-values-and-screening)
  * [2.2 Cloud Probability](#22-cloud-probability)
  * [2.3 Reading the Data](#23-reading-the-data)
- [3 Gridded Product Generation](#3-gridded-product-generation)
- [4 Known Limitations](#4-known-limitations)
- [5 References](#5-references)
- [6 Acronyms](#6-acronyms)

## 1 Introduction

### 1.1 Identification

This document describes the file structure and datasets provided in the EMIT L3 Gridded Reflectance data product. The algorithms and data content are described briefly in this guide, with the purpose of providing the user with sufficient information about the content and structure of the data files to access and use the data, in addition to understanding the uncertainties involved in the products. Full algorithmic detail is given in the EMIT L3 Gridded Reflectance ATBD.

### 1.2 Overview

The EMIT Project delivers space-based measurements of surface mineralogy of the Earth's arid dust source regions. These measurements are used to initialize Earth System Models (ESM) of the dust cycle, which describe the generation, lofting, transport, and deposition of mineral dust. Earth System Models incorporate the dust cycle to estimate the impacts of mineral dust on the optical and radiative properties of the atmosphere, and a variety of environmental and ecological processes. EMIT on the ISS makes measurements over the sunlit Earth's surface in the range of ±52° latitude. EMIT-based maps of the fractional cover of surface classes is an essential product needed for analysis of the relative abundance of source minerals to address the prime mission science questions, as well as supporting additional science and applications uses.

The EMIT instrument is a Dyson imaging spectrometer that uses contiguous spectroscopic measurements in the visible to short wavelength infrared region of the spectrum to resolve absorption features of dust-forming minerals. From the instrument's focal plane array, on-board avionics reads out raw detector counts at 1.6 Gbps, then digitizes and stores this data to a high-speed Solid-State Data Recorder (SSDR). From there, the avionics software reads the raw uncompressed data, packages this data into frames of 32 instrument lines, screens for cloudy pixels within the frames, and performs a lossless 4:1 compression of the frame's science data before storing the processed, compressed data back onto the SSDR. The data is later read from the SSDR, wrapped in CCSDS packets and then formatted as ethernet packets for transmission over the International Space Station (ISS) network and downlinked to the EMIT Instrument Operation System (IOS). Once on the ground, the EMIT IOS delivers the raw ethernet data to the SDS where Level 0 processing removes the Huntsville Operations and Support Center (HOSC) ethernet headers, groups CCSDS packet streams by APID, and sorts them by course and fine time.

Level 1B processing produces calibrated radiance together with per-pixel geolocation (`loc`) and observation geometry (`obs`) files in the instrument (non-orthorectified) frame. Level 2A processing inverts the radiance to estimated surface reflectance, reflectance uncertainty, and an atmospheric state, also in the instrument frame. A separate Level 2A Mask product provides cloud, cirrus, water, and spacecraft flags together with a SpecTf cloud probability.

The L3 Gridded Reflectance product resamples those L2A and L1B results onto a regular geographic (WGS-84 latitude/longitude) grid. No new retrieval is performed at L3: the reflectance, uncertainty, atmospheric state, and observation geometry values are those produced upstream, placed onto a map grid by nearest-neighbor assignment and screened using the L2A mask.

The EMIT L3 Gridded Reflectance products are delivered as NetCDF-4 files, with quicklooks as PNG files.

### 1.3 File Formats

#### 1.3.1 Metadata Structure

EMIT is operating from the ISS, orbiting Earth approx. 16 times in a 24-hour day period. EMIT starts and stops data recording based on a surface coverage acquisition mask. The top-level metadata identifier for EMIT data is an orbit, representing a single rotation of the ISS around Earth. Within an orbit, a period of continuous data acquisition is called an orbit segment. An orbit contains multiple orbit segments, where each orbit segment can cover up to thousands of kilometers down-track, depending on the acquisition mask map. Each orbit segment is subsequently chunked into granules of 1280 lines down-track called scenes. The last scene in an orbit segment is merged into the one before, making the last scene to be between 1280 and 2560 lines down-track. Scenes, also referred to as "granules", can be downloaded as NetCDF files, and are identified by a date-time string in the file name.

#### 1.3.2 Data Products

The "EMIT L3 ......" collection (EMITL3RFL) contains three NetCDF files per granule together with a quicklook PNG (Browse), as described in Table 1-1.

**Table 1-1:** EMITL3RFL collection file list and naming convention

| File                                                          | Description             |
| ------------------------------------------------------------- | ----------------------- |
| `EMIT_L3_RFL_<VVV>_<YYYYMMDDTHHMMSS>.nc`      | Gridded reflectance     |
| `EMIT_L3_RFLUNCERT_<VVV>_<YYYYMMDDTHHMMSS>.nc`| Gridded uncertainty     |
| `EMIT_L3_OBS_<VVV>_<YYYYMMDDTHHMMSS>.nc`      | Gridded observation data|
| `EMIT_L3_RFL_<VVV>_<YYYYMMDDTHHMMSS>.png`     | Browse                  |

`<VVV>` gives the product version number, e.g., 002

`<YYYYMMDDTHHMMSS>` is a time stamp, e.g., 20220101T083015

The structure of the gridded reflectance file is described in Table 1-2.

**Table 1-2:** EMITL3RFL NetCDF File Structure

| Group                    | Field Name                    | Type    | Dimensions          | Units         | Comments                                                    |
| ------------------------ | ----------------------------- | ------- | ------------------- | ------------- | ----------------------------------------------------------- |
| Root                     | reflectance                   | float32 | (bands, lat, lon)   | unitless      | Hemispherical-Directional Reflectance Factor, fraction of 1  |
| Root                     | lon                           | float64 | (lon)               | degrees_east  | Longitude of pixel centers (WGS-84)                         |
| Root                     | lat                           | float64 | (lat)               | degrees_north | Latitude of pixel centers (WGS-84)                          |
| Root                     | crs                           | char    | scalar              | n/a           | CRS definition; holds `spatial_ref`, `crs_wkt`, `GeoTransform` |
| Root                     | bands                         | float32 | (bands)             | nm            | Wavelength centers, duplicated at root for xarray/QGIS      |
| sensor\_band\_parameters | wavelengths                   | float32 | (bands)             | nm            | Wavelength centers                                          |
| sensor\_band\_parameters | fwhm                          | float32 | (bands)             | nm            | Full width at half maximum                                  |
| sensor\_band\_parameters | good\_wavelengths             | uint8   | (bands)             | unitless      | 1 = usable, 0 = not usable (bad band list)                  |
| state\_variables         | aerosol\_optical\_thickness   | float32 | (lat, lon)          | unitless      | AOD at 550 nm from the L2A state file (`AOT550`)            |
| state\_variables         | water\_vapor                  | float32 | (lat, lon)          | g/cm^2        | Column water vapor from the L2A state file (`H2OSTR`)       |
| masks                    | cloud\_probability            | float32 | (lat, lon)          | unitless      | Gridded `SpecTf-Cloud Probability` from the L2A Mask product |

The uncertainty file follows the same structure, with a single `reflectance_uncertainty` variable (float32, `(bands, lat, lon)`, unitless) in place of `reflectance`, and without the `state_variables` and `masks` groups. Uncertainty is reported as one standard deviation, per channel, with covariance ignored.

The observation file contains the layers listed in Table 1-3, all float32 on `(lat, lon)`.

**Table 1-3:** `observation_parameters` layers in the gridded observation file

| Variable              | Units    | Description                                                                     |
| --------------------- | -------- | ------------------------------------------------------------------------------- |
| path\_length          | m        | Sensor-to-ground distance                                                       |
| to\_sensor\_azimuth   | degrees  | 0 to 360 degrees clockwise from N                                               |
| to\_sensor\_zenith    | degrees  | 0 to 90 degrees from zenith                                                     |
| to\_sun\_azimuth      | degrees  | 0 to 360 degrees clockwise from N                                               |
| to\_sun\_zenith       | degrees  | 0 to 90 degrees from zenith                                                     |
| solar\_phase          | degrees  | Angle between to-sensor and to-sun vectors in the principal plane               |
| slope                 | degrees  | Local surface slope derived from a DEM                                          |
| aspect                | degrees  | Local surface aspect, 0 to 360 degrees clockwise from N                         |
| cosine\_i             | unitless | Apparent local illumination factor from DEM slope/aspect and to-sun vector, -1 to 1 |
| utc\_time             | hours    | Decimal hours for mid-line pixels                                               |
| earth\_sun\_distance  | AU       | Earth-sun distance                                                              |

All 2D and 3D data variables carry `grid_mapping = "crs"`, so CF-aware readers will pick up the projection automatically.

#### 1.3.3 Grid Definition

Output files are on a regular grid in EPSG:4326 (geographic, WGS-84). The default pixel size is 0.00055 degrees, chosen to approximate the 60 m native EMIT ground sampling distance. Grid extents are snapped outward to whole multiples of the pixel size, so that grids from adjacent scenes at the same pixel size align exactly.

Because the grid is angular, the ground size of a pixel is constant north-south and scales with the cosine of latitude east-west. At the default pixel size this is approximately 61 m north-south everywhere, and east-west approximately 61 m at the equator, 53 m at 30°, and 38 m at 52°. Pixels are square in degrees, not on the ground. At higher latitudes the across-track grid is finer than the input sampling, so input pixels are replicated across track without adding information.

The transform is written to `crs.GeoTransform` in GDAL order, referencing the upper-left corner of the upper-left pixel. The `lon` and `lat` coordinate variables give pixel centers, offset by half a pixel from that corner. Latitude decreases with increasing row index.

#### 1.3.4 Storage and Compression

All gridded data variables share the same storage settings. Values are quantized with `least_significant_digit` set to 5, retaining five decimal places to improve compressibility. Variables are chunked at 10 x 256 x 256 in the `(bands, lat, lon)` case and 256 x 256 in the two-dimensional `(lat, lon)` case, and are written with zlib compression at complevel 1. The chunk shape spans a block of ten contiguous channels over a 256 by 256 spatial tile, which suits both spectral access at a point and spatial access within a band. The fill value is -9999 throughout.

### 1.4 Product Availability

The EMIT L3 Gridded Reflectance products are planned for distribution through the NASA Land Processes Distributed Active Archive Center (LP DAAC, <https://lpdaac.usgs.gov/>) and NASA Earthdata (<https://earthdata.nasa.gov/>).

## 2 Working with the Data

### 2.1 Fill Values and Screening

The reserved fill value throughout the L3 products is **-9999** (floating point), set as the NetCDF `_FillValue` on all gridded data variables. Fill values occur in three situations:

1. Pixels that fall outside the scene footprint on the output grid.
2. Pixels flagged by the L2A cloud screening mask, which is gridded alongside the science data and then applied to the output.
3. Pixels carrying -9999 from upstream L1B/L2A processing, for example onboard cloud screening.

Before it is applied, the gridded screening mask is filtered to drop small flagged regions, which are potential false positives.

The same filtered mask is applied to the reflectance, atmospheric state, and observation parameter layers, so the extent of valid data is common across those products.

Users must screen for -9999 before any arithmetic. Because gridding is nearest-neighbor, a fill value in the input propagates as a fill value in the output and is never averaged into a neighboring pixel.

Separately from fill screening, the `sensor_band_parameters/good_wavelengths` array flags channels where reflectance is not usable. These are the deep atmospheric water vapor absorption regions, nominally 1325–1435 nm and 1770–1962 nm, where surface reflectance cannot be meaningfully retrieved. Values in these channels are present in the file but should not be interpreted. Most users will want to subset to `good_wavelengths == 1` before analysis.

Reflectance is reported as a fraction relative to 1, not scaled to integer counts and not percent.

### 2.2 Cloud Probability

The `masks/cloud_probability` layer in the reflectance file carries the gridded SpecTf cloud probability, in the range 0 to 1, where higher values indicate greater confidence that a cloud is present. It is provided unthresholded and is independent of the screening applied in Section 2.1, so users whose cloud sensitivity requirements differ from the delivered screening can threshold it themselves.

The probability expresses the model's confidence that a pixel is a cloud pixel. Cells outside the footprint carry -9999. See the EMIT L2A Mask ATBD for the model description and calibration.

### 2.3 Reading the Data

Because the files are CF-compliant NetCDF-4 with a `crs` variable and a `GeoTransform` attribute, they can be opened directly by GDAL:

```bash
gdalinfo NETCDF:"EMIT_L3_RFL_002_20220101T083015_2530101_007.nc":reflectance
gdal_translate -b 36 -b 24 -b 12 -a_nodata -9999 \
  NETCDF:"EMIT_L3_RFL_002_20220101T083015_2530101_007.nc":reflectance rgb.tif
```

The same files open in QGIS as a multiband raster, and in Python via `xarray.open_dataset`, where the root-level `bands` variable supplies wavelength coordinates.

Statistics attributes (`STATISTICS_MINIMUM`, `STATISTICS_MAXIMUM`, `STATISTICS_MEAN`, `STATISTICS_STDDEV`, `STATISTICS_VALID_PERCENT`) are written on the `reflectance` and `reflectance_uncertainty` variables to control default display stretches in GIS software. These are fixed nominal values intended for rendering, not computed per-granule statistics, and should not be used for analysis.

## 3 Gridded Product Generation

The L3 Gridded Reflectance product is generated by `output_conversion.py` in the [emit-sds-l3rfl](https://github.com/emit-sds/emit-sds-l3rfl) repository. The script takes as input L2A reflectance, reflectance uncertainty, atmospheric state, and mask files, together with L1B location and observation files, all in ENVI format, and writes the three NetCDF files and the browse PNG.

## 5 References

- EMIT L2A Estimated Surface Reflectance and Uncertainty ATBD: <https://github.com/emit-sds/emit-sds-l2a/blob/develop/docs/EMITL2A_ATBD.md>
- EMIT L2A Mask ATBD: <https://github.com/emit-sds/emit-sds-masks/blob/develop/docs/EMIT_L2A_Mask_ATBD.md>
- EMIT L2A Mask User Guide: <https://github.com/emit-sds/emit-sds-masks/blob/develop/docs/EMIT_L2A_Mask_User_Guide.md>
- EMIT L3 Gridded Reflectance code repository: <https://github.com/emit-sds/emit-sds-l3rfl>
- ISOFIT: <https://github.com/isofit/isofit>

## 6 Acronyms

| Acronym | Definition                                       |
| ------- | ------------------------------------------------ |
| AOD     | Aerosol Optical Depth                            |
| ATBD    | Algorithm Theoretical Basis Document             |
| CF      | Climate and Forecast (metadata conventions)      |
| DAAC    | Distributed Active Archive Center                |
| DEM     | Digital Elevation Model                          |
| EMIT    | Earth Surface Mineral Dust Source Investigation  |
| ENVI    | Environment for Visualizing Images               |
| ESM     | Earth System Model                               |
| GDAL    | Geospatial Data Abstraction Library              |
| GSD     | Ground Sampling Distance                         |
| HDRF    | Hemispherical-Directional Reflectance Factor     |
| ISS     | International Space Station                      |
| LP DAAC | Land Processes Distributed Active Archive Center |
| NCEI    | National Centers for Environmental Information   |
| PNG     | Portable Network Graphics                        |
| SDS     | Science Data System                              |
| SpecTf  | Spectral Transformer Model                       |
| WGS-84  | World Geodetic System 1984                       |