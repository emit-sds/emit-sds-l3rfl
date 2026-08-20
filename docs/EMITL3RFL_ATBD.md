JPL D-TBD
ATBD-EMIT-L3RFL

# Earth Surface Mineral dust source InvesTigation (EMIT)

## EMIT L3 Gridded Estimated Surface Reflectance and Uncertainty

### Algorithm Theoretical Basis Document

**Version:** 0.1
**Release Date:** TBD

Jet Propulsion Laboratory
California Institute of Technology
Pasadena, California 91109

**Change Log**

| Version | Date       | Comments      |
| ------- | ---------- | ------------- |
| 0.1     | 2026-08-14 | Initial Draft |


## Table of Contents

- [1. Key Team Members](#1-key-team-members)
- [2. Historical Context and Background on the EMIT Mission and its Instrumentation](#2-historical-context-and-background-on-the-emit-mission-and-its-instrumentation)
- [3. Algorithm Rationale](#3-algorithm-rationale)
- [4. Algorithm Implementation](#4-algorithm-implementation)
  * [4.1 Input Data](#41-input-data)
  * [4.2 Theoretical Description](#42-theoretical-description)
    + [4.2.1 Grid Definition](#421-grid-definition)
    + [4.2.2 Nearest-Neighbor Assignment](#422-nearest-neighbor-assignment)
    + [4.2.3 Footprint Determination and Fill](#423-footprint-determination-and-fill)
    + [4.2.4 CLoud Masking](#424-masking)
      - [4.2.4.1 Connected-Component Filtering of the Cloud Screening Mask](#4241-connected-component-filtering-of-the-screening-mask)
    + [4.2.5 Bad Band Handling](#425-bad-band-handling)
    + [4.2.6 Browse Image Generation](#426-browse-image-generation)
  * [4.3 Practical Considerations](#43-practical-considerations)
  * [4.4 Output Data](#44-output-data)
- [5. Calibration, uncertainty characterization and propagation, and validation](#5-calibration-uncertainty-characterization-and-propagation-and-validation)
- [6. Constraints and Limitations](#6-constraints-and-limitations)
- [7. Code Repository and References](#7-code-repository-and-references)

## **1. Key Team Members**

The EMIT Science Data System team at the Jet Propulsion Laboratory, California Institute of Technology, developed and maintains the L3 gridding stage. The algorithm outputs it uses are documented in the EMIT L2A Reflectance ATBD and the EMIT L2A Mask ATBD. 

## **2. Historical Context and Background on the EMIT Mission and its Instrumentation**

Mineral dust aerosols originate as soil particles lifted into the atmosphere by wind erosion. Mineral dust created by human activity makes a large contribution to the uncertainty of direct radiative forcing by anthropogenic aerosols (USGCRP and IPCC) and is a prominent aerosol constituent around the globe. Dust radiative forcing is highly dependent on its mineral-specific absorption properties, and the current range of iron oxide abundance in dust source models translates into a large range of values, even changing the sign of the forcing predicted by Earth System Models (Li et al., 2021). NASA selected the Earth Surface Mineral Dust Source Investigation (EMIT) to close this knowledge gap, launching an imaging spectrometer to the International Space Station to directly measure and map the soil mineral composition of critical dust-forming regions worldwide.

EMIT data products are organized into levels. Level 1B provides calibrated radiance at sensor along with per-pixel geolocation, observation geometry, and a geolocation lookup table. Level 2A inverts radiance to estimated surface reflectance and uncertainty using an Optimal Estimation approach, and produces the atmospheric state and data masks. Both are delivered in the instrument (non-orthorectified) frame, with orthorectification left to the user via the geolocation lookup table. The L3 stage described here produces the corresponding map-projected representation.

## **3. Algorithm Rationale**

EMIT L1B and L2A products are delivered in the instrument frame with an accompanying geolocation lookup table (GLT). This preserves the native sampling of the instrument: every spectrum in the delivered cube is a measured spectrum, with no resampling applied. Users working spatially, or combining EMIT with other gridded datasets, apply the GLT following EMIT conventions before beginning their analysis.

The L3 stage provides the gridded form directly. The product opens in standard geospatial software with georeferencing already applied.

Datasets are resampled by nearest-neighbor assignment. Each output cell takes the value of the single closest input pixel, so every output spectrum is an unmodified measured spectrum.

The output grid is defined in geographic coordinates (EPSG:4326) at a fixed angular pixel size of 0.00055 degrees, and is snapped to multiples of that pixel size. A global geographic grid gives every granule a consistent and directly comparable definition.

## **4. Algorithm Implementation**

### **4.1 Input Data**

While EMIT data products delivered to the DAAC use DAAC formatting conventions, the L3 gridding stage operates on the internal representation: binary data cubes with detached human-readable ASCII header files. The formatting convention adheres to the ENVI standard, accessible (Jul 2026) at <https://www.nv5geospatialsoftware.com/docs/ENVIHeaderFiles.html>. The header files consist of data fields in equals-sign-separated pairs and describe the layout of the file. The specific input files needed for the L3 stage are:

**1. A surface reflectance file** containing estimated spectral surface reflectance for every pixel, size [rows x cols x channels] in single-precision floating point, in the non-orthorectified instrument frame. The ENVI header supplies `wavelength`, `fwhm`, and, where present, `bbl` (bad band list).

**2. A reflectance uncertainty file** containing predicted per-channel uncertainty as one standard deviation, with the same dimensions and frame as the reflectance file.

**3. An atmospheric state file** containing the spatially smooth atmospheric solution, size [rows x cols x N] in single-precision floating point. Two channels are used, identified by name from the `band names` header field: `AOT550` (aerosol optical depth at 550 nm) and `H2OSTR` (column water vapor, g/cm²). Additional channels may be present and are ignored.

**4. A mask file** containing the L2A mask flags. Two bands are read by index: the `SpecTf-Cloud flag` band to mask the output products, and the `SpecTf-Cloud Probability` band, which is gridded and written to the reflectance product (see Sections 4.2.4 and 4.4).

**4. A location file** containing information about the specific location of each pixel. The channels include:
1. Longitude in decimal degrees, east of zero, WGS-84 datum
2. Latitude in decimal degrees, WGS-84 datum
3. Elevation in meters above mean sea level

**VI. An observation file** containing information about the observation geometry for every pixel. The channels include:
1. Path length - the direct geometric distance from the sensor to the location on the surface of the Earth, as defined by a Digital Elevation Model (DEM)
2. To-sensor azimuth, in decimal degrees, at the surface
3. To-sensor zenith, in decimal degrees, at the surface
4. To-sun azimuth, in decimal degrees, at the surface
5. To-sun zenith, in decimal degrees, at the surface
6. Phase angle in degrees, representing the angular difference between incident and observation rays
7. Terrain slope in degrees as determined from DEMs
8. Terrain aspect in degrees, as determined from DEMs
9. The cosine of the solar incidence angle relative to the surface normal
10. UTC time

"Bad data" at the periphery outside the field of view, or masked as a result of onboard cloud masking or instrument error, is assigned the reserved (floating point) NODATA value -9999 in the input data. In addition to these per-acquisition files, the stage takes a product version string and a software delivery build number, which are recorded in the output global metadata for provenance.

### **4.2 Theoretical Description**

The L3 reflectance algorithm resamples instrument-frame data onto a regular geographic grid. It is a geometric operation only, values are not combined, altered, or recomputed. The steps are grid definition, nearest-neighbor assignment, footprint determination, masking with connected-component filtering, and output writing, with browse image generation as a byproduct of the gridded reflectance.

#### 4.2.1 Grid Definition

Let $\lambda$ and $\phi$ denote the longitude and latitude arrays read from channels 1 and 2 of the location file, each of shape [rows x cols]. For a target angular pixel size $p$ (default $p = 0.00055°$, approximately 61 m at the equator), the grid bounds are snapped outward to whole multiples of $p$:

$$\lambda_{min} = p \left\lfloor \frac{\min(\lambda)}{p} \right\rfloor, \qquad \lambda_{max} = p \left\lceil \frac{\max(\lambda)}{p} \right\rceil$$

$$\phi_{min} = p \left\lfloor \frac{\min(\phi)}{p} \right\rfloor, \qquad \phi_{max} = p \left\lceil \frac{\max(\phi)}{p} \right\rceil$$

Snapping to multiples of $p$ ensures separate granules co-register. Output dimensions follow:

$$n_{lines} = \left\lceil \frac{\phi_{max} - \phi_{min}}{p} \right\rceil, \qquad n_{cols} = \left\lceil \frac{\lambda_{max} - \lambda_{min}}{p} \right\rceil$$

Output cell centers are placed at

$$\lambda_{j} = \lambda_{min} + j\,p, \qquad \phi_{i} = \phi_{max} - i\,p$$

for row index $i$ and column index $j$, giving latitude that decreases with increasing row index in the usual raster convention. The corresponding GDAL geotransform, which references the upper-left corner of the upper-left pixel rather than its center, is

$$\left( \lambda_{min} - \tfrac{p}{2},\; p,\; 0,\; \phi_{max} + \tfrac{p}{2},\; 0,\; -p \right)$$

and is written to the `crs` variable along with the EPSG:4326 well-known text, semi-major axis, and inverse flattening.

#### 4.2.2 Nearest-Neighbor Assignment

The input pixel centers are indexed with a k-d tree, which supports efficient nearest-neighbor lookup over the irregular geolocation of the instrument frame. For each output cell center, the single nearest input pixel is found by Euclidean distance in degree space, and the result is mapped back to an input row and column pair. This index array is computed once and reused for every band of every output file, so reflectance, uncertainty, state variables, mask layers, and observation layers are guaranteed to be mutually consistent: an output pixel draws all of its values from the same input pixel.

#### 4.2.3 Footprint Determination and Fill

A nearest neighbor exists for every output cell, including cells far outside the acquisition. Two independent criteria determine whether a cell is retained.

First, the scene footprint is constructed as a closed polygon by walking the perimeter of the input geolocation array: the first row, then the last column, then the last row reversed, then the first column reversed. This traces the true swath boundary, including its curvature, rather than a bounding box. Output cell centers are tested for containment in this polygon.

Second, a distance criterion admits cells that fall marginally outside the polygon but still lie within half of the input sampling distance of a real pixel, which prevents erosion of a pixel-wide strip along the scene edge. The input sampling distance is estimated as the coarser of the mean spacings along the two axes of the input array.

An output cell is retained if it is inside the footprint polygon or if its nearest-neighbor distance is no greater than half the input sampling distance; otherwise it is set to the fill value -9999. This same fill mask is applied to every band of every output product, including the gridded cloud screening mask of Section 4.2.4.

#### 4.2.4 Cloud Masking

The `SpecTf-Cloud Flag` band of the L2A mask is gridded through the same nearest-neighbor assignment as the science data, using the shared index array of Section 4.2.2, and the resulting gridded mask is then filtered and applied to the output products. Gridding the mask rather than applying it in the instrument frame means the flag and the spectrum it describes are transported together, and the filtering step of Section 4.2.4.1 operates on the same spatial sampling in which the product is delivered.

The screening band is selected by a zero-based index into the mask file (default 5, i.e. the sixth band), and any non-zero value is treated as flagged. Because mask band ordering differs between mask product versions, the index appropriate to the mask version being processed must be supplied; the band list for a given version is given in the EMIT L2A Mask User Guide.

Cells outside the acquisition footprint carry the fill value -9999 in the gridded mask, which is non-zero and therefore flagged. These cells are already fill in every output product, so the outcome is unchanged, but it means the flagged set includes the exterior of the swath as well as the screened pixels within it.

##### 4.2.4.1 Connected-Component Filtering of the Cloud Screening Mask

The gridded cloud screening mask is not applied directly. It is first filtered to drop small flagged regions, which are potential false positives.

First, the mask is labeled into connected components under 4-connectivity, so that two flagged cells belong to the same component only if they share an edge. Unflagged background is excluded from the test that follows. Second, the mask is eroded with a three-by-three structuring element of all ones, so that a cell survives erosion only if it and all eight of its neighbors are flagged. A component is retained in full if any one of its cells survives erosion, and is discarded otherwise.

##### 4.2.4.2 Application to the Output Products

The filtered mask is applied after gridding. Reflectance is gridded band by band and the flagged cells of every band are then set to -9999. The atmospheric state variables and each observation parameter are gridded individually and set to -9999 at the flagged cells in the same way. The identical mask is used in all cases, so a cell that is fill in one layer is fill in all of them, and the spatial extent of valid data is common to the reflectance, state, and observation products.

#### 4.2.5 Cloud Probability Layer

The SpecTf cloud probability band of the L2A mask is gridded and written to the reflectance product as `masks/cloud_probability`, a two-dimensional layer on the output `lat` and `lon` grid. The band is selected by a zero-based index (default 4), and the same nearest-neighbor index array is used, so the probability at an output cell is the probability of the same input pixel that supplied the spectrum at that cell.

The layer is carried through unmodified and allows users whose cloud sensitivity requirements differ from the screening applied here to threshold the probability themselves. Cells outside the footprint carry the fill value -9999. Interpretation of the probability, including its calibration and the meaning of intermediate values, is described in the EMIT L2A Mask ATBD.

#### 4.2.6 Bad Band Handling

The `good_wavelengths` array marks channels where retrieved reflectance is not usable, principally the deep water vapor absorption regions where surface reflectance cannot be recovered. It is taken from the `bbl` field of the input ENVI header when present. For early data in which `bbl` was not populated, it is reconstructed from the wavelength grid by flagging channels in 1325–1435 nm and 1770–1962 nm as unusable and all others as usable.

The flagged channels are retained in the output rather than removed, so that band indices remain stable across granules regardless of which convention produced the flags. The array is advisory: users are expected to subset on it, and no values are altered on its basis.

#### 4.2.7 Browse Image Generation

A quicklook RGB is generated from the gridded and masked reflectance. Bands nearest 660, 550, and 440 nm are selected by minimizing $|\lambda_{band} - \lambda_{target}|$, so that band selection adapts automatically to the wavelength grid of the L1B processing epoch. Fill values are converted to NaN and excluded from the statistics. A single lower and upper bound is taken as the 2nd and 98th percentile computed jointly across all three bands, and applied identically to each. The result is scaled to 8-bit and written as a PNG.

### **4.3 Practical Considerations**

**Pixel size.** The grid is angular, so ground distance is constant north-south and scales with $\cos\phi$ east-west. The default 0.00055° gives approximately 61 m north-south everywhere, and east-west approximately 61 m at the equator, 53 m at 30°, and 38 m at 52°. Output cells are square in degrees, not on the ground.

### **4.4 Output Data**

The specific output files from the L3 Reflectance stage are:

**I. A gridded reflectance file**, containing the `reflectance` variable on dimensions (bands, lat, lon) in single-precision floating point, units of reflectance as a fraction relative to 1. It also carries the `sensor_band_parameters` group (wavelengths, fwhm, good_wavelengths), the `state_variables` group (aerosol_optical_thickness, water_vapor), the `masks` group (cloud_probability, the gridded SpecTf cloud probability on dimensions (lat, lon)), coordinate variables `lon` and `lat` in double precision, a root-level `bands` variable duplicating the wavelength centers for reader convenience, and the `crs` variable.

**II. A gridded reflectance uncertainty file**, containing the `reflectance_uncertainty` variable on dimensions (bands, lat, lon), one standard deviation per channel with covariance ignored, together with the same band parameter, coordinate, and CRS variables.

**III. A gridded observation file**, containing the `observation_parameters` group with eleven 2D layers on dimensions (lat, lon): path length, to-sensor azimuth and zenith, to-sun azimuth and zenith, solar phase, slope, aspect, cosine i, UTC time, and Earth-sun distance.

**IV. A browse image**, a 3-band 8-bit RGB PNG.

| Output file          | Format                                                          | Interpretation                                        |
| -------------------- | --------------------------------------------------------------- | ----------------------------------------------------- |
| Gridded reflectance  | bands x lat x lon, 32-bit float, NetCDF-4, EPSG:4326            | Hemispherical-directional reflectance factor          |
| Gridded uncertainty  | bands x lat x lon, 32-bit float, NetCDF-4, EPSG:4326            | Reflectance uncertainty (one standard deviation)      |
| Gridded observations | lat x lon per layer, 32-bit float, NetCDF-4, EPSG:4326          | Observation geometry and timing                       |
| Browse               | lat x lon x 3, 8-bit unsigned, PNG                              | Stretched RGB quicklook                               |

*Table 1: Output files*

All data variables carry `_FillValue = -9999` and `grid_mapping = "crs"`. Global attributes are inherited from the corresponding input ENVI metadata and augmented with a product title, a summary, and `ncei_template_version = "NCEI_NetCDF_Grid_Template_v2.0"`.

All gridded data variables share the same storage settings. Values are quantized with `least_significant_digit` set to 5, retaining five decimal places to improve compressibility. Variables are chunked at 10 x 256 x 256 in the (bands, lat, lon) case and 256 x 256 in the two-dimensional (lat, lon) case, and are written with zlib compression at complevel 1. The chunk shape spans a block of ten contiguous channels over a 256 by 256 spatial tile, which suits both spectral access at a point and spatial access in a band.

## **5. Calibration, uncertainty characterization and propagation, and validation**

No retrieval is performed at L3, and therefore no new uncertainty is introduced radiometrically. The `reflectance_uncertainty` values are the posterior standard deviations produced by the L2A Optimal Estimation retrieval, carried through the same nearest-neighbor gather as the reflectance itself and referring to the same input pixel. Their characterization and validation are described in the EMIT L2A ATBD, Section 4.

The uncertainty specific to this stage is geometric. It has two components. The first is the geolocation accuracy of the L1B `loc` product, which the gridding inherits, an error in the reported position of an input pixel becomes an error in where its spectrum is placed on the grid. The second is the resampling displacement, the distance between an output cell center and the center of the input pixel assigned to it. For a nearest-neighbor assignment with output sampling comparable to input sampling, the expected displacement is a fraction of a pixel, with a worst case of approximately half the input sampling distance.

## **6. Code Repository and References**

### 6.1 Repository

The L3 gridding code is open source under the Apache 2.0 license and available at:
> <https://github.com/emit-sds/emit-sds-l3rfl>

Supporting NetCDF conversion and metadata utilities are provided by `emit_utils`:
> <https://github.com/emit-sds/emit-utils>
