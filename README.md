# YOIFITS

This project provides support for OI-FITS (optical interferometry data format) in
[Yorick](http://yorick.github.com/).

To use this software you must have Yorick (of course) and the
[Yeti](https://github.com/emmt/Yeti) plugin installed.

## Installation

Copy `oifits.i` in `${Y_SITE}/i/` or in your Yorick directory.
Then from Yorick:
```
    #include "oifits.i"
```


## Usage

### Dealing with an existing OI-FITS file

First, load OI-FITS data from a file:
```
    ws = oifits_load(filename);
```
where `ws` is the *handle* to manipulate the contents of the OI-FITS file and
`filename` is the name of the OI-FITS file.

OI-FITS data is stored in a collection of FITS header data units (HDU).  This
structure is reproduced by the plugin where a datablock correspond to a FITS
HDU.  In order to access the OI-FITS contents, you have to select a specific
datablock.  For instance, you can loop over all the datablocks of the handle
`ws` with:
```
    for (db = oifits_first(ws); db; db = oifits_next(ws, db)) {
        ...;
    }
```
where `db` is another kind of handle but to a specific datablock this time.

There are many functions to access the contents of a given datablock.  Again,
the OI-FITS structure of a HDU is reproduced and the general syntax to query a
specific field is:
```
    oifits_get_FIELDNAME(ws, db);
```
where `FIELDNAME` is the name of the FITS keyword or OI-FITS column.  For
instance, to query the effective wavelength of an **OI_WAVELENGTH**,
**OI_VIS**, **OI_VIS2** or **OI_T3** datablock `db`, just do:
```
    oifits_get_eff_wave(ws, db);
```
The reasons that you have to provide both `ws` and `db` are that some
information (as the effective wavelength) of the datablock `db` may be stored
in an other datablock of the same OI-FITS handle `ws` and Yorick does not allow
circular references.

* `oifits_get_revn`      - get revision number of the table definition

To query the contents of **OI_VIS**, **OI_VIS2** or **OI_T3** datablocks:
* `oifits_get_date_obs`  - get UTC start date of observations
* `oifits_get_time`      - get UTC time of observation (s)
* `oifits_get_mjd`       - get Modified Julian Date
* `oifits_get_int_time`  - get integration time (s)
* `oifits_get_sta_index` - get station numbers contributing to the data
* `oifits_get_flag`      - get flags
* `oifits_get_visamp`    - get visibility amplitude
* `oifits_get_visamperr` - get error in visibility amplitude
* `oifits_get_visphi`    - get visibility phase (deg)
* `oifits_get_visphierr` - get error in visibility phase (deg)
* `oifits_get_vis2data`  - get squared visibility
* `oifits_get_vis2err`   - get error in squared visibility
* `oifits_get_t3amp`     - get triple-product amplitude
* `oifits_get_t3amperr`  - get error in triple product amplitude
* `oifits_get_t3phi`     - get triple-product phase (deg)
* `oifits_get_t3phierr`  - get error in triple product phase (deg)
* `oifits_get_ucoord`    - get u coordinate of the data (m)
* `oifits_get_vcoord`    - get v coordinate of the data (m)
* `oifits_get_u1coord`   - get u coordinate of baseline AB of the triangle (m)
* `oifits_get_v1coord`   - get v coordinate of baseline AB of the triangle (m)
* `oifits_get_u2coord`   - get u coordinate of baseline BC of the triangle (m)
* `oifits_get_v2coord`   - get v coordinate of baseline BC of the triangle (m)

To query attributes for **OI_ARRAY** data block:
* `oifits_get_arrname`   - get identifier of corresponding OI_ARRAY
* `oifits_get_frame`     - get coordinate frame
* `oifits_get_arrayx`    - get X coordinate of array center (m)
* `oifits_get_arrayy`    - get Y coordinate of array center (m)
* `oifits_get_arrayz`    - get Z coordinate of array center (m)
* `oifits_get_tel_name`  - get telescope name
* `oifits_get_sta_name`  - get station name
* `oifits_get_sta_index` - get station number
* `oifits_get_diameter`  - get element diameter (m)
* `oifits_get_staxyz`    - get station coordinate relative to array center (m)

To query attributes for **OI_TARGET** data block:
* `oifits_get_target_id` - get index number of target(s)
* `oifits_get_target`    - get target names(s)
* `oifits_get_raep0`     - get R.A. at mean equinox (deg)
* `oifits_get_decep0`    - get decl. at mean equinox (deg)
* `oifits_get_equinox`   - get equinox
* `oifits_get_ra_err`    - get error in R.A. at mean equinox (deg)
* `oifits_get_dec_err`   - get error in decl. at mean equinox (deg)
* `oifits_get_sysvel`    - get Systemic radial velocity (m/s)
* `oifits_get_veltyp`    - get reference for radial velocity ("LSR", "GEOCENTR", etc.)
* `oifits_get_veldef`    - get definition of radial velocity ("OPTICAL", "RADIO")
* `oifits_get_pmra`      - get proper motion in R.A. (deg/yr)
* `oifits_get_pmdec`     - get proper motion in decl. (deg/yr)
* `oifits_get_pmra_err`  - get error of proper motion in R.A. (deg/yr)
* `oifits_get_pmdec_err` - get error of proper motion in decl. (deg/yr)
* `oifits_get_parallax`  - get parallax (deg)
* `oifits_get_para_err`  - get error in parallax (deg)
* `oifits_get_spectyp`   - get spectral type

To query attributes for **OI_WAVELENGTH** data block:
* `oifits_get_insname`   - get identifier of corresponding OI_WAVELENGTH
* `oifits_get_eff_wave`  - get effective wavelength of channel (m)
* `oifits_get_eff_band`  - get effective bandpass of channel (m)
