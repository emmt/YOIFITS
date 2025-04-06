# YOIFITS

This project provides support for OI-FITS (optical interferometry data format) in
[Yorick](http://github.com/LLNL/yorick/).

[Installation](#installation) is explained below.


## Usage

If properly installed, the software exploits the *auto-load* facility of Yorick so that
there are no needs to `#include "oifits.i"`.


### Dealing with an existing OI-FITS file

First, load OI-FITS data from a file:

```c
ws = oifits_load(filename);
```

where `ws` is the *handle* to manipulate the contents of the OI-FITS file and `filename`
is the name of the OI-FITS file.

OI-FITS data is stored in a collection of FITS header data units (HDU). This structure is
reproduced by the plugin where a data-block corresponds to a FITS HDU. In order to access
the OI-FITS contents, you have to select a specific data-block. For instance, you can loop
over all the data-blocks of the handle `ws` with:

```c
for (db = oifits_first(ws); db; db = oifits_next(ws, db)) {
    ...;
}
```

where `db` is another kind of handle but to a specific data-block this time.

In the loop, the function `oifits_is_data` can be used to check whether a given data-block
contains interferometric data (true for `OI_VIS`, `OI_VIS2` and `OI_T3` data-blocks).

There are many functions to access the contents of a given data-block. Again, the OI-FITS
structure of a HDU is reproduced and the general syntax to query a specific field is:

```c
oifits_get_FIELDNAME(ws, db);
```

where `FIELDNAME` is the name of the OI-FITS field which is the same as the corresponding
FITS keyword or column in lowercase letters and with non alphanumerical characters
replaced by an underscore `_`. For instance, to query the effective wavelength, just do:

```c
oifits_get_eff_wave(ws, db);
```

where `db` is an `OI_WAVELENGTH`, `OI_VIS`, `OI_VIS2` or `OI_T3` data-block of OI-FITS
instance `ws`. The reasons that you have to provide both `ws` and `db` are that some
information (as the effective wavelength) of the data-block `db` may be stored in an other
data-block of the same OI-FITS instance `ws` and Yorick does not allow circular
references.

The following is applicable to all data-blocks:
* `oifits_get_revn`: get revision number of the table definition

To query the fields of `OI_VIS`, `OI_VIS2` or `OI_T3` data-blocks:
* `oifits_get_date_obs`: get UTC start date of observations
* `oifits_get_time`: get UTC time of observation (s)
* `oifits_get_mjd`: get Modified Julian Date
* `oifits_get_int_time`: get integration time (s)
* `oifits_get_sta_index`: get station numbers contributing to the data
* `oifits_get_flag`: get flags
* `oifits_get_visamp`: get visibility amplitude
* `oifits_get_visamperr`: get error in visibility amplitude
* `oifits_get_visphi`: get visibility phase (deg)
* `oifits_get_visphierr`: get error in visibility phase (deg)
* `oifits_get_vis2data`: get squared visibility
* `oifits_get_vis2err`: get error in squared visibility
* `oifits_get_t3amp`: get triple-product amplitude
* `oifits_get_t3amperr`: get error in triple product amplitude
* `oifits_get_t3phi`: get triple-product phase (deg)
* `oifits_get_t3phierr`: get error in triple product phase (deg)
* `oifits_get_ucoord`: get u coordinate of the data (m)
* `oifits_get_vcoord`: get v coordinate of the data (m)
* `oifits_get_u1coord`: get u coordinate of baseline AB of the triangle (m)
* `oifits_get_v1coord`: get v coordinate of baseline AB of the triangle (m)
* `oifits_get_u2coord`: get u coordinate of baseline BC of the triangle (m)
* `oifits_get_v2coord`: get v coordinate of baseline BC of the triangle (m)

To query fields of `OI_ARRAY` data block:
* `oifits_get_arrname`: get identifier of corresponding OI_ARRAY
* `oifits_get_frame`: get coordinate frame
* `oifits_get_arrayx`: get X coordinate of array center (m)
* `oifits_get_arrayy`: get Y coordinate of array center (m)
* `oifits_get_arrayz`: get Z coordinate of array center (m)
* `oifits_get_tel_name`: get telescope name
* `oifits_get_sta_name`: get station name
* `oifits_get_sta_index`: get station number
* `oifits_get_diameter`: get element diameter (m)
* `oifits_get_staxyz`: get station coordinate relative to array center (m)

To query fields of `OI_TARGET` data block:
* `oifits_get_target_id`: get index number of target(s)
* `oifits_get_target`: get target names(s)
* `oifits_get_raep0`: get R.A. at mean equinox (deg)
* `oifits_get_decep0`: get decl. at mean equinox (deg)
* `oifits_get_equinox`: get equinox
* `oifits_get_ra_err`: get error in R.A. at mean equinox (deg)
* `oifits_get_dec_err`: get error in decl. at mean equinox (deg)
* `oifits_get_sysvel`: get Systemic radial velocity (m/s)
* `oifits_get_veltyp`: get reference for radial velocity ("LSR", "GEOCENTR",
  etc.)
* `oifits_get_veldef`: get definition of radial velocity ("OPTICAL", "RADIO")
* `oifits_get_pmra`: get proper motion in R.A. (deg/yr)
* `oifits_get_pmdec`: get proper motion in decl. (deg/yr)
* `oifits_get_pmra_err`: get error of proper motion in R.A. (deg/yr)
* `oifits_get_pmdec_err`: get error of proper motion in decl. (deg/yr)
* `oifits_get_parallax`: get parallax (deg)
* `oifits_get_para_err`: get error in parallax (deg)
* `oifits_get_spectyp`: get spectral type

To query fields of `OI_WAVELENGTH` data block:
* `oifits_get_insname`: get identifier of corresponding OI_WAVELENGTH
* `oifits_get_eff_wave`: get effective wavelength of channel (m)
* `oifits_get_eff_band`: get effective bandpass of channel (m)


### Creating a new OI-FITS file

In order to create a new OI-FITS file, you first create a new OI-FITS handle in Yorick,
then you populate it with data-blocks and, finally, you save it to the disk.

To create a new OI-FITS instance:

```c
ws = oifits_new();
```

To add data-blocks:

```c
oifits_insert, ws, db1, db2, ...;
```

where `db1`, `db2`, *etc.* are OI-FITS data-blocks which have been freshly created (see
below) or which are borrowed from another OI-FITS instance.

To create a new data-block, the general syntax is:

```c
db = oifits_new_DBTYPE(key1=val1, key2=val2, ...);
```

where `DBTYPE` is the data-block type (`target`, `array`, `wavelength`, `spectrum`, `vis`,
`vis2` or `t3`) and all fields of the data-block are passed by keyword. Note that all
fields must be specified. See the individual documentation of the data-block constructors
to figure out which fields are required. For instance:

```c
help, oifits_new_vis2;
```

Optionally, the OI-FITS instance to which insert the new data-block can be specified in a
data-block constructor:

```c
oifits_new_DBTYPE, ws, key1=val1, key2=val2, ...;
```

is the same as:

```c
oifits_insert, ws, oifits_new_DBTYPE(key1=val1, key2=val2, ...);
```

There may be any number of data-blocks in an OI-FITS instance and they may be inserted at
any time and in any order. After inserting the data-blocks, it is necessary to make sure
that the internals of the OI-FITS instance are consistent (otherwise some functionalities
may not work as expected). Updating internal information is done by:

```c
oifits_update, ws;
```

Saving an OI-FITS instance to a file is as simple as:

```c
oifits_save, ws, filename;
```

where `filename` is the name of the OI-FITS file. The `oifits_save` subroutine has a
number of keywords to add comments, history records or to allow overwriting an existing
file (which is forbidden by default).


### Simulated Data

If the contents of an OI-FITS data-block correspond to model data, you may simulate
additive noise with:

```c
out = oifits_add_noise(ws, method, level);
```

which has the effect of adding noise to all `OI_VIS`, `OI_VIS2` or `OI_T3` data-blocks of
the OI-FITS instance `ws` and returns a new OI-FITS instance. The other arguments are as
follows:

* `method = 1` or `"generate"` to add noise to noiseless data; the standard deviation of
  the noise is taken from the contents of `ws`; in this case, the `level` argument must be
  undefined or omitted.

* `method = 2` or `"snr"` to add noise to noiseless data; the standard deviation of noise
  is computed to achieve a signal-to-noise ratio equal to the value of `level`.

* `method = 3` or `"amplify"` to add noise to noisy data so that the standard deviation of
  total noise (existing one plus added one) is multiplied by the value of `level` which
  must be greater or equal one. The standard deviation of the noise prior to the
  amplification is taken from the contents of `ws`.

The operation can also be carried out *in-place* by calling `oifits_add_noise` as a
subroutine:

```c
oifits_add_noise, ws, method, level;
```


## Installation

### Installation with EasyYorick

Installation of YOIFITS by [EasyYorick](https://github.com/emmt/EasyYorick) is fully
supported. Assuming you have installed EasyYorick, you just have to execute:

```sh
ypkg install yorick yeti yoifits
```

which should install Yorick, Yeti (if not yet installed) and YOIFITS.

To upgrade to the last master version:

```sh
ypkg upgrade yoifits
```


### Manual installation

1. You must have [Yorick](http://github.com/LLNL/yorick/) and
   [Yeti](https://github.com/emmt/Yeti) installed on your machine.

2. Unpack the [software code](https://github.com/emmt/YOIFITS/archive/master.zip)
   somewhere or clone the Git repository:

   ```sh
   git clone https://github.com/emmt/YOIFITS.git yoifits
   ```

   if you want/prefer to use HTTPS, or:

   ```sh
   git clone git@github.com:emmt/YOIFITS.git yoifits
   ```

   if you want/prefer to use SSH. Any of these commands creates a local GIT repository
   named `yoifits`.


3. Configure for compilation. There are two possibilities (the first one is recommended):

   - For an **out-of-place build**, create a dedicated build directory, say `$BUILD_DIR`,
     go to the build directory and run the configuration script:

     ```sh
     mkdir -p $BUILD_DIR
     cd $BUILD_DIR
     $SRC_DIR/configure
     ```

     where `$SRC_DIR` is the path to the source directory of the plug-in code. To see the
     configuration options, type:

     ```sh
     $SRC_DIR/configure --help
     ```

   - For an **in-place build**, go to the source directory, say `$SRC_DIR`, of the plug-in
     code and run the configuration script:

     ```sh
     cd $SRC_DIR
     ./configure
     ```

     To see the configuration options, type:

     ```sh
     ./configure --help
     ```

4. Compile the code:

   ```sh
   make clean
   make
   ```

5. Install YOIFITS (you must have write access granted to Yorick directories):

   ```sh
   make install
   ```
