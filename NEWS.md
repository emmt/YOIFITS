# Changes in `YOIFITS`

This file describes the most important changes in `YOIFITS`, a Yorick interface to the
OI-FITS data format.

This file structure is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec).

## Version 1.2.0 (2025-04-07)

### Added

- Function `oifits_save` can create a new OI-FITS file if destination is specified as a
  string (as before) or can append to an existing FITS file if destination is specified as
  a FITS handle.

## Version 1.1.0 (2025-04-06)

### Changed

None.

### Added

- Optional columns `NS_MODEL_$x` in OIFITS data-blocks `OI_VIS`, `OI_VIS2`, and `OI_T3`
  with model of the data in column `$x` can be loaded/saved from/to OIFITS files. In the
  object storing OIFITS data, these columns correspond to members with the same name in
  lower case letters with the `NS_` prefix suppressed.

- New routines: `oifits_select_target`, `oifits_warn`, `oifits_list_targets`, and
  `oifits_list_instruments`.

- Automatically fix missing `OI_REVN` in `OI_FLUX` for Gravity data.

- Optionally tolerate a small discrepancy of values when comparing tables.

### Fixed

- Fix getting units in OI-FLUX.

- Simplify and improve `oifits_merge`
