# Changes in `YOIFITS`

This file describes the most important changes in `YOIFITS`, a Yorick interface to the
OI-FITS data format.

This file structure is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec).

## Unreleased

### Changed

None.

### Added

- New routines: `oifits_select_target`, `oifits_warn`, `oifits_list_targets`, and
  `oifits_list_instruments`.

- Automatically fix missing `OI_REVN` in `OI_FLUX` for Gravity data.

- Optionally tolerate a small discrepancy of values when comparing tables.

### Fixed

- Fix getting units in OI-FLUX.

- Simplify and improve `oifits_merge`
