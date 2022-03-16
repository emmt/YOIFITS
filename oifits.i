/*
 * oifits.i -
 *
 * Reading and writing of optical/IR interferometric data files in OI-FITS
 * format for Yorick + Yeti.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2005 Éric Thiébaut, Clémentine Béchet, Julien Salmon.
 * Copyright (C) 2006-2018 Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *-----------------------------------------------------------------------------
 */

/* OI-FITS is based on fits.i */
if (! is_func(fits_open)) include, "fits.i", 1;

/* OI-FITS requires hash tables provided by Yeti: */
if (! is_func(h_new)) include, "yeti.i", 1;

/*
 * IMPORTANT
 * =========
 *
 * OI-FITS public interface only consists in functions (with names starting
 * with "oifits_") and opaque objects.  The caller should never directly
 * access to internals of theses opaque objects otherwise this could give
 * raise to conflicts with the future evolution of the OI-FITS interface.
 *
 *
 * IMPLEMENTATION NOTES
 * ====================
 *
 * A set of OI-FITS data is memorized as a hash table, the 'master', with
 * children hash tables, one for each data block (e.g. OI_TARGET,
 * OI_WAVELENGTH, OI_VIS, ...).
 *
 * Structure of the 'master' hash table:
 *   master.__revn         [integer] revision number
 *   master.__counter      [integer] data block counter
 *   master.__array        [string] names of OI_ARRAY data blocks
 *   master.__corr         [string] names of OI_CORR data blocks
 *   master.__data         [string] names of OI_VIS, OI_VIS2 and OI_T3 data
 *                                  blocks
 *   master.__errmsg       [string] pending error message(s)
 *   master.__first        [string] name of first data block
 *   master.__db#          [hash]   data blocks (# is an integer)
 *   master.__target       [string] name of OI_TARGET data block
 *   master.__wavelength   [string] names of OI_WAVELENGTH data blocks
 *   master.__clean        [logical] true if update of internals not needed
 *
 * Structure of a data block:
 *   db.__nbaselines   [integer] number of observed baselines
 *   db.__nwavelengths [integer] number of spectral channels
 *   db.__is_data  [integer] this data block contains measurements
 *   db.__ins      [string]  name ("__db#') of related OI_WAVELENGTH in master
 *   db.__arr      [string]  name ("__db#') of related OI_ARRAY in master
 *   db.__corr     [string]  name ("__db#') of related OI_CORR in master
 *   db.__self     [string]  name ("__db#') of this data block in its parent
 *   db.__index    [integer] index of this data block
 *   db.__type     [integer] data block type (see oifits_get_type)
 *   db.__next     [string]  name of next data block in same parent
 *   db.__class    [string]  data block class: "TARGET", or "VIS2", ...
 *   db.KEY        column/keyword value (except OI_REVN -> revn)
 *   db.KEY_units  units for KEY
 *
 * FITS BINTABLE are read as array of pointers:
 *
 *   <-------- NCOLS ---------------->
 *  ^
 *  |       |              |
 *  |       |<-- NCELLS -->|
 * NROWS    |              |
 *  |       |              |
 *  |       |              |
 *  |       |              |
 *  V       |              |
 *
 *
 * REFERENCES
 * ==========
 * [1] Pauls, T. A.; Young, J. S.; Cotton, W. D. & Monnier, J. D.: "A Data
 *     Exchange Standard for Optical (Visible/IR) Interferometry",
 *     Publications of the Astronomical Society of Pacific vol. 117,
 *     pp. 1255-1262 (2005).
 */

/*---------------------------------------------------------------------------*/

/* Some constants. */
local OIFITS_PI, OIFITS_MICRON, OIFITS_DEGREE, OIFITS_ARCSECOND;
/* DOCUMENT OIFITS_PI             = 3.1415.....
         or OIFITS_MICRON         = micron to meter conversion factor
         or OIFITS_DEGREE         = degree to radian conversion factor
         or OIFITS_ARCSECOND      = arcsecond to radian conversion factor

   SEE ALSO: oifits_new.
 */
OIFITS_PI = 3.141592653589793238462643383279503;
OIFITS_DEGREE = OIFITS_PI/180.0;
OIFITS_ARCSECOND = OIFITS_DEGREE/3600.0;
OIFITS_MICRON = 1e-6;
OIFITS_NANOMETER = 1e-9;

func oifits_new(nil)
/* DOCUMENT oifits_new()
     Creates a new OI-FITS master object.

   SEE ALSO: oifits_load, oifits_save, oifits_insert.
 */
{
  if (! is_void(nil)) error, "expecting nil argument";
  return h_new(__counter = 0, __clean = 1n, __first = "");
}

/* steps:
 *   1. figure out NBASELINES and NWAVELENGTHS
 *   2. create default (zero) structure
 *   3. fill with provided data
 *   4. complete missing data (e.g., ucoord, vcoord)
 */

func oifits_insert(master, ..)
/* DOCUMENT oifits_insert, master, db1, db2, ...;
         or oifits_insert(master, db1, db2, ...);

     Inserts new OI-FITS data-blocks DB1, DB2 etc into OI-FITS instance
     MASTER.  DBn are data-blocks as returned by one of the oifits_new_*
     constructors (which see).  DBn can also be data-blocks from another
     OI-FITS instance, in this case the contents of the source data-block
     is cloned into a new data-block since OI-FITS instances cannot share
     data.  When called as a function, the returned value is MASTER.


   EXAMPLE:
     Create a new OI-FITS instance and save it into a file:
      master = oifits_new();
      target_db = oifits_new_target(target = ...);
      wavelength_db = oifits_new_wavelength(eff_wave = ...);
      oifits_insert, master, db1;
      oifits_insert, master, db2;
      ...;
      oifits_update, master;
      oifits_save, master, filename;

   SEE ALSO: oifits_get_type, oifits_load, oifits_save, oifits_new,
             oifits_new_target, oifits_new_array, oifits_new_wavelength,
             oifits_new_vis, oifits_new_vis2, oifits_new_t3. */
{
  local db;

  counter = master.__counter;
  if (is_void(counter)) {
    counter = 0L;
  }

  while (more_args()) {
    eq_nocopy, db, next_arg();
    if (! is_hash(db) || ! h_has(db, "__class") || ! h_has(db, "__type")
        || ! h_has(db, "__is_data") || ! h_has(db, "__nbaselines")
        || ! h_has(db, "__nwavelengths")) {
      error, ("invalid input data-block"
              + " (use one of the oifits_new_* constructors)");
    }
    db = oifits_unlink_datablock(db);

    /* Get name for new data block and insert in master. */
    do {
      self = swrite(format="__db%d", ++counter);
    } while (h_has(master, self));
    h_set, db, __next = master.__first, __self = self, __index=counter;
    h_set, master, __counter = counter, self, db, __first = self, __clean = 0n;
  }

  return master;
}

func oifits_unlink_datablock(db)
/* DOCUMENT idb = oifits_unlink_datablock(db);
     Make sure that OI-FITS data-block DB is not owned.  DB itself is returned
     if it is not owned, otherwise a fresh copy of DB contents is returned.

   SEE ALSO: oifits_new_target.
 */
{
  if (h_has(db, "__self")) {
    /* Input data-block is owned, make a private copy. */
    cpy = h_new();
    _oifits_copy_member, cpy, db, "__type";
    _oifits_copy_member, cpy, db, "__class";
    _oifits_copy_member, cpy, db, "__is_data";
    _oifits_copy_member, cpy, db, "__nbaselines";
    _oifits_copy_member, cpy, db, "__nwavelengths";
    for (key = h_first(db); key; key = h_next(db, key)) {
      if (strpart(key, 1:2) != "__") {
        _oifits_copy_member, cpy, db, key;
      }
    }
    return cpy;
  } else {
    /* Input data-block is not owned, return a new reference to it. */
    return db;
  }
}

func oifits_update(master, errmode=, revn=, force=, quiet=)
/* DOCUMENT oifits_update(master)

     Update internals of OI-FITS main object MASTER and return it.

     If keyword ERRMODE is true, any inconsistency will be printed out and an
     error will be raised.  The default is ERRMODE=1.

     If keyword QUIET is true, warning messages are suppressed.

     Keyword REVN can be set to the value of the new revision number of
     OI-FITS format.

     Keyword FORCE can be set true to force recomputation of internal cache
     and links even if it is not necessary.

     This function _must_ be called after any modification of MASTER
     (e.g. oifits_insert).  The finalization consists in:

      - (re)building chained list of data-blocks (master.__first,
        db.__next, ...)
      - (re)building lists of OI_ARRAY, OI_TARGET, ..., data-blocks
      - linking data to related wavelength and array blocks
      - supplying default FLAG
      - (NOT YET DONE) compute missing u-v coordinates

   SEE ALSO: oifits_new, oifits_insert, oifits_load.
 */
{
  /* Fast shortcut. */
  if (master.__clean && ! force) {
    return master;
  }

  /* Transfer stack of internal error messages. */
  extern _oifits_error_stack;
  local _oifits_error_count;
  _oifits_error_count = 0;
  _oifits_error_stack = h_pop(master, "__errmsg");
  if (is_void(errmode)) errmode = 1n;

  /* Extract and re-order list of OI-FITS datablocks, then rebuild the list
     of all datablocks. */
  key_list = h_keys(master);
  if (is_void(key_list)) return master;
  key_list = key_list(where(strgrep("^__db[0-9]+$", key_list)(2,) >= 0));
  if ((ndb = numberof(key_list)) < 1) return master;
  index_list = array(long, ndb);
  if (sread(key_list, format="__db%d", index_list) != ndb) {
    error, "bug or corrupted opaque OI-FITS object";
  }
  i = sort(index_list);
  index_list = index_list(i);
  key_list = key_list(i);

  /* First pass to get revision numbers and build chained list of
     datablocks. */
  counter = 0;
  prev = [];
  revn_list = array(long, ndb);
  for (i = 1; i <= ndb; ++i) {
    key = key_list(i);
    db = h_get(master, key);
    if (! db) {
      /* should never happens, however... */
      h_pop, master, key;
      continue;
    }
    if (index_list(i) != ++counter) {
      /* rename data-block */
      newkey = swrite(format="__db%d", counter);
      if (h_has(master, newkey)) error, "this should not happens";
      h_set, master, newkey, h_pop(master, key);
      h_set, db, __index = counter, __self = newkey;
      eq_nocopy, key, newkey;
    }
    if (is_hash(prev)) {
      h_set, prev, __next=key;
    } else {
      h_set, master, __first=key;
    }
    h_pop, db, "__next";
    prev = db;
    /* store revision number */
    r = db.revn;
    if (! is_integer(r) || ! is_scalar(r) || r < 1 || r > _OIFITS_REVN_MAX) {
      if (_oifits_error("bad/missing OI_REVN for data block %d", db.__index)) {
        _oifits_report_error;
      }
    } else {
      revn_list(i) = r;
    }
  }

  /* Check new revision number. */
  if (is_void(revn)) {
    revn = max(revn_list);
  } else if (! is_integer(revn) || ! is_scalar(revn) || revn < 1) {
    error, "bad value for REVN keyword";
  } else if (revn > _OIFITS_REVN_MAX) {
    error, swrite(format="unsupported version of OI-FITS format (REVN = %d)",
                  revn);
  } else if (revn < max(revn_list)) {
    error, "decreasing revision number is not allowed";
  } else {
    revn = long(revn);
  }

  /* Build other links and fix datablocks. */
  data_list = [];
  array_list = [];
  target_list = [];
  corr_list = [];
  wavelength_list = [];
  ins_table = h_new();
  arr_table = h_new();
  corr_table = h_new();
  for (key = master.__first; key; key = db.__next) {
    db = master(key);

    /* Get OI-FITS type and instanciate class-specific datablock members. */
    type = db.__type;
    if (is_integer(type) && is_scalar(type) &&
        1 <= type && type <= numberof(_OIFITS_DATABLOCK_CLASS)) {
      type = long(type);
      class = _OIFITS_DATABLOCK_CLASS(type);
    } else {
      if (_oifits_error("bad OI-FITS type for data block %d", db.__index)) {
        _oifits_report_error;
      }
      continue;
    }
    is_data = (type == OIFITS_TYPE_VIS ||
               type == OIFITS_TYPE_VIS2 ||
               type == OIFITS_TYPE_T3 ||
               type == OIFITS_TYPE_FLUX);
    h_set, db, __self=key, __type=type, __class=class, __is_data=is_data;

    /* Check INSNAME and build table of OI_WAVELENGTH datablocks. */
    if (is_data || type == OIFITS_TYPE_WAVELENGTH) {
      insname = db.insname;
      if (! is_string(insname) || ! is_scalar(insname)) {
        if (_oifits_error("bad/missing INSNAME for data block %d",
                          db.__index)) {
          _oifits_report_error;
        }
        continue;
      }
      insname = oifits_fix_name(insname);
      h_set, db, insname=insname;
      if (type == OIFITS_TYPE_WAVELENGTH) {
        if (h_has(ins_table, insname)) {
          if (_oifits_error("too many OI_WAVELENGTH blocks with INSNAME = '%s'",
                            insname)) {
            _oifits_report_error;
          }
        } else {
          h_set, ins_table, insname, key;
        }
      }
    }

    /* Check ARRNAME and build table of OI_ARRAY datablocks. */
    if (is_data || type == OIFITS_TYPE_ARRAY) {
      arrname = db.arrname;
      if (is_string(arrname) && is_scalar(arrname)) {
        /* remove trailing spaces and register ARRNAME */
        arrname = oifits_fix_name(arrname);
        h_set, db, arrname=arrname;
      } else if (type == OIFITS_TYPE_ARRAY || ! is_void(arrname)) {
        if (_oifits_error("bad/missing ARRNAME for data block %d",
                          db.__index)) {
          _oifits_report_error;
        }
        continue;
      }
      if (type == OIFITS_TYPE_ARRAY) {
        if (h_has(arr_table, arrname)) {
          if (_oifits_error("too many OI_ARRAY blocks with ARRNAME = '%s'",
                            arrname)) {
            _oifits_report_error;
          }
        } else {
          h_set, arr_table, arrname, key;
        }
      }
    }

    /* Check CORRNAME and build table of OI_CORR datablocks. */
    if (is_data || type == OIFITS_TYPE_CORR) {
      corrname = db.corrname;
      if (is_string(corrname) && is_scalar(corrname)) {
        /* remove trailing spaces and register CORRNAME */
        corrname = oifits_fix_name(corrname);
        h_set, db, corrname=corrname;
      } else if (type == OIFITS_TYPE_CORR || ! is_void(corrname)) {
        if (_oifits_error("bad/missing CORRNAME for data block %d",
                          db.__index)) {
          _oifits_report_error;
        }
        continue;
      }
      if (type == OIFITS_TYPE_CORR) {
        if (h_has(corr_table, corrname)) {
          if (_oifits_error("too many OI_CORR blocks with CORRNAME = '%s'",
                            corrname)) {
            _oifits_report_error;
          }
        } else {
          h_set, corr_table, corrname, key;
        }
      }
    }

    /* Grow list of target/data. */
    if (type == OIFITS_TYPE_TARGET) {
      grow, target_list, key;
    } else if (type == OIFITS_TYPE_ARRAY) {
      grow, array_list, key;
    } else if (type == OIFITS_TYPE_WAVELENGTH) {
      grow, wavelength_list, key;
    } else if (type == OIFITS_TYPE_CORR) {
      grow, corr_list, key;
    } else {
      grow, data_list, key;
    }
  }

  /* Deal with having more than one TARGET datablock. */
  n = numberof(target_list);
  if (n == 0) {
    if (_oifits_error("missing mandatory OI_TARGET datablock")) {
      _oifits_report_error;
    }
  } else if (n > 1) {
    if (_oifits_error("too many OI_TARGET datablocks")) {
      _oifits_report_error;
    }
  }

  /* Fix internal links. */
  nil = [];
  for (key = master.__first; key; key = db.__next) {
    db = master(key);
    insname = db.insname;
    ins_lnk = (insname ? ins_table(db.insname) : nil);
    arrname = db.arrname;
    arr_lnk = (arrname ? arr_table(db.arrname) : nil);
    corrname = db.corrname;
    corr_lnk = (corrname ? corr_table(db.corrname) : nil);
    h_set, db, __ins=ins_lnk, __arr=arr_lnk, __corr=corr_lnk;
    if (is_void(ins_lnk) && (db.__is_data ||
                             db.__type == OIFITS_TYPE_WAVELENGTH) &&
        _oifits_error("bad/missing INSNAME for data block %d", db.__index)) {
      _oifits_report_error;
    }
    if (db.__is_data && ! is_void(ins_lnk) &&
        master(ins_lnk).__nwavelengths != db.__nwavelengths) {
      if (_oifits_error("inconsistent number of wavelengths for data block %d",
                        db.__index)) {
        _oifits_report_error;
      }
      h_pop, db, "__ins";
    }
    if (!quiet && is_void(arr_lnk) &&
        (db.__is_data || db.__type == OIFITS_TYPE_ARRAY)) {
      oifits_warn, swrite(format="no ARRNAME for data block %d", db.__index);
    }
  }

  /* Instanciate members of MASTER. */
  h_set, master,
    __revn = revn,
    __array = array_list,
    __wavelength = wavelength_list,
    __data = data_list,
    __target = target_list,
    __counter = counter,
    __clean = 1n;
  if (_oifits_error_count) {
    h_set, master, "__errmsg", _oifits_error_stack;
    if (errmode) {
      write, format="*** ERROR: %s\n", _oifits_error_stack;
      error, "too many errors";
    }
  }
  return master;
}

func oifits_merge(.., quiet=, dest=, atol=, rtol=)
/* DOCUMENT data = oifits_merge(...);
     Merge OI-FITS files (given their names) or OI-DATA.  For now, only data
     with the same OI_TARGET and OI_ARRAY data blocks can be merged.

     Keyword DEST can be specified with an existing OI-FITS instance into which
     merge all inputs.  By default an empty destination is created.

     Keywords ATOL and RTOL can be used to specify the absolute and relative
     tolerances when comparing numerical values.

   SEE ALSO: oifits_new, oifits_load, oifits_save, oifits_is_approx.
 */
{
  target_db = [];
  arr_table = h_new();
  ins_table = h_new();
  corr_table = h_new();
  inspol_table = h_new();
  if (is_void(dest)) {
    dest = oifits_new();
  } else {
    /* Register pre-existing entries that may be linked by data-blocks. */
    for (db = oifits_first(dest); db; db = oifits_next(dest, db)) {
      type = oifits_get_type(db);
      if (type == OIFITS_TYPE_TARGET) {
        if (is_void(target_db)) {
          target_db = db;
        } else {
          error, "more than one OI_TARGET in destination";
        }
      } else if (type == OIFITS_TYPE_ARRAY) {
        name = oifits_fix_name(oifits_get_arrname(dest, db));
        if (h_has(arr_table, name)) {
          error, ("more than one OI_ARRAY named \"" + name
                  + "\" in destination");
        } else {
          h_set, arr_table, name, db;
        }
      } else if (type == OIFITS_TYPE_INSPOL) {
        name = oifits_fix_name(oifits_get_arrname(dest, db));
        if (h_has(inspol_table, name)) {
          error, ("more than one OI_INSPOL named \"" + name
                  + "\" in destination");
        } else {
          h_set, inspol_table, name, db;
        }
      } else if (type == OIFITS_TYPE_WAVELENGTH) {
        name = oifits_fix_name(oifits_get_insname(dest, db));
        if (h_has(ins_table, name)) {
          error, ("more than one OI_WAVELENGTH named \"" + name
                  + "\" in destination");
        } else {
          h_set, ins_table, name, db;
        }
      } else if (type == OIFITS_TYPE_CORR) {
        name = oifits_fix_name(oifits_get_corrname(dest, db));
        if (h_has(ins_table, name)) {
          error, ("more than one OI_CORR named \"" + name
                  + "\" in destination");
        } else {
          h_set, corr_table, name, db;
        }
      }
    }
  }

  nargs = 0;
  name_list = [];
  name_index = 1;
  while (1n) {
    /* Peek next argument of next file name in list. */
    if (name_index < numberof(name_list)) {
      filename = name_list(++name_index);
      src = oifits_load(filename, quiet=quiet);
    } else if (more_args()) {
      local src;
      eq_nocopy, src, next_arg();
      ++nargs;
      if (is_string(src)) {
        eq_nocopy, name_list, src;
        name_index = 0;
        continue;
      }
      filename = swrite(format="ARG%d", nargs);
    } else {
      return oifits_update(dest, quiet=quiet);
    }
    for (db = oifits_first(src); db; db = oifits_next(src, db)) {
      type = oifits_get_type(db);
      if (type == OIFITS_TYPE_TARGET) {
        if (is_void(target_db)) {
          target_db = oifits_unlink_datablock(db);
          oifits_insert, dest, target_db;
        } else {
          msg = _oifits_compare_tables(target_db, db, atol, rtol);
          if (msg) {
            error, swrite(format="cannot merge OI_TARGET from \"%s\" (%s)", filename, msg);
          }
        }
      } else if (type == OIFITS_TYPE_ARRAY) {
        name = oifits_fix_name(oifits_get_arrname(src, db));
        if (h_has(arr_table, name)) {
          msg = _oifits_compare_tables(h_get(arr_table, name), db, atol, rtol);
          if (msg) {
            error, swrite(format="cannot merge OI_ARRAY \"%s\" from \"%s\" (%s)",
                          name, filename, msg);
          }
        } else {
          new_db = oifits_unlink_datablock(db);
          h_set, new_db, arrname = name;
          h_set, arr_table, name, new_db;
          oifits_insert, dest, new_db;
        }
      } else if (type == OIFITS_TYPE_INSPOL) {
        name = oifits_fix_name(oifits_get_arrname(src, db));
        if (h_has(inspol_table, name)) {
          msg = _oifits_compare_tables(h_get(inspol_table, name), db,
                                       atol, rtol);
          if (msg) {
            error, swrite(format="cannot merge OI_INSPOL \"%s\" from \"%s\" (%s)",
                          name, filename, msg);
          }
        } else {
          new_db = oifits_unlink_datablock(db);
          h_set, new_db, arrname = name;
          h_set, inspol_table, name, new_db;
          oifits_insert, dest, new_db;
        }
      } else if (type == OIFITS_TYPE_WAVELENGTH) {
        name = oifits_fix_name(oifits_get_insname(src, db));
        if (h_has(ins_table, name)) {
          msg = _oifits_compare_tables(h_get(ins_table, name), db, atol, rtol);
          if (msg) {
            error, swrite(format="cannot merge OI_WAVELENGTH \"%s\" from \"%s\" (%s)",
                          name, filename, msg);
          }
        } else {
          new_db = oifits_unlink_datablock(db);
          h_set, new_db, insname = name;
          h_set, ins_table, name, new_db;
          oifits_insert, dest, new_db;
        }
      } else if (type == OIFITS_TYPE_CORR) {
        name = oifits_fix_name(oifits_get_corrname(src, db));
        if (h_has(corr_table, name)) {
          msg = _oifits_compare_tables(h_get(corr_table, name), db, atol, rtol);
          if (msg) {
            error, swrite(format="cannot merge OI_CORR \"%s\" from \"%s\" (%s)",
                          name, filename, msg);
          }
        } else {
          new_db = oifits_unlink_datablock(db);
          h_set, new_db, corrname = name;
          h_set, corr_table, name, new_db;
          oifits_insert, dest, new_db;
        }
      } else if (type == OIFITS_TYPE_VIS || type == OIFITS_TYPE_VIS2 ||
                 type == OIFITS_TYPE_T3 || type == OIFITS_TYPE_FLUX) {
        oifits_insert, dest, db;
      }
    }
  }
}

local oifits_selected_target;
if (is_void(oifits_selected_target)) oifits_selected_target = string();
func oifits_select_target(src, pattern)
/* DOCUMENT oifits_select_target(master, pattern);
         or oifits_select_target, master, pattern;

     Makes an OI-FITS instance with a single target.  `master` is the source
     OI-FITS instance and `pattern` is a glob-style target name.  The result
     may be empty if the target is not found.  When called as a subroutine, the
     modification is done in-place; otherwise, a new instance is returned.

     External variable `oifits_selected_target` can be used to determine the
     name of the selected target.

   SEE ALSO: oifits_new, oifits_load, oifits_save.
 */
{
  extern oifits_selected_target;
  oifits_selected_target = string();
  oifits_update, src;
  dst = oifits_new();
  get_rows = _oifits_get_rows;

  /* First pass: Search the numerical identifier of the target. */
  src_target_id = -1;
  dst_target_id = 1;
  arrnames = h_new();
  insnames = h_new();
  corrnames = h_new();
  for (db = oifits_first(src); db; db = oifits_next(src, db)) {
    type = oifits_get_type(db);
    if (type != OIFITS_TYPE_TARGET) {
      continue;
    }
    if (src_target_id != -1) {
      error, "more than one OI_TARGET data-block in OI-FITS instance";
    }
    target = strtrim(oifits_get_target(src, db));
    j = where(strglob(pattern, target, case=0, path=3, esc=0));
    n = numberof(j);
    if (n != 1) {
      if (n == 0) {
        break;
      }
      error, "too many matching targets";
    }
    j = j(1);
    oifits_selected_target = target = target(j);
    src_target_id = oifits_get_target_id(src, db)(j);
    oifits_new_target, dst,
      revn = oifits_get_revn(src, db),
      target_id = dst_target_id,
      target = target,
      raep0 = get_rows(oifits_get_raep0(src, db), j),
      decep0 = get_rows(oifits_get_decep0(src, db), j),
      equinox = get_rows(oifits_get_equinox(src, db), j),
      ra_err = get_rows(oifits_get_ra_err(src, db), j),
      dec_err = get_rows(oifits_get_dec_err(src, db), j),
      sysvel = get_rows(oifits_get_sysvel(src, db), j),
      veltyp = get_rows(oifits_get_veltyp(src, db), j),
      veldef = get_rows(oifits_get_veldef(src, db), j),
      pmra = get_rows(oifits_get_pmra(src, db), j),
      pmdec = get_rows(oifits_get_pmdec(src, db), j),
      pmra_err = get_rows(oifits_get_pmra_err(src, db), j),
      pmdec_err = get_rows(oifits_get_pmdec_err(src, db), j),
      parallax = get_rows(oifits_get_parallax(src, db), j),
      para_err = get_rows(oifits_get_para_err(src, db), j),
      spectyp = get_rows(oifits_get_spectyp(src, db), j),
      category = get_rows(oifits_get_category(src, db), j);
  }
  if (src_target_id < 0) {
    error, ("no OI_TARGET data-block in OI-FITS instance matching \""
            + pattern + "\"");
  }

  /* Second pass: Append all relevant datablocks. */
  for (db = oifits_first(src); db; db = oifits_next(src, db)) {
    type = oifits_get_type(db);
    if (type == OIFITS_TYPE_VIS || type == OIFITS_TYPE_VIS2 ||
        type == OIFITS_TYPE_T3 || type == OIFITS_TYPE_FLUX) {
      j = where(oifits_get_target_id(src, db) == src_target_id);
      n = numberof(j);
      if (n < 1) continue;
      arrname = oifits_get_arrname(src, db);
      h_set, arrnames, arrname, 1n;
      insname = oifits_get_insname(src, db);
      h_set, insnames, insname, 1n;
      corrname = oifits_get_corrname(src, db);
      if (! is_void(corrname)) {
        h_set, corrnames, corrname, 1n;
      }
      if (type == OIFITS_TYPE_VIS) {
        oifits_new_vis, dst,
          revn = oifits_get_revn(src, db),
          date_obs = oifits_get_date_obs(src, db),
          arrname = arrname,
          insname = insname,
          corrname = corrname,
          amptyp = oifits_get_amptyp(src, db),
          phityp = oifits_get_phityp(src, db),
          amporder = oifits_get_amporder(src, db),
          phiorder = oifits_get_phiorder(src, db),
          target_id = array(dst_target_id, n),
          time = get_rows(oifits_get_time(src, db), j),
          mjd = get_rows(oifits_get_mjd(src, db), j),
          int_time = get_rows(oifits_get_int_time(src, db), j),
          visamp = get_rows(oifits_get_visamp(src, db), j),
          visamperr = get_rows(oifits_get_visamperr(src, db), j),
          corrindx_visamp = get_rows(oifits_get_corrindx_visamp(src, db), j),
          visphi = get_rows(oifits_get_visphi(src, db), j),
          visphierr = get_rows(oifits_get_visphierr(src, db), j),
          corrindx_visphi = get_rows(oifits_get_corrindx_visphi(src, db), j),
          visrefmap = get_rows(oifits_get_visrefmap(src, db), j),
          rvis = get_rows(oifits_get_rvis(src, db), j),
          rviserr = get_rows(oifits_get_rviserr(src, db), j),
          corrindx_rvis = get_rows(oifits_get_corrindx_rvis(src, db), j),
          ivis = get_rows(oifits_get_ivis(src, db), j),
          iviserr = get_rows(oifits_get_iviserr(src, db), j),
          corrindx_ivis = get_rows(oifits_get_corrindx_ivis(src, db), j),
          ucoord = get_rows(oifits_get_ucoord(src, db), j),
          vcoord = get_rows(oifits_get_vcoord(src, db), j),
          sta_index = get_rows(oifits_get_sta_index(src, db), j),
          flag = get_rows(oifits_get_flag(src, db), j);
      } else if (type == OIFITS_TYPE_VIS2) {
        oifits_new_vis2, dst,
          revn = oifits_get_revn(src, db),
          date_obs = oifits_get_date_obs(src, db),
          arrname = arrname,
          insname = insname,
          corrname = corrname,
          target_id = array(dst_target_id, n),
          time = get_rows(oifits_get_time(src, db), j),
          mjd = get_rows(oifits_get_mjd(src, db), j),
          int_time = get_rows(oifits_get_int_time(src, db), j),
          vis2data = get_rows(oifits_get_vis2data(src, db), j),
          vis2err = get_rows(oifits_get_vis2err(src, db), j),
          corrindx_vis2data = get_rows(oifits_get_corrindx_vis2data(src, db), j),
          ucoord = get_rows(oifits_get_ucoord(src, db), j),
          vcoord = get_rows(oifits_get_vcoord(src, db), j),
          sta_index = get_rows(oifits_get_sta_index(src, db), j),
          flag = get_rows(oifits_get_flag(src, db), j);
      } else if (type == OIFITS_TYPE_T3) {
        oifits_new_t3, dst,
          revn = oifits_get_revn(src, db),
          date_obs = oifits_get_date_obs(src, db),
          arrname = arrname,
          insname = insname,
          corrname = corrname,
          target_id = array(dst_target_id, n),
          time = get_rows(oifits_get_time(src, db), j),
          mjd = get_rows(oifits_get_mjd(src, db), j),
          int_time = get_rows(oifits_get_int_time(src, db), j),
          t3amp = get_rows(oifits_get_t3amp(src, db), j),
          t3amperr = get_rows(oifits_get_t3amperr(src, db), j),
          corrindx_t3amp = get_rows(oifits_get_corrindx_t3amp(src, db), j),
          t3phi = get_rows(oifits_get_t3phi(src, db), j),
          t3phierr = get_rows(oifits_get_t3phierr(src, db), j),
          corrindx_t3phi = get_rows(oifits_get_corrindx_t3phi(src, db), j),
          u1coord = get_rows(oifits_get_u1coord(src, db), j),
          v1coord = get_rows(oifits_get_v1coord(src, db), j),
          u2coord = get_rows(oifits_get_u2coord(src, db), j),
          v2coord = get_rows(oifits_get_v2coord(src, db), j),
          sta_index = get_rows(oifits_get_sta_index(src, db), j),
          flag = get_rows(oifits_get_flag(src, db), j);
      } else if (type == OIFITS_TYPE_FLUX) {
        oifits_new_flux, dst,
          revn = oifits_get_revn(src, db),
          date_obs = oifits_get_date_obs(src, db),
          arrname = arrname,
          insname = insname,
          corrname = corrname,
          fov = oifits_get_fov(src, db),
          fovtype = oifits_get_fovtype(src, db),
          calstat = oifits_get_calstat(src, db),
          target_id = array(dst_target_id, n),
          mjd = get_rows(oifits_get_mjd(src, db), j),
          int_time = get_rows(oifits_get_int_time(src, db), j),
          fluxdata = get_rows(oifits_get_fluxdata(src, db), j),
          fluxdata_units = oifits_get_fluxdata_units(src, db),
          fluxerr = get_rows(oifits_get_fluxerr(src, db), j),
          fluxerr_units = oifits_get_fluxerr_units(src, db),
          corrindx_fluxdata = get_rows(oifits_get_corrindx_fluxdata(src, db), j),
          sta_index = get_rows(oifits_get_sta_index(src, db), j),
          flag = get_rows(oifits_get_flag(src, db), j);
      }
    }
  }

  /* Third pass: Deal with non-data entries. */
  for (db = oifits_first(src); db; db = oifits_next(src, db)) {
    type = oifits_get_type(db);
    if (type == OIFITS_TYPE_WAVELENGTH &&
        insnames(oifits_get_insname(src, db))) {
      oifits_new_wavelength, dst,
        revn = oifits_get_revn(src, db),
        insname = oifits_get_insname(src, db),
        eff_wave = oifits_get_eff_wave(src, db),
        eff_band = oifits_get_eff_band(src, db);
    } else if (type == OIFITS_TYPE_ARRAY &&
               arrnames(oifits_get_arrname(src, db))) {
      oifits_new_array, dst,
        revn = oifits_get_revn(src, db),
        arrname = oifits_get_arrname(src, db),
        frame = oifits_get_frame(src, db),
        arrayx = oifits_get_arrayx(src, db),
        arrayy = oifits_get_arrayy(src, db),
        arrayz = oifits_get_arrayz(src, db),
        tel_name = oifits_get_tel_name(src, db),
        sta_name = oifits_get_sta_name(src, db),
        sta_index = oifits_get_sta_index(src, db),
        diameter = oifits_get_diameter(src, db),
        staxyz = oifits_get_staxyz(src, db),
        fov = oifits_get_fov(src, db),
        fovtype = oifits_get_fovtype(src, db);
    } else if (type == OIFITS_TYPE_INSPOL &&
               arrnames(oifits_get_arrname(src, db))) {
      j = where(oifits_get_target_id(src, db) == src_target_id);
      n = numberof(j);
      if (n < 1) continue;
      oifits_new_inspol, dst,
        revn = oifits_get_revn(src, db),
        date_obs = oifits_get_date_obs(src, db),
        npol = oifits_get_npol(src, db),
        arrname = oifits_get_arrname(src, db),
        orient = oifits_get_orient(src, db),
        model = oifits_get_model(src, db),
        target_id = array(dst_target_id, n),
        insname = get_rows(oifits_get_insname(src, db), j),
        mjd_obs = get_rows(oifits_get_mjd_obs(src, db), j),
        mjd_end = get_rows(oifits_get_mjd_end(src, db), j),
        jxx = get_rows(oifits_get_jxx(src, db), j),
        jyy = get_rows(oifits_get_jyy(src, db), j),
        jxy = get_rows(oifits_get_jxy(src, db), j),
        jyx = get_rows(oifits_get_jyx(src, db), j),
        sta_index = get_rows(oifits_get_sta_index(src, db), j);
    } else if (type == OIFITS_TYPE_CORR &&
               corrnames(oifits_get_corrname(src, db))) {
      oifits_new_flux, dst,
        revn = oifits_get_revn(src, db),
        corrname = oifits_get_corrname(src, db),
        ndata = oifits_get_ndata(src, db),
        iindx = oifits_get_iindx(src, db),
        jindx = oifits_get_jindx(src, db),
        corr = oifits_get_corr(src, db);
    }
  }

  if (am_subroutine()) {
    keys = h_keys(src);
    for (k = numberof(keys); k >= 1; --k) {
      h_pop, src, keys(k);
    }
    keys = h_keys(dst);
    for (k = numberof(keys); k >= 1; --k) {
      key = keys(k);
      h_set, src, key, h_pop(dst, key);
    }
    oifits_update, src, force=1n;
  }
  return dst;
}

func _oifits_get_rows(arg, i)
{
  return is_void(arg) ? [] : arg(i,..);
}

/*---------------------------------------------------------------------------*/
/* CONSTRUCTORS FOR ALL DATA-BLOCK TYPES */

/* IMPORTANT: all units are the one of OI-FITS standard. */

func oifits_new_target(master,
                       revn=,
                       target_id=,
                       target=,
                       raep0=,
                       decep0=,
                       equinox=,
                       ra_err=,
                       dec_err=,
                       sysvel=,
                       veltyp=,
                       veldef=,
                       pmra=,
                       pmdec=,
                       pmra_err=,
                       pmdec_err=,
                       parallax=,
                       para_err=,
                       spectyp=,
                       category=)
/* DOCUMENT db = oifits_new_target(key1 = value1, ...);
         or db = oifits_new_target(master, key1 = value1, ...);
         or oifits_new_target, master, key1 = value1, ...;

     Creates a new OI-FITS data-block of type OI_TARGET.  If MASTER is
     non-nil, it must be an existing OI-FITS master instance to store the new
     data-block (in this case the returned value can be ignored).  All members
     of DB are specified by keywords (KEY1 = VALUE1, KEY2 = VALUE2, etc).

     As of revision 2 of OI-FITS standard, the members of an OI_TARGET
     data-block are:

     Keyword    Units   Description
     -------------------------------------------------------------------------
     revn               revision number (default is last version)
     target_id          index number
     target             target name
     raep0      deg     RA at mean equinox
     decep0     deg     DEC at mean equinox
     equinox    yr      equinox
     ra_err     deg     error in RA at mean equinox
     dec_err    deg     error in DEC at mean equinox
     sysvel     m/s     systemic radial velocity
     veltyp             reference for radial velocity ("RADIO" or "OPTICAL")
     veldef             definition of radial velocity ("SR", "HELIOCEN",
                        "BARYCENT", "GEOCENTR", "TOPOCENT", or "UNKNOWN")
     pmra       deg/yr  proper motion in RA
     pmdec      deg/yr  proper motion in DEC
     pmra_err   deg/yr  error of proper motion in RA
     pmdec_err  deg/yr  error of proper motion in DEC
     parallax   deg     parallax
     para_err   deg     error in parallax
     spectyp            spectral type
     category           "CAL" or "SCI"
     -------------------------------------------------------------------------

   SEE ALSO: oifits_insert, oifits_get, oifits_new,
             oifits_new_array, oifits_new_wavelength,
             oifits_new_vis, oifits_new_vis2, oifits_new_t3. */
{
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_stop;
  db = _oifits_datablock_builder(OIFITS_TYPE_TARGET,
                                 h_new(revn = (is_void(revn) ? 2 : revn),
                                       target_id = target_id,
                                       target = target,
                                       raep0 = raep0,
                                       decep0 = decep0,
                                       equinox = equinox,
                                       ra_err = ra_err,
                                       dec_err = dec_err,
                                       sysvel = sysvel,
                                       veltyp = veltyp,
                                       veldef = veldef,
                                       pmra = pmra,
                                       pmdec = pmdec,
                                       pmra_err = pmra_err,
                                       pmdec_err = pmdec_err,
                                       parallax = parallax,
                                       para_err = para_err,
                                       spectyp = spectyp,
                                       category = category));
  if (numberof(_oifits_error_stack)) {
    _oifits_report_error;
  }
  if (master) {
    oifits_insert, master, db;
  }
  return db;
}

func oifits_new_array(master,
                      revn=,
                      arrname=,
                      frame=,
                      arrayx=,
                      arrayy=,
                      arrayz=,
                      tel_name=,
                      sta_name=,
                      sta_index=,
                      diameter=,
                      staxyz=,
                      fov=,
                      fovtype=)
/* DOCUMENT db = oifits_new_array(key1 = value1, ...);
         or db = oifits_new_array(master, key1 = value1, ...);
         or oifits_new_array, master, key1 = value1, ...;

     Creates a new OI-FITS data-block of type OI_ARRAY (position of
     telescopes).  If MASTER is non-nil, it must be an existing OI-FITS master
     instance to store the new data-block (in this case the returned value can
     be ignored).  All members of DB are specified by keywords (KEY1 = VALUE1,
     KEY2 = VALUE2, etc).

     As of revision 2 of OI-FITS standard, the members of an OI_ARRAY
     data-block are:

     Keyword    Units   Description
     ---------------------------------------------------------------
     revn               revision number (default is last version)
     arrname            array name for cross-referencing
     frame              coordinate frame
     arrayx     m       array center X-coordinate
     arrayy     m       array center Y-coordinate
     arrayz     m       array center Z-coordinate
     tel_name           telescope name
     sta_name           station name
     sta_index          station index
     diameter   m       element diameter
     staxyz     m       station coordinates relative to array center
     fov        arcsec  photometric field of view
     fovtype            model for FOV: "FWHM" or "RADIUS"
     ---------------------------------------------------------------

   SEE ALSO: oifits_insert, oifits_get, oifits_new, oifits_new_target,
             oifits_new_wavelength, oifits_new_vis, oifits_new_vis2,
             oifits_new_t3. */
{
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_stop;
  db = _oifits_datablock_builder(OIFITS_TYPE_ARRAY,
                                 h_new(revn = (is_void(revn) ? 2 : revn),
                                       arrname = arrname,
                                       frame = frame,
                                       arrayx = arrayx,
                                       arrayy = arrayy,
                                       arrayz = arrayz,
                                       tel_name = tel_name,
                                       sta_name = sta_name,
                                       sta_index = sta_index,
                                       diameter = diameter,
                                       staxyz = staxyz,
                                       fov = fov,
                                       fovtype = fovtype));
  if (numberof(_oifits_error_stack)) {
    _oifits_report_error;
  }
  if (master) {
    oifits_insert, master, db;
  }
  return db;
}

func oifits_new_wavelength(master,
                           revn=,
                           insname=,
                           eff_wave=,
                           eff_band=)
/* DOCUMENT db = oifits_new_array(key1 = value1, ...);
         or db = oifits_new_array(master, key1 = value1, ...);
         or oifits_new_array, master, key1 = value1, ...;

     Creates a new OI-FITS data-block of type OI_WAVELENGTH (spectral
     channels).  If MASTER is non-nil, it must be an existing OI-FITS
     master instance to store the new data-block (in this case the returned
     value can be ignored).  All members of DB are specified by keywords
     (KEY1 = VALUE1, KEY2 = VALUE2, etc).

     As of revision 2 of OI-FITS standard, the members of an OI_WAVELENGTH
     data-block are:

     Keyword    Units   Description
     ------------------------------------------------------------
     revn               revision number (default is last version)
     insname            name of detector for cross-referencing
     eff_wave   m       effective wavelength of channel
     eff_band   m       effective bandpass of channel
     ------------------------------------------------------------

   SEE ALSO: oifits_insert, oifits_get, oifits_new, oifits_new_array,
             oifits_new_target, oifits_new_vis, oifits_new_vis2,
             oifits_new_t3. */
{
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_stop;
  db = _oifits_datablock_builder(OIFITS_TYPE_WAVELENGTH,
                                 h_new(revn = (is_void(revn) ? 2 : revn),
                                       insname = insname,
                                       eff_wave = eff_wave,
                                       eff_band = eff_band));
  if (numberof(_oifits_error_stack)) {
    _oifits_report_error;
  }
  if (master) {
    oifits_insert, master, db;
  }
  return db;
}

local oifits_new_vis;
/* DOCUMENT db = oifits_new_vis(key1 = value1, ...);
         or db = oifits_new_vis(master, key1 = value1, ...);
         or oifits_new_vis, master, key1 = value1, ...;

     Creates a new OI-FITS data-block of type OI_VIS (complex visibility).  If
     MASTER is non-nil, it must be an existing OI-FITS master instance to
     store the new data-block (in this case the returned value can be
     ignored).  All members of DB are specified by keywords (KEY1 = VALUE1,
     KEY2 = VALUE2, etc).

     As of revision 2 of OI-FITS standard, the members of an OI_VIS data-block
     are:

     Keyword       Units  Description
     --------------------------------------------------------------------------
     revn                 revision number (default is last version)
     date_obs             UTC start date of observations
     arrname              name of corresponding OI_ARRAY table
     insname              name of corresponding OI_WAVELENGTH table
     corrname             name of corresponding OI_CORR table
     amptyp               'absolute', 'differential', or 'correlated flux'
     phityp               'absolute', or 'differential'",
     amporder             polynomial fit order for differential chromatic amplitudes
     phiorder             polynomial fit order for differential chromatic phases
     target_id            target number as index into OI_TARGET table
     time            s    UTC time of observation
     mjd             day  modified Julian Day
     int_time        s    integration time
     visamp               complex visibility amplitude
     visamperr            error in complex visibility amplitude
     corrindx_visamp      index into correlation matrix for 1st VISAMP element
     visphi          deg  complex visibility phase
     visphierr       deg  error in complex visibility phase
     corrindx_visphi      index into correlation matrix for 1st VISPHI element
     visrefmap            boolean matrix indicating which spectral channels were
                          taken as reference for differential chromatic visibility
                          computation
     rvis                 real part of complex coherent flux
     rviserr              error on RVIS
     corrindx_rvis        index into correlation matrix for 1st RVIS element
     ucoord          m    U coordinate of the data
     vcoord          m    V coordinate of the data
     sta_index            station numbers contributing to the data
     flag                 flag
     --------------------------------------------------------------------------

   SEE ALSO: oifits_insert, oifits_get, oifits_new, oifits_new_array,
             oifits_new_target, oifits_new_wavelength, oifits_new_vis2,
             oifits_new_t3. */
func oifits_new_vis(master,
                    revn=,
                    date_obs=,
                    arrname=,
                    insname=,
                    corrname=,
                    amptyp=,
                    phityp=,
                    amporder=,
                    phiorder=,
                    target_id=,
                    time=,
                    mjd=,
                    int_time=,
                    visamp=,
                    visamperr=,
                    corrindx_visamp=,
                    visphi=,
                    visphierr=,
                    corrindx_visphi=,
                    visrefmap=,
                    rvis=,
                    rviserr=,
                    corrindx_rvis=,
                    ivis=,
                    iviserr=,
                    corrindx_ivis=,
                    ucoord=,
                    vcoord=,
                    sta_index=,
                    flag=)
{
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_stop;
  db = _oifits_datablock_builder(OIFITS_TYPE_VIS,
                                 h_new(revn = (is_void(revn) ? 2 : revn),
                                       date_obs = date_obs,
                                       arrname = arrname,
                                       insname = insname,
                                       corrname = corrname,
                                       amptyp = amptyp,
                                       phityp = phityp,
                                       amporder = amporder,
                                       phiorder = phiorder,
                                       target_id = target_id,
                                       time = time,
                                       mjd = mjd,
                                       int_time = int_time,
                                       visamp = visamp,
                                       visamperr = visamperr,
                                       corrindx_visamp = corrindx_visamp,
                                       visphi = visphi,
                                       visphierr = visphierr,
                                       corrindx_visphi = corrindx_visphi,
                                       visrefmap = visrefmap,
                                       rvis = rvis,
                                       rviserr = rviserr,
                                       corrindx_rvis = corrindx_rvis,
                                       ivis = ivis,
                                       iviserr = iviserr,
                                       corrindx_ivis = corrindx_ivis,
                                       ucoord = ucoord,
                                       vcoord = vcoord,
                                       sta_index = sta_index,
                                       flag = flag));
  if (numberof(_oifits_error_stack)) {
    _oifits_report_error;
  }
  if (master) {
    oifits_insert, master, db;
  }
  return db;
}

local oifits_new_vis2;
/* DOCUMENT db = oifits_new_vis2(key1 = value1, ...);
         or db = oifits_new_vis2(master, key1 = value1, ...);
         or oifits_new_vis2, master, key1 = value1, ...;

     Creates a new OI-FITS data-block of type OI_VIS2 (squared visibility
     data).  If MASTER is non-nil, it must be an existing OI-FITS master
     instance to store the new data-block (in this case the returned value can
     be ignored).  All members of DB are specified by keywords (KEY1 = VALUE1,
     KEY2 = VALUE2, etc).

     As of revision 2 of OI-FITS standard, the members of an OI_VIS2
     data-block are:

     Keyword          Units  Description
     --------------------------------------------------------------------------
     revn                    revision number (default is last version)
     date_obs                UTC start date of observations
     arrname                 name of corresponding OI_ARRAY table
     insname                 name of corresponding OI_WAVELENGTH table
     corrname                name of corresponding OI_CORR table
     target_id               target number as index into OI_TARGET table
     time               s    UTC time of observation
     mjd                day  modified Julian Day
     int_time           s    integration time
     vis2data                squared visibility
     vis2err                 error in squared visibility
     corrindx_vis2data       index into correlation matrix for 1st VIS2DATA element
     ucoord             m    U coordinate of the data
     vcoord             m    V coordinate of the data
     sta_index               station numbers contributing to the data
     flag                    flag
     --------------------------------------------------------------------------

   SEE ALSO: oifits_insert, oifits_get, oifits_new, oifits_new_array,
             oifits_new_target, oifits_new_wavelength, oifits_new_vis,
             oifits_new_t3. */
func oifits_new_vis2(master,
                     revn=,
                     date_obs=,
                     arrname=,
                     insname=,
                     corrname=,
                     target_id=,
                     time=,
                     mjd=,
                     int_time=,
                     vis2data=,
                     vis2err=,
                     corrindx_vis2data=,
                     ucoord=,
                     vcoord=,
                     sta_index=,
                     flag=)
{
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_stop;
  db = _oifits_datablock_builder(OIFITS_TYPE_VIS2,
                                 h_new(revn = (is_void(revn) ? 2 : revn),
                                       date_obs = date_obs,
                                       arrname = arrname,
                                       insname = insname,
                                       corrname = corrname,
                                       target_id = target_id,
                                       time = time,
                                       mjd = mjd,
                                       int_time = int_time,
                                       vis2data = vis2data,
                                       vis2err = vis2err,
                                       corrindx_vis2data = corrindx_vis2data,
                                       ucoord = ucoord,
                                       vcoord = vcoord,
                                       sta_index = sta_index,
                                       flag = flag));
  if (numberof(_oifits_error_stack)) {
    _oifits_report_error;
  }
  if (master) {
    oifits_insert, master, db;
  }
  return db;
}

local oifits_new_t3;
/* DOCUMENT db = oifits_new_t3(key1 = value1, ...);
         or db = oifits_new_t3(master, key1 = value1, ...);
         or oifits_new_t3, master, key1 = value1, ...;

     Creates a new OI-FITS data-block of type OI_T3 (bispectrum or triple
     correlation data).  If MASTER is non-nil, it must be an existing OI-FITS
     master instance to store the new data-block (in this case the returned
     value can be ignored).  All members of DB are specified by keywords (KEY1
     = VALUE1, KEY2 = VALUE2, etc).

     As of revision 2 of OI-FITS standard, the members of an OI_T3 data-block
     are:

     Keyword       Units  Description
     ------------------------------------------------------------------------
     revn                 revision number (default is last version)
     date_obs             UTC start date of observations
     arrname              name of corresponding OI_ARRAY table
     insname              name of corresponding OI_WAVELENGTH table
     corrname             name of corresponding OI_CORR table
     target_id            target number as index into OI_TARGET table
     time            s    UTC time of observation
     mjd             day  modified Julian Day
     int_time        s    integration time
     t3amp                triple product amplitude
     t3amperr             error in triple product amplitude
     corrindx_t3amp       index into correlation matrix for 1st T3AMP element
     t3phi           deg  triple product phase
     t3phierr        deg  error in triple product phase
     corrindx_t3phi       index into correlation matrix for 1st T3PHI element
     u1coord         m    U coordinate of baseline AB of the triangle
     v1coord         m    V coordinate of baseline AB of the triangle
     u2coord         m    U coordinate of baseline BC of the triangle
     v2coord         m    V coordinate of baseline BC of the triangle
     sta_index            station numbers contributing to the data
     flag                 flag
     ------------------------------------------------------------------------

   SEE ALSO: oifits_insert, oifits_get, oifits_new, oifits_new_array,
             oifits_new_target, oifits_new_wavelength, oifits_new_vis,
             oifits_new_vis2. */
func oifits_new_t3(master,
                   revn=,
                   date_obs=,
                   arrname=,
                   insname=,
                   corrname=,
                   target_id=,
                   time=,
                   mjd=,
                   int_time=,
                   t3amp=,
                   t3amperr=,
                   corrindx_t3amp=,
                   t3phi=,
                   t3phierr=,
                   corrindx_t3phi=,
                   u1coord=,
                   v1coord=,
                   u2coord=,
                   v2coord=,
                   sta_index=,
                   flag=)
{
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_stop;
  db = _oifits_datablock_builder(OIFITS_TYPE_T3,
                                 h_new(revn = (is_void(revn) ? 2 : revn),
                                       date_obs = date_obs,
                                       arrname = arrname,
                                       insname = insname,
                                       corrname = corrname,
                                       target_id = target_id,
                                       time = time,
                                       mjd = mjd,
                                       int_time = int_time,
                                       t3amp = t3amp,
                                       t3amperr = t3amperr,
                                       corrindx_t3amp = corrindx_t3amp,
                                       t3phi = t3phi,
                                       t3phierr = t3phierr,
                                       corrindx_t3phi = corrindx_t3phi,
                                       u1coord = u1coord,
                                       v1coord = v1coord,
                                       u2coord = u2coord,
                                       v2coord = v2coord,
                                       sta_index = sta_index,
                                       flag = flag));
  if (numberof(_oifits_error_stack)) {
    _oifits_report_error;
  }
  if (master) {
    oifits_insert, master, db;
  }
  return db;
}

local oifits_new_flux;
/* DOCUMENT db = oifits_new_flux(key1 = value1, ...);
         or db = oifits_new_flux(master, key1 = value1, ...);
         or oifits_new_flux, master, key1 = value1, ...;

     Creates a new OI-FITS data-block of type OI_FLUX (spectral energy
     distribution).  If MASTER is non-nil, it must be an existing OI-FITS
     master instance to store the new data-block (in this case the returned
     value can be ignored).  All members of DB are specified by keywords (KEY1
     = VALUE1, KEY2 = VALUE2, etc).

     As of revision 2 of OI-FITS standard, the members of an OI_FLUX
     data-block are:

     Keyword    Units   Description
     -------------------------------------------------------------------------
     revn               revision number (default is last version)
     date_obs           UTC start date of observations
     insname            name of corresponding detector
     arrname            name of corresponding array
     corrname           name of corresponding OI_CORR table
     fov        arcsec  area of sky over which flux is integrated
     fovtype            model for FOV: "FWHM" or "RADIUS"
     calstat            "C": spectrum is calibrated, "U": uncalibrated
     target_id          target number as index into OI_TARGET table
     mjd        day     modified Julian Day
     int_time   s       integration time
     flux               target spectral energy distribution
     fluxerr            flux error
     corrindx_fluxdata  index into correlation matrix for 1st FLUXDATA element
     sta_index          station numbers contributing to the data
     flag               flag
     -------------------------------------------------------------------------

   SEE ALSO: oifits_insert, oifits_get, oifits_new, oifits_new_array,
             oifits_new_target, oifits_new_wavelength. */
func oifits_new_flux(master,
                     revn=,
                     date_obs=,
                     insname=,
                     arrname=,
                     corrname=,
                     fov=,
                     fovtype=,
                     calstat=,
                     target_id=,
                     mjd=,
                     int_time=,
                     fluxdata=,
                     fluxdata_units=,
                     fluxerr=,
                     fluxerr_units=,
                     corrindx_fluxdata=,
                     sta_index=,
                     flag=)
{
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_stop;
  db = _oifits_datablock_builder(OIFITS_TYPE_FLUX,
                                 h_new(revn = (is_void(revn) ? 1 : revn),
                                       date_obs = date_obs,
                                       insname = insname,
                                       arrname = arrname,
                                       corrname = corrname,
                                       fov = fov,
                                       fovtype = fovtype,
                                       calstat = calstat,
                                       target_id = target_id,
                                       mjd = mjd,
                                       int_time = int_time,
                                       fluxdata = fluxdata,
                                       fluxdata_units = (is_void(fluxdata_units) ? "UNKNOWN" : fluxdata_units),
                                       fluxerr = fluxerr,
                                       fluxerr_units = (is_void(fluxerr_units) ? "UNKNOWN" : fluxerr_units),
                                       corrindx_fluxdata = corrindx_fluxdata,
                                       sta_index = sta_index,
                                       flag = flag));
  if (numberof(_oifits_error_stack)) {
    _oifits_report_error;
  }
  if (master) {
    oifits_insert, master, db;
  }
  return db;
}

local oifits_new_inspol;
/* DOCUMENT db = oifits_new_inspol(key1 = value1, ...);
         or db = oifits_new_inspol(master, key1 = value1, ...);
         or oifits_new_inspol, master, key1 = value1, ...;

     Creates a new OI-FITS data-block of type OI_INSPOL (beam polarization).
     If MASTER is non-nil, it must be an existing OI-FITS master instance to
     store the new data-block (in this case the returned value can be
     ignored).  All members of DB are specified by keywords (KEY1 = VALUE1,
     KEY2 = VALUE2, etc).

     As of revision 2 of OI-FITS standard, the members of an OI_INSPOL
     data-block are:

     Keyword    Units  Description
     ---------------------------------------------------------------------
     revn              revision number (default is last version)
     date_obs          UTC start date of observations
     npol              number of polarisation types in this table
     arrname           name of corresponding OI_ARRAY table
     orient            orientation of the Jones matrix, could be "NORTH"
                       (for on-sky orientation), or "LABORATORY"
     model             method for estimating the Jones matrix
     target_id         target number as index into OI_TARGET table
     insname           name of corresponding OI_WAVELENGTH table
     mjd_obs     day   modified Julian Day, start of time lapse
     mjd_end     day   modified Julian Day, end of time lapse
     jxx               complex Jones matrix component along X axis
     jyy               complex Jones matrix component along Y axis
     jxy               complex Jones matrix component between X and Y axis
     jyx               complex Jones matrix component between Y and X axis
     sta_index         station number for the Jones matrix
     ---------------------------------------------------------------------

   SEE ALSO: oifits_insert, oifits_get, oifits_new, oifits_new_array,
             oifits_new_target, oifits_new_wavelength. */
func oifits_new_inspol(master,
                       revn=,
                       date_obs=,
                       npol=,
                       arrname=,
                       orient=,
                       model=,
                       target_id=,
                       mjd_obs=,
                       mjd_end=,
                       jxx=,
                       jyy=,
                       jxy=,
                       jyx=,
                       sta_index=)
{
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_stop;
  db = _oifits_datablock_builder(OIFITS_TYPE_INSPOL,
                                 h_new(revn = (is_void(revn) ? 1 : revn),
                                       date_obs = date_obs,
                                       npol = npol,
                                       arrname = arrname,
                                       orient = orient,
                                       model = model,
                                       target_id = target_id,
                                       insname = insname,
                                       mjd_obs = mjd_obs,
                                       mjd_end = mjd_end,
                                       jxx = jxx,
                                       jyy = jyy,
                                       jxy = jxy,
                                       jyx = jyx,
                                       sta_index = sta_index));
  if (numberof(_oifits_error_stack)) {
    _oifits_report_error;
  }
  if (master) {
    oifits_insert, master, db;
  }
  return db;
}

local oifits_new_corr;
/* DOCUMENT db = oifits_new_corr(key1 = value1, ...);
         or db = oifits_new_corr(master, key1 = value1, ...);
         or oifits_new_corr, master, key1 = value1, ...;

     Creates a new OI-FITS data-block of type OI_CORR (correlations).
     If MASTER is non-nil, it must be an existing OI-FITS master instance to
     store the new data-block (in this case the returned value can be
     ignored).  All members of DB are specified by keywords (KEY1 = VALUE1,
     KEY2 = VALUE2, etc).

     As of revision 2 of OI-FITS standard, the members of an OI_CORR
     data-block are:

     Keyword  Units  Description
     ----------------------------------------------------------
     revn            revision number (default is last version)
     corrname        name of correlated data set
     ndata           number of correlated data
     iindx           first index of correlation matrix element
     jindx           second index of correlation matrix element
     corr            matrix element at (iindx,jindx)
     ----------------------------------------------------------

   SEE ALSO: oifits_insert, oifits_get, oifits_new, oifits_new_array,
             oifits_new_target, oifits_new_wavelength. */
func oifits_new_corr(master,
                     revn=,
                     corrname=,
                     ndata=,
                     iindx=,
                     jindx=,
                     corr=)
{
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_stop;
  db = _oifits_datablock_builder(OIFITS_TYPE_CORR,
                                 h_new(revn = (is_void(revn) ? 1 : revn),
                                       corrname = corrname,
                                       ndata = ndata,
                                       iindx = iindx,
                                       jindx = jindx,
                                       corr = corr));
  if (numberof(_oifits_error_stack)) {
    _oifits_report_error;
  }
  if (master) {
    oifits_insert, master, db;
  }
  return db;
}

/* NOTE: To simplify the code of the instance builder, we heavily rely on
   class definition tables.  This means that consistency of these tables is
   assumed and must be asserted at initialization time when _oifits_init is
   called. */

func _oifits_datablock_builder(type, src, extname, hdu)
/**DOCUMENT db = _oifits_datablock_builder(type, tbl);
         or db = _oifits_datablock_builder(type, fh, extname, hdu);

     Returns an instanciated OI-FITS data-block with given TYPE (see
     oifits_get_type) and with contents taken from hash-table TBL or read from
     OI-FITS file handle FH.  In the later case, EXTNAME and HDU are the FITS
     extension name and header-data unit number.

     This private routine is intended to centralize most of the checkings to
     insure consistent definitions of OI-FITS data-blocks.

     In case of error, an empty result is returned.  However whether the
     function immediately returns or not depends on the _oifits_on_error
     function which can be redefined (e.g. to use built-in "error" function)
     prior to calling _oifits_datablock_builder.  The stack of error messages
     is stored in external variable _oifits_error_stack.


  SEE ALSO: oifits_insert, oifits_new_target, oifits_new_array,
            oifits_new_wavelength, oifits_new_vis, oifits_new_vis2,
            oifits_new_t3. */
{
  /* Clear error stack. */
  extern _oifits_error_stack;
  _oifits_error_count = 0;
  _oifits_error_stack = [];

  /* Get revision number. */
  reading = (! is_void(extname));
  if (reading) {
    revn = fits_get(src, "OI_REVN");
  } else {
    revn = h_get(src, "revn");
  }
  if (is_void(revn)) {
    if (type == OIFITS_TYPE_FLUX) {
      oifits_warn, "missing OI_REVN has been fixed in OI_FLUX";
      revn = 1;
      if (!reading) h_set, src, revn=revn;
    }
  }
  if (_oifits_get_integer_scalar(revn) ||
      revn < 1 || revn > _OIFITS_REVN_MAX) {
    _oifits_error, "unrecognized or missing revision of the OI-FITS format";
    return;
  }

  /* Get datablock type. */
  if (_oifits_get_integer_scalar(type) || type < 1
      || type > numberof(_OIFITS_DATABLOCK_CLASS)) {
     _oifits_error, "bad OI-FITS data block type";
    return;
  }
  class = _OIFITS_DATABLOCK_CLASS(type);
  if (! symbol_exists(_oifits_classdef_name(class, revn))) {
    error, swrite("invalid revision number for OI_%s", class);
  }
  is_data = (type == OIFITS_TYPE_VIS ||
             type == OIFITS_TYPE_VIS2 ||
             type == OIFITS_TYPE_T3 ||
             type == OIFITS_TYPE_FLUX);

  /* If datablock contents is to be read from OI-FITS file, read all
     the columns of the binary table and build a hash table for fast
     linking of the column index given its name.*/
  if (reading) {
    ptr = fits_read_bintable(src);
    ttype = h_new();
    for (i = numberof(ptr); i >= 1; --i) {
      name = fits_get(src, swrite(format="TTYPE%d", i));
      if (structof(name) == string) {
        h_set, ttype, oifits_fix_name(name), i;
      }
    }
  }

  /* Setup new hash table to store the OI-FITS datablock. */
  db = h_new(__class = class, __type = type, __is_data = is_data,
             revn = revn);

  /* Instanciate data-block object from header part of HDU. */
  local table, value; /* for parsing class definition tables */
  eq_nocopy, table, _oifits_classdef_header(class, revn);
  for (j = numberof(table); j >= 2 /* omit 1st line, i.e. OI_REVN */; --j) {
    /* Parse table entry. */
    entry = table(j);
    member = entry.member;
    keyword = entry.keyword;
    units = entry.units;
    comment = entry.comment;
    multiplier = entry.multiplier;
    ctype = entry.ctype;
    optional = entry.optional;

    /* Prepare for error messages. */
    if (reading) {
      what = swrite(format="keyword '%s' in %s (HDU %d)",
                    keyword, extname, hdu);
    } else {
      what = swrite(format="member '%s' for %s",
                    member, class);
    }

    /* Get value of keyword and check its data type and dimension
       list before memorizing it. */
    if (reading) {
      value = fits_get(src, keyword);
    } else {
      eq_nocopy, value, h_get(src, member);
    }
    if (! is_scalar(value)) {
      if (is_void(value)) {
        if (! optional && _oifits_error("missing mandatory %s", what)) {
          return;
        } else {
          continue;
        }
      } else {
        if (_oifits_error("expecting a scalar for %s", what)) {
          return;
        } else {
          continue;
        }
      }
    }
    if (_oifits_datablock_builder_fix_type(ctype, value)) {
      if (_oifits_error("bad data type for %s", what)) {
        return;
      } else {
        continue;
      }
    }
    h_set, db, member, value;
  }

  /* Instanciate data-block object from data (table) part of HDU. */
  local table, value; /* for parsing class definition tables */
  nrows = 0; /* we don't know yet the number of input rows */
  nwavelengths = 0; /* number of spectral bandwidths */
  eq_nocopy, table, _oifits_classdef_column(class, revn);
  missing = array(int,  numberof(table));
  for (j = 1; j <= numberof(table); ++j) { // <=== BEWARE ORDER IS IMPORTANT!!!
    /* Parse table entry. */
    entry = table(j);
    member = entry.member;
    keyword = entry.keyword;
    multiplier = entry.multiplier;
    ctype = entry.ctype;
    optional = entry.optional;

    /* Prepare for error messages. */
    if (reading) {
      what = swrite(format="column '%s' in %s (HDU %d)", keyword, extname, hdu);
    } else {
      what = swrite(format="member '%s' for %s", member, class);
    }

    /* Get value of field. */
    if (reading) {
      i = ttype(keyword);
      if (is_void(i) || is_void(*ptr(i))) {
        missing(j) = 1n;
        continue;
      }
      eq_nocopy, value, *ptr(i);
    } else {
      eq_nocopy, value, h_get(src, member);
      if (is_void(value)) {
        missing(j) = 1n;
        continue;
      }
    }

    /* Check and fix data type. */
    if (_oifits_datablock_builder_fix_type(ctype, value)) {
      if (_oifits_error("bad data type for %s", what)) {
        return;
      } else {
        continue;
      }
    }

    /* Check correctness of dimension list. */
    dims = dimsof(value);
    rank = dims(1);
    nrows1 = (rank >= 1 ? dims(2) : 1);
    ncells = numberof(value)/nrows1;
    if (nrows == 0) {
      nrows = nrows1;
    } else if (nrows1 != nrows) {
      if (_oifits_error("bad number of rows for %s", what)) {
        return;
      } else {
        continue;
      }
    }
    if (multiplier < 0) {
      /* Each cell has abs(multiplier) dimensions each of length the number of
         spectral channels. */
      multiplier = -multiplier;
      if (nwavelengths == 0) {
        nwavelengths = lround(ncells^(1.0/multiplier));
      }
      if (nwavelengths == 1) {
        /* Apply Yorick rules for unit trailing dimensions and add them if
           they are missing. */
        while (rank < 1 + multiplier) {
          value = value(.., -);
          ++rank;
        }
        dims = dimsof(value);
      }
      err = (rank != 1 + multiplier || anyof(dims(3:) != nwavelengths));
      if (err && keyword == "FLAG" && ncells == 1) {
        /* Hack for AMBER data ;-( */
        err = 0n;
        for (p = 1; p <= multiplier; ++p) {
          value = value(.., -:1:nwavelengths);
        }
      }
    } else {
      err = (rank > 2 || (is_string(value) ? ncells != 1
                          : ncells != multiplier));
    }
    if (err) {
      if (_oifits_error("bad cell dimensions for %s", what)) {
        return;
      } else {
        continue;
      }
    }

    /* Memorize datablock member. */
    if (h_has(db, member)) {
      if (_oifits_error("multiple definitions for %s", what)) {
        return;
      } else {
        continue;
      }
    }
    h_set, db, member, value;

    /* Get units for that column/member. */
    member_units = member + "_units";
    if (reading) {
      given_units = fits_get(src, swrite(format="TUNIT%d", i));
    } else {
      given_units = h_get(src, member_units);
    }
    if (is_void(given_units)) {
      given_units = entry.units;
    }
    if (structof(given_units) == string && given_units != "-") {
      h_set, db, member_units, given_units;
    } else {
      h_pop, db, member_units;
    }
  }

  /* Fix number of baselines and spectral channels. */
  if (is_data) {
    nbaselines = nrows;
  } else {
    nbaselines = 0;
    if (type == OIFITS_TYPE_WAVELENGTH) {
      nwavelengths = nrows;
    }
  }
  h_set, db, __nbaselines = nbaselines, __nwavelengths = nwavelengths;

  /* Deal with missing members (must be done at the end). */
  if (anyof(missing)) {
    missing = where(missing);
    eq_nocopy, table, _oifits_classdef_column(db.__class, revn);
    for (j = numberof(missing); j >= 1; --j) {
      k = missing(j);
      member = table(k).member;
      optional = table(k).optional;
      if (member == "flag" && nbaselines >= 1 && nwavelengths >= 1) {
        h_set, db, flag = array(0n, nbaselines, nwavelengths);
      } else if (! optional){
        if (reading) {
          keyword = table(k).keyword;
          what = swrite(format="column '%s' in %s (HDU %d)",
                        keyword, extname, hdu);
        } else {
          what = swrite(format="member '%s' for %s",
                        member, class);
        }
        if (_oifits_error("missing mandatory %s", what)) {
          return;
        } else {
          continue;
        }
      }
    }
  }
  if (_oifits_error_count == 0) {
    return db;
  }
}

func _oifits_datablock_builder_fix_type(ctype, &value)
{
  s = structof(value);
  if (ctype == _OIFITS_CTYPE_REAL) {
    if (s == double) {
      return 0;
    }
    if (s == float || s == long || s == int ||
        s == short || s == char) {
      value = double(value);
      return 0;
    }
    return -1;
  }
  if (ctype == _OIFITS_CTYPE_INTEGER) {
    if (s == long) {
      return 0;
    }
    if (s == int || s == short || s == char) {
      value = long(value);
      return 0;
    }
    return -1;
  }
  if (ctype == _OIFITS_CTYPE_STRING) {
    return (s == string ? 0 : -1);
  }
  if (ctype == _OIFITS_CTYPE_LOGICAL) {
    return (s == int ? 0 : -1);
  }
  if (ctype == _OIFITS_CTYPE_COMPLEX) {
    if (s == complex) {
      return 0;
    }
    if (s == double || s == float || s == long || s == int ||
        s == short || s == char) {
      value = complex(value);
      return 0;
    }
    return -1;
  }
  return -1;
}

/*---------------------------------------------------------------------------*/
/* LOADING/SAVING AN OI-FITS FILE */

func oifits_load(filename, quiet=, errmode=)
/* DOCUMENT master = oifits_load(filename);
     Returns an opaque handle with the contents of OI-FITS data file FILENAME.

     Keyword ERRMODE can be set to specify error reporting in oifits_update;
     the default is ERRMODE=1.

   SEE ALSO: oifits_save, oifits_new, oifits_insert, oifits_update.
*/
{
  /* Setup for error management. */
  local _oifits_error_stack;
  _oifits_on_error = _oifits_on_error_push;
  if (is_void(errmode)) errmode = 1n;

  /* Load the data from the FITS file. */
  fh = fits_open(filename);
  master = oifits_new();
  counter = 0;
  while (! fits_eof(fh)) {
    hdu = fits_current_hdu(fh);
    xtension = fits_get_xtension(fh);
    nbytes = fits_get_data_size(fh);
    skip = 1n;
    if (xtension == "BINTABLE") {
      extname = fits_get(fh, "EXTNAME");
      if (structof(extname) == string) {
        extname = oifits_fix_name(extname);
        type = _OIFITS_TYPE_TABLE(extname);
        if (! is_void(type)) {
          skip = 0n;
          if (! quiet) {
            write, format="loading %s extension data from unit %d\n",
              extname, hdu;
          }
          db = _oifits_datablock_builder(type, fh, extname, hdu);
          if (! is_void(db)) {
            oifits_insert, master, db;
          } else if (errmode && numberof(_oifits_error_stack)) {
            _oifits_report_error;
          }
        }
      }
    }
    if (skip && ! quiet) {
      write, format="skipping %s %s extension in unit %d\n",
        (nbytes ? "non-empty" : "empty"), xtension, hdu;
    }
    fits_next_hdu, fh;
  }
  if (numberof(_oifits_error_stack)) {
    h_set, master, errmsg=_oifits_error_stack;
  }
  return oifits_update(master, errmode=errmode, quiet=1n);
}

func oifits_save(master, filename, revn=, overwrite=,
                 comment=, history=)
/* DOCUMENT oifits_save, master, filename;

     Write optical interferometric data contained in instance MASTER into
     FILENAME in OI-FITS format.  That is to say the OI_TARGET (only one
     otherwise an error is produced) and eventual OI_ARRAY, OI_WAVELENGTH,
     OI_VIS, OI_VIS2 and OI_T3 data.

     Keyword REVN can be used to specify the revision number of the OI-FITS

   SEE ALSO: oifits_write_hashtable_as_bintable.
*/
{
  local table, value;

  /* Assert consistency of MASTER. */
  if (! is_hash(master)) error, " expected hash table for first argument";
  oifits_update, master, errmode=1, revn=revn;

  /* Create new FITS file with empty primary HDU. */
  grow, comment, "FITS (Flexible Image Transport System) format defined in Astronomy and Astrophysics Supplement Series v44/p363, v44/p371, v73/p359, v73/p365", "Contact the NASA Science Office of Standards and Technology for the FITS Definition document #100 and other FITS information.", "Binary Extensions conform to OI-DATA standard for exchange of optical interferometry data currently described at http://www.mrao.cam.ac.uk/~jsy1001/exchange/.", "This file was created by Yeti/Yorick code from Centre de Recherche Astrophysique de Lyon (CRAL).";
  fh = fits_create(filename, bitpix=8, overwrite=overwrite, extend=1,
                   history=history, comment=comment);
  fits_write_header, fh;

  /* Create a BINTABLE HDU for every OI_FITS datablocks. */
  db_list = master.__target;
  grow, db_list, master.__array, master.__wavelength, master.__data;
  for (i = 1 ; i <= numberof(db_list) ; ++i) {
    key = db_list(i);
    db = master(key);
    class = db.__class;
    extname = "OI_" + class;

    fits_new_bintable, fh;
    fits_set, fh, "EXTNAME", extname;

    /* Write header description of the BINTABLE. */
    revn = oifits_get_revn(master, db);
    eq_nocopy, table, _oifits_classdef_header(class, revn);
    n = numberof(table);
    for (j = 1; j <= n; ++j) {
      entry = table(j);
      eq_nocopy, value, h_get(db, entry.member);
      if (is_void(value)) {
        continue;
      }
      fits_set, fh, entry.keyword, value, entry.comment;
    }

    /* Write column description of the BINTABLE. */
    eq_nocopy, table, _oifits_classdef_column(class, revn);
    n = numberof(table);
    ptr = array(pointer, n);
    ncols = 0;
    for (j = 1; j <= n; ++j) {
      entry = table(j);
      eq_nocopy, value, h_get(db, entry.member);
      if (is_void(value)) {
        continue;
      }
      ++ncols;
      nth = swrite(format="%d", ncols);
      ptr(ncols) = &value;

      /* TTYPE# */
      fits_set, fh, "TTYPE"+nth, entry.keyword, entry.comment;

      /* TFORM# */
      multiplier = entry.multiplier;
      letter = entry.letter;
      if (multiplier > 0) {
        ncells = multiplier;
      } else if (is_string(value)) {
        ncells = max(strlen(value));
      } else {
        dims = dimsof(value);
        ncells = (dims(1) >= 2 ? numberof(value)/dims(2) : 1);
      }
      fits_set, fh, "TFORM"+nth, swrite(format="%d%s", ncells, letter);

      /* TUNIT# */
      units = h_get(db, entry.member + "_units");
      if (! is_string(units) || ! is_scalar(units)) {
        if (! is_void(units)) {
          oifits_warn, swrite(format="bad units for '%s'", member);
        }
        units = entry.units;
      }
      if (strlen(units) && units != "-") {
        fits_set, fh, "TUNIT"+nth, units;
      }
    }
    if (ncols > 0) {
      fits_write_bintable, fh, ptr(1:ncols);
    }
  }
  fits_close, fh;
}


/*---------------------------------------------------------------------------*/
/* ACCESSING CONTENTS OF OI-FITS OBJECTS */

local oifits_first, oifits_last, oifits_next;
/* DOCUMENT oifits_first(master);
         or oifits_last(master);
         or oifits_next(master, db);
     Get first, last or next datablock in OI-FITS handle MASTER.  Useful to run
     across all datablocks.  For instance:

       for (db = oifits_first(master); db; db = oifits_next(master, db)) {
         ...;
       }

    Beware that getting the last datablock involves travelling all the chain
    of datablocks.

   SEE ALSO h_new, h_keys. */

func oifits_first(master) { return h_get(master, master.__first); }

func oifits_next(master, db)
{
  return (is_void(db) || is_void((key = db.__next))) ? [] : h_get(master, key);
}

func oifits_last(master)
{
  db = oifits_first(master);
  while (db) {
    key = db.__next;
    if (is_void(key)) break;
    next = h_get(master, key);
    if (is_void(next)) break;
    db = next;
  }
  return db;
}

func oifits_is_data(db) { return db.__is_data; }
/* DOCUMENT oifits_is_data(db)
     Return whether OI-FITS datablock DB contains data (VIS, VIS2 or T3).
   SEE ALSO: oifits_get_type.
 */

local oifits_get_amporder, oifits_get_amptyp, oifits_get_arrayx, oifits_get_arrayy, oifits_get_arrayz, oifits_get_arrname, oifits_get_calstat, oifits_get_category, oifits_get_corr, oifits_get_corrindx_fluxdata, oifits_get_corrindx_ivis, oifits_get_corrindx_rvis, oifits_get_corrindx_t3amp, oifits_get_corrindx_t3phi, oifits_get_corrindx_vis2data, oifits_get_corrindx_visamp, oifits_get_corrindx_visphi, oifits_get_corrname, oifits_get_date_obs, oifits_get_dec_err, oifits_get_decep0, oifits_get_diameter, oifits_get_equinox, oifits_get_flag, oifits_get_fluxdata, oifits_get_fluxdata_units, oifits_get_fluxerr, oifits_get_fluxerr_units, oifits_get_fov, oifits_get_fovtype, oifits_get_frame, oifits_get_iindx, oifits_get_insname, oifits_get_int_time, oifits_get_ivis, oifits_get_iviserr, oifits_get_jindx, oifits_get_jxx, oifits_get_jxy, oifits_get_jyx, oifits_get_jyy, oifits_get_mjd, oifits_get_mjd_end, oifits_get_mjd_obs, oifits_get_model, oifits_get_ndata, oifits_get_npol, oifits_get_orient, oifits_get_para_err, oifits_get_parallax, oifits_get_phiorder, oifits_get_phityp, oifits_get_pmdec, oifits_get_pmdec_err, oifits_get_pmra, oifits_get_pmra_err, oifits_get_ra_err, oifits_get_raep0, oifits_get_revn, oifits_get_rvis, oifits_get_rviserr, oifits_get_spectyp, oifits_get_sta_index, oifits_get_sta_name, oifits_get_staxyz, oifits_get_sysvel, oifits_get_t3amp, oifits_get_t3amperr, oifits_get_t3phi, oifits_get_t3phierr, oifits_get_target, oifits_get_target_id, oifits_get_tel_name, oifits_get_time, oifits_get_u1coord, oifits_get_u2coord, oifits_get_ucoord, oifits_get_v1coord, oifits_get_v2coord, oifits_get_vcoord, oifits_get_veldef, oifits_get_veltyp, oifits_get_vis2data, oifits_get_vis2err, oifits_get_visamp, oifits_get_visamperr, oifits_get_visphi, oifits_get_visphierr, oifits_get_visrefmap;
local oifits_get_eff_wave, oifits_get_eff_band;
/* DOCUMENT oifits_get_...(master, db);

     This functions query a given attribute of OI-FITS data block DB
     in OI-FITS handle MASTER:

     oifits_get_date_obs  - get UTC start date of observations
     oifits_get_time      - get UTC time of observation (s)
     oifits_get_mjd       - get Modified Julian Date
     oifits_get_int_time  - get integration time (s)
     oifits_get_sta_index - get station numbers contributing to the data
     oifits_get_flag      - get flags
     oifits_get_visamp    - get visibility amplitude
     oifits_get_visamperr - get error in visibility amplitude
     oifits_get_visphi    - get visibility phase (deg)
     oifits_get_visphierr - get error in visibility phase (deg)
     oifits_get_vis2data  - get squared visibility
     oifits_get_vis2err   - get error in squared visibility
     oifits_get_t3amp     - get triple-product amplitude
     oifits_get_t3amperr  - get error in triple product amplitude
     oifits_get_t3phi     - get triple-product phase (deg)
     oifits_get_t3phierr  - get error in triple product phase (deg)
     oifits_get_ucoord    - get u coordinate of the data (m)
     oifits_get_vcoord    - get v coordinate of the data (m)
     oifits_get_u1coord   - get u coordinate of baseline AB of the triangle (m)
     oifits_get_v1coord   - get v coordinate of baseline AB of the triangle (m)
     oifits_get_u2coord   - get u coordinate of baseline BC of the triangle (m)
     oifits_get_v2coord   - get v coordinate of baseline BC of the triangle (m)


     Query attributes for OI_ARRAY data block:

     oifits_get_arrname   - get identifier of corresponding OI_ARRAY
     oifits_get_revn      - get revision number of the table definition
     oifits_get_frame     - get coordinate frame
     oifits_get_arrayx    - get X coordinate of array center (m)
     oifits_get_arrayy    - get Y coordinate of array center (m)
     oifits_get_arrayz    - get Z coordinate of array center (m)
     oifits_get_tel_name  - get telescope name
     oifits_get_sta_name  - get station name
     oifits_get_sta_index - get station number
     oifits_get_diameter  - get element diameter (m)
     oifits_get_staxyz    - get station coordinate relative to array center (m)
     oifits_get_fov       - get size of photometric field of view (arcsec)
     oifits_get_fovtype   - get model for FOV: "FWHM" or "RADIUS"


     Query attributes for OI_TARGET data block:

     oifits_get_target_id - get index number of target(s)
     oifits_get_target    - get target names(s)
     oifits_get_raep0     - get R.A. at mean equinox (deg)
     oifits_get_decep0    - get decl. at mean equinox (deg)
     oifits_get_equinox   - get equinox
     oifits_get_ra_err    - get error in R.A. at mean equinox (deg)
     oifits_get_dec_err   - get error in decl. at mean equinox (deg)
     oifits_get_sysvel    - get Systemic radial velocity (m/s)
     oifits_get_veltyp    - get reference for radial velocity (LSR, GEOCENTR, etc.)
     oifits_get_veldef    - get definition of radial velocity (OPTICAL, RADIO)
     oifits_get_pmra      - get proper motion in R.A. (deg/yr)
     oifits_get_pmdec     - get proper motion in decl. (deg/yr)
     oifits_get_pmra_err  - get error of proper motion in R.A. (deg/yr)
     oifits_get_pmdec_err - get error of proper motion in decl. (deg/yr)
     oifits_get_parallax  - get parallax (deg)
     oifits_get_para_err  - get error in parallax (deg)
     oifits_get_spectyp   - get spectral type
     oifits_get_category  - get category of target: "CAL" or "SCI"


     Query attributes for OI_WAVELENGTH data block:

     oifits_get_insname   - get identifier of corresponding OI_WAVELENGTH
     oifits_get_eff_wave  - get effective wavelength of channel (m)
     oifits_get_eff_band  - get effective bandpass of channel (m)


   SEE ALSO: oifits_new.
*/
func oifits_get_amporder(master, db)          { return db.amporder; }
func oifits_get_amptyp(master, db)            { return db.amptyp; }
func oifits_get_arrayx(master, db)            { return db.arrayx; }
func oifits_get_arrayy(master, db)            { return db.arrayy; }
func oifits_get_arrayz(master, db)            { return db.arrayz; }
func oifits_get_arrname(master, db)           { return db.arrname; }
func oifits_get_calstat(master, db)           { return db.calstat; }
func oifits_get_category(master, db)          { return db.category; }
func oifits_get_corr(master, db)              { return db.corr; }
func oifits_get_corrindx_fluxdata(master, db) { return db.corrindx_fluxdata; }
func oifits_get_corrindx_ivis(master, db)     { return db.corrindx_ivis; }
func oifits_get_corrindx_rvis(master, db)     { return db.corrindx_rvis; }
func oifits_get_corrindx_t3amp(master, db)    { return db.corrindx_t3amp; }
func oifits_get_corrindx_t3phi(master, db)    { return db.corrindx_t3phi; }
func oifits_get_corrindx_vis2data(master, db) { return db.corrindx_vis2data; }
func oifits_get_corrindx_visamp(master, db)   { return db.corrindx_visamp; }
func oifits_get_corrindx_visphi(master, db)   { return db.corrindx_visphi; }
func oifits_get_corrname(master, db)          { return db.corrname; }
func oifits_get_date_obs(master, db)          { return db.date_obs; }
func oifits_get_dec_err(master, db)           { return db.dec_err; }
func oifits_get_decep0(master, db)            { return db.decep0; }
func oifits_get_diameter(master, db)          { return db.diameter; }
func oifits_get_equinox(master, db)           { return db.equinox; }
func oifits_get_flag(master, db)              { return db.flag; }
func oifits_get_fluxdata(master, db)          { return db.fluxdata; }
func oifits_get_fluxdata_units(master, db)    { return db.fluxdata_units; }
func oifits_get_fluxerr(master, db)           { return db.fluxerr; }
func oifits_get_fluxerr_units(master, db)     { return db.fluxerr_units; }
func oifits_get_fov(master, db)               { return db.fov; }
func oifits_get_fovtype(master, db)           { return db.fovtype; }
func oifits_get_frame(master, db)             { return db.frame; }
func oifits_get_iindx(master, db)             { return db.iindx; }
func oifits_get_insname(master, db)           { return db.insname; }
func oifits_get_int_time(master, db)          { return db.int_time; }
func oifits_get_ivis(master, db)              { return db.ivis; }
func oifits_get_iviserr(master, db)           { return db.iviserr; }
func oifits_get_jindx(master, db)             { return db.jindx; }
func oifits_get_jxx(master, db)               { return db.jxx; }
func oifits_get_jxy(master, db)               { return db.jxy; }
func oifits_get_jyx(master, db)               { return db.jyx; }
func oifits_get_jyy(master, db)               { return db.jyy; }
func oifits_get_mjd(master, db)               { return db.mjd; }
func oifits_get_mjd_end(master, db)           { return db.mjd_end; }
func oifits_get_mjd_obs(master, db)           { return db.mjd_obs; }
func oifits_get_model(master, db)             { return db.model; }
func oifits_get_ndata(master, db)             { return db.ndata; }
func oifits_get_npol(master, db)              { return db.npol; }
func oifits_get_orient(master, db)            { return db.orient; }
func oifits_get_para_err(master, db)          { return db.para_err; }
func oifits_get_parallax(master, db)          { return db.parallax; }
func oifits_get_phiorder(master, db)          { return db.phiorder; }
func oifits_get_phityp(master, db)            { return db.phityp; }
func oifits_get_pmdec(master, db)             { return db.pmdec; }
func oifits_get_pmdec_err(master, db)         { return db.pmdec_err; }
func oifits_get_pmra(master, db)              { return db.pmra; }
func oifits_get_pmra_err(master, db)          { return db.pmra_err; }
func oifits_get_ra_err(master, db)            { return db.ra_err; }
func oifits_get_raep0(master, db)             { return db.raep0; }
func oifits_get_revn(master, db)              { return db.revn; }
func oifits_get_rvis(master, db)              { return db.rvis; }
func oifits_get_rviserr(master, db)           { return db.rviserr; }
func oifits_get_spectyp(master, db)           { return db.spectyp; }
func oifits_get_sta_index(master, db)         { return db.sta_index; }
func oifits_get_sta_name(master, db)          { return db.sta_name; }
func oifits_get_staxyz(master, db)            { return db.staxyz; }
func oifits_get_sysvel(master, db)            { return db.sysvel; }
func oifits_get_t3amp(master, db)             { return db.t3amp; }
func oifits_get_t3amperr(master, db)          { return db.t3amperr; }
func oifits_get_t3phi(master, db)             { return db.t3phi; }
func oifits_get_t3phierr(master, db)          { return db.t3phierr; }
func oifits_get_target(master, db)            { return db.target; }
func oifits_get_target_id(master, db)         { return db.target_id; }
func oifits_get_tel_name(master, db)          { return db.tel_name; }
func oifits_get_time(master, db)              { return db.time; }
func oifits_get_u1coord(master, db)           { return db.u1coord; }
func oifits_get_u2coord(master, db)           { return db.u2coord; }
func oifits_get_ucoord(master, db)            { return db.ucoord; }
func oifits_get_v1coord(master, db)           { return db.v1coord; }
func oifits_get_v2coord(master, db)           { return db.v2coord; }
func oifits_get_vcoord(master, db)            { return db.vcoord; }
func oifits_get_veldef(master, db)            { return db.veldef; }
func oifits_get_veltyp(master, db)            { return db.veltyp; }
func oifits_get_vis2data(master, db)          { return db.vis2data; }
func oifits_get_vis2err(master, db)           { return db.vis2err; }
func oifits_get_visamp(master, db)            { return db.visamp; }
func oifits_get_visamperr(master, db)         { return db.visamperr; }
func oifits_get_visphi(master, db)            { return db.visphi; }
func oifits_get_visphierr(master, db)         { return db.visphierr; }
func oifits_get_visrefmap(master, db)         { return db.visrefmap; }

func oifits_get_eff_wave(master, db)
{
  return _oifits_get_link_member(master, db, "__ins", "eff_wave");
}

func oifits_get_eff_band(master, db)
{
  return _oifits_get_link_member(master, db, "__ins", "eff_band");
}

func _oifits_get_link_member(master, db, link, key)
{
  if (! master.__clean) {
    oifits_update, master;
  }
  return master(db(link))(key);
}

func _oifits_push(master, key, value)
{
  h_set, master, key, grow(h_pop(master, key), value);
}

func oifits_list_targets(master)
/* DOCUMENT oifits_list_targets(master);
         or oifits_list_targets, master;

     When called as a function, yields a list of target names found in MiRA
     instance `master`.  When called as a subroutine, print the numerical
     identifier and the name of each target.

   SEE ALSO: oifits_new, oifits_list_instruments, oifits_get_target,
             oifits_get_target_id. */
{
  local result;
  verbose = am_subroutine();
  if (verbose) {
    write, format="%s:\n", "Target(s)";
  }
  for (db = oifits_first(master); db; db = oifits_next(master, db)) {
    type = oifits_get_type(db);
    if (type != OIFITS_TYPE_TARGET) {
      continue;
    }
    names = strtrim(oifits_get_target(master, db), 2);
    if (verbose) {
      ids = oifits_get_target_id(master, db);
      write, format="  %3d: \"%s\"\n", ids, names;
    } else {
      grow, result, names;
    }
  }
  return result;
}

func oifits_list_instruments(master)
/* DOCUMENT oifits_list_instruments(master);
         or oifits_list_instruments, master;

     When called as a function, yields a list of instrument names found in MiRA
     instance `master`.  When called as a subroutine, print the name and the
     spectral range of each instrument.

   SEE ALSO: oifits_new, oifits_list_targets, oifits_get_insname,
             oifits_get_eff_wave. */
{
  local wave, result;
  verbose = am_subroutine();
  if (verbose) {
    write, format="%s:\n", "Instrument(s)";
  }
  for (db = oifits_first(master); db; db = oifits_next(master, db)) {
    type = oifits_get_type(db);
    if (type != OIFITS_TYPE_WAVELENGTH) {
      continue;
    }
    name = strtrim(oifits_get_insname(master, db), 2);
    if (verbose) {
      eq_nocopy, wave, oifits_get_eff_wave(master, db);
      write, format="  \"%s\" (%.1fnm - %.1fnm)\n", name,
        min(wave)/OIFITS_NANOMETER,
        max(wave)/OIFITS_NANOMETER;
    } else {
      grow, result, name;
    }
  }
  return result;
}

/*---------------------------------------------------------------------------*/
/* ERROR MANAGEMENT */

/* _oifits_datablock_builder calls _oifits_error to format the
   error message and which calls _oifits_on_error to manage
   the error.  The result is returned to _oifits_datablock_builder,
   0 means continue (unless it is a fatal error), -1 means return
   immediately.  In this way it is possible to customize the behaviour
   of _oifits_datablock_builder depending on the objectives:

   - reading a file mostly for testing: no error is considered as fatal,
     all errors are stored into an extern variable

   - usual behaviour: immediately report an error.
*/

local _oifits_on_error;
func _oifits_on_error_stop(message)
{
  extern _oifits_error_stack;
  _oifits_error_stack = message;
  return -1;
}

func _oifits_on_error_push(message)
{
  extern _oifits_error_stack;
  grow, _oifits_error_stack, message;
  return 0;
}

func _oifits_on_error_warn(message)
{
  oifits_warn, message;
  return 0;
}

func _oifits_on_error_ignore(message)
{
  return 0;
}

_oifits_on_error = _oifits_on_error_stop;

func _oifits_error(message, ..)
{
  local a1, a2, a3, a4, a5, a6, a7, a8, a9;
  if (more_args()) {
    eq_nocopy, a1, next_arg();
    if (more_args()) {
      eq_nocopy, a2, next_arg();
      if (more_args()) {
        eq_nocopy, a3, next_arg();
        if (more_args()) {
          eq_nocopy, a4, next_arg();
          if (more_args()) {
            eq_nocopy, a5, next_arg();
            if (more_args()) {
              eq_nocopy, a6, next_arg();
              if (more_args()) {
                eq_nocopy, a7, next_arg();
                if (more_args()) {
                  eq_nocopy, a8, next_arg();
                  if (more_args()) {
                    eq_nocopy, a9, next_arg();
                    if (more_args()) {
                      error, "too many arguments";
                    } else {
                      message = swrite(format=message, a1, a2, a3, a4, a5, a6, a7, a8, a9);
                    }
                  } else {
                    message = swrite(format=message, a1, a2, a3, a4, a5, a6, a7, a8);
                  }
                } else {
                  message = swrite(format=message, a1, a2, a3, a4, a5, a6, a7);
                }
              } else {
                message = swrite(format=message, a1, a2, a3, a4, a5, a6);
              }
            } else {
              message = swrite(format=message, a1, a2, a3, a4, a5);
            }
          } else {
            message = swrite(format=message, a1, a2, a3, a4);
          }
        } else {
          message = swrite(format=message, a1, a2, a3);
        }
      } else {
        message = swrite(format=message, a1, a2);
      }
    } else {
      message = swrite(format=message, a1);
    }
  }
  ++_oifits_error_count;
  return _oifits_on_error(message);
}

func _oifits_report_error
{
  extern _oifits_error_stack;
  n = numberof(_oifits_error_stack);
  if (n == 0) {
    message = "no error messages!";
  } else {
    message = _oifits_error_stack(1);
    if (n > 1) {
      message += sum("\n" + _oifits_error_stack(2:n));
    }
  }
  error, message;
}
errs2caller, _oifits_report_error;

func oifits_warn(msg)
/* DOCUMENT oifits_warn, msg;
     Print a warning message.

   SEE ALSO: write, error.
 */
{
  write, format="WARNING: %s\n", msg;
}

func oifits_clear_error(master)
/* DOCUMENT msg = oifits_clear_error(master);
     Clear error messages in MASTER returning the former error messages if any.

   SEE ALSO: oifits_get_error.
 */
{
  return h_pop(master, "__errmsg");
}

func oifits_get_error(master)
/* DOCUMENT msg = oifits_get_error(master);
     Return any error messages in MASTER.

   SEE ALSO: oifits_clear_error.
 */
{
  return h_get(master, "__errmsg");
}

/*---------------------------------------------------------------------------*/
/* SIMULATE DATA */

func oifits_add_noise(master, method, level)
/* DOCUMENT oifits_add_noise(master, method, level);
         or oifits_add_noise, master, method, level;
         or oifits_add_noise(master, method);
         or oifits_add_noise, master, method;

     Add noise to the measurements in MASTER.  The possible methods for that
     are:

       METHOD = 1 or "generate" to add noise to noiseless data; the standard
           deviation of the noise is taken from the contents of MASTER; in this
           case, the LEVEL argument must be nil or omitted.

       METHOD = 2 or "snr" to add noise to noiseless data; the standard
           deviation of noise is computed to achieve a signal-to-noise ratio
           equal to the value of LEVEL.

       METHOD = 3 or "amplify" to add noise to noisy data so that the standard
           deviation of total noise (existing one plus added one) is multiplied
           by the value of LEVEL which must be greater or equal one.  The
           standard deviation of the noise prior to the amplification is taken
           from the contents of MASTER.

   SEE ALSO: oifits_random_normal.
 */
{
  /* Setup noisification method. */
  if (method == 1 || method == "generate") {
    /* Noiseless data: add noise with given standard deviation. */
    if (! is_void(level)) {
      error, "expected nil argument";
    }
    method = 1;
  } else if (method == 2 || method == "snr") {
    /* Noiseless data: add noise to achieve given SNR. */
    if (min(level) <= 0.0) {
      error, "signal to noise ratio must be strictly greater than zero";
    }
    method = 2;
  } else if (method == 3 || method == "amplify") {
    /* Noisy data: amplify noise. */
    if (min(level) < 1.0) {
      error, "noise amplification factor must be greater or equal one";
    }
    method = 3;
  } else {
    error, "bad value for noisification method";
  }

  /* Setup random generator. */
  randgen = oifits_random_normal;

  /* Avoid side effects. */
  if (! am_subroutine()) master = oifits_clone(master, 1n);

  /* Add noise to all data blocks. */
  for (db = oifits_first(master); db ; db = oifits_next(master, db)) {
    type = oifits_get_type(db);
    if (type == OIFITS_TYPE_VIS) {
      local visamp;    eq_nocopy, visamp,    db.visamp;
      local visamperr; eq_nocopy, visamperr, db.visamperr;
      local visphi;    eq_nocopy, visphi,    db.visphi;
      local visphierr; eq_nocopy, visphierr, db.visphierr;
      if (method == 1) {
        /* Noiseless data: add noise with given standard deviation. */
        h_set, db,
          visamp = randgen(visamp, visamperr),
          visphi = randgen(visphi, visphierr);
      } else if (method == 2) {
        /* Noiseless data: add noise to achieve given SNR. */
        visamperr = (1.0/level(1))*visamp;
        visphierr = array((OIFITS_DEGREE/level(0)), dimsof(visphi));
        h_set, db,
          visamp = randgen(visamp, visamperr), visamperr = visamperr,
          visphi = randgen(visphi, visphierr), visphierr = visphierr;
      } else if (method == 3) {
        /* Noisy data: amplify noise. */
        h_set, db,
          visamp = randgen(visamp, sqrt(level(1)^2 - 1.0)*visamperr),
          visamperr = level(1)*visamperr,
          visphi = randgen(visphi, sqrt(level(0)^2 - 1.0)*visphierr),
          visphierr = level(0)*visphierr;
      }
    } else if (type == OIFITS_TYPE_VIS2) {
      local vis2data; eq_nocopy, vis2data, db.vis2data;
      local vis2err;  eq_nocopy, vis2err,  db.vis2err;
      if (method == 1) {
        /* Noiseless data: add noise with given standard deviation. */
        h_set, db, vis2data = randgen(vis2data, vis2err);
      } else if (method == 2) {
        /* Noiseless data: add noise to achieve given SNR. */
        vis2err = (1.0/level(1))*vis2data;
        h_set, db, vis2data = randgen(vis2data, vis2err), vis2err = vis2err;
      } else if (method == 3) {
        /* Noisy data: amplify noise. */
        h_set, db,
          vis2data = randgen(vis2data, sqrt(level(1)^2 - 1.0)*vis2err),
          vis2err = level(1)*vis2err;
      }
    } else if (type == OIFITS_TYPE_T3) {
      local t3amp;    eq_nocopy, t3amp,    db.t3amp;
      local t3amperr; eq_nocopy, t3amperr, db.t3amperr;
      local t3phi;    eq_nocopy, t3phi,    db.t3phi;
      local t3phierr; eq_nocopy, t3phierr, db.t3phierr;
      if (method == 1) {
        /* Noiseless data: add noise with given standard deviation. */
        h_set, db,
          t3amp = randgen(t3amp, t3amperr),
          t3phi = randgen(t3phi, t3phierr);
      } else if (method == 2) {
        /* Noiseless data: add noise to achieve given SNR.  Assume that each
           phase closure is the sum of 3 independent phase data. */
        t3amperr = (1.0/level(1))*t3amp;
        t3phierr = array((sqrt(3.0)*OIFITS_DEGREE/level(0)), dimsof(t3phi));
        h_set, db,
          t3amp = randgen(t3amp, t3amperr), t3amperr = t3amperr,
          t3phi = randgen(t3phi, t3phierr), t3phierr = t3phierr;
      } else if (method == 3) {
        /* Noisy data: amplify noise. */
        h_set, db,
          t3amp = randgen(t3amp, sqrt(level(1)^2 - 1.0)*t3amperr),
          t3amperr = level(1)*t3amperr,
          t3phi = randgen(t3phi, sqrt(level(0)^2 - 1.0)*t3phierr),
          t3phierr = level(0)*t3phierr;
      }
    }
  }
  return master;
}

func oifits_random_normal(a, b)
/* DOCUMENT oifits_random_normal(stdev)
         or oifits_random_normal(mean, stdev)
     Return an array of pseudo random values folling Gaussian distribution
     centerd at MEAN (assumed to be zero in the first case) and standard
     deviation given by STDEV.

   SEE ALSO: random_n.
 */
{
  if (is_void(b)) {
    return random_n(dimsof(a))*a;
  } else {
    return a + random_n(dimsof(a, b))*b;
  }
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func oifits_fix_name(s)
/* DOCUMENT oifits_fix_name(str)
     Return (array of) string(s) STR converted into uppercase letters
     and with trailing spaces stripped.

   SEE ALSO: strtrim, strcase, strupper. */
{
  s = strcase(1, strtrim(s, 2));
  if (is_array((j = where(!s)))) s(j) = ""; // fix nil strings
  return s;
}

func oifits_clone(obj, copy_array)
/* DOCUMENT oifits_clone(obj)
         or oifits_clone(obj, 0/1)
     Returns a clone of object OBJ.  If second argument is true, array
     members found in OBJ are duplicated in its clone; otherwise, such
     members just get referenced in the clone.  The cloning is recursive.

   SEE ALSO: h_new, _lst. */
{
  if (is_array(obj)) {
    if (is_pointer(obj)) {
      write, "FIXME: shall do something special with array of pointers";
    }
    if (copy_array) {
      copy = obj;
      return copy;
    }
    return obj;
  }
  if (is_hash(obj)) {
    clone = h_new();
    for (key = h_first(obj); key; key = h_next(obj, key)) {
      h_set, clone, key, oifits_clone(h_get(obj, key), copy_array);
    }
    return clone;
  }
  if (is_list(obj)) {
    write, "FIXME: shall do something special with lists";
  }
  return obj;
}

func oifits_is_approx(x, y, atol, rtol)
/* DOCUMENT oifits_is_approx(x, y, atol, rtol);

     yields whether X and Y are approximately the same.  Argument ATOL and RTOL
     specify absolute and relative tolerances.  If omitted they are assumed to
     be zero which measn an exact approximation.

   SEE ALSO: allof
 */
{
  if (is_void(atol) || atol <= 0) {
    if (is_void(rtol) || rtol <= 0) {
      return allof(x == y);
    } else {
      return allof(abs(x - y) <= rtol*max(abs(x), abs(y)));
    }
  } else {
    if (is_void(rtol) || rtol <= 0) {
      return (max(abs(x - y)) <= atol);
    } else {
      return allof(abs(x - y) <= atol + rtol*max(abs(x), abs(y)));
    }
  }
}

/*---------------------------------------------------------------------------*/
/* PARSING OF ARGUMENTS */

local _oifits_get_integer_scalar, _oifits_get_integer_vector;
local _oifits_get_integer_array, _oifits_get_real_scalar;
local _oifits_get_real_vector, _oifits_get_real_array;
local _oifits_get_string_scalar, _oifits_get_string_vector;
local _oifits_get_string_array;
/* DOCUMENT _oifits_get_integer_scalar(arg, def)
         or _oifits_get_integer_vector(arg, def)
         or _oifits_get_integer_array(arg, def)
         or _oifits_get_real_scalar(arg, def)
         or _oifits_get_real_vector(arg, def)
         or _oifits_get_real_array(arg, def)
         or _oifits_get_string_scalar(arg, def)
         or _oifits_get_string_vector(arg, def)
         or _oifits_get_string_array(arg, def)

     These functions make ARG into a scalar, a vector or an array of given
     type (long for integer and double for real).  If ARG is nil, its value
     is given by DEF.  Conversion, if required, is done in-place (only
     works if ARG is a symbol not an expression).  These functions return
     0 upon success; -1 otherwise.

   SEE ALSO: is_array, is_scalar, is_vector, structof, oi_check_arg.
 */

func _oifits_get_string_scalar(&arg, def)
{
  if (is_void(arg)) arg = def;
  return (is_scalar(arg) && is_string(arg) ? 0 : -1);
}

func _oifits_get_string_vector(&arg, def)
{
  if (is_void(arg)) arg = def;
  return (is_vector(arg) && is_string(arg) ? 0 : -1);
}

func _oifits_get_string_array(&arg, def)
{
  if (is_void(arg)) arg = def;
  return (is_string(arg) ? 0 : -1);
}

func _oifits_get_integer_scalar(&arg, def)
{
  if (is_void(arg)) arg = def;
  if (is_scalar(arg)) {
    if ((s = structof(arg)) == long) return 0;
    if (s == int || s == short || s == char) {
      arg = long(arg);
      return 0;
    }
  }
  return -1;
}

func _oifits_get_integer_vector(&arg, def)
{
  if (is_void(arg)) arg = def;
  if (is_vector(arg)) {
    if ((s = structof(arg)) == long) return 0;
    if (s == int || s == short || s == char) {
      arg = long(arg);
      return 0;
    }
  }
  return -1;
}

func _oifits_get_integer_array(&arg, def)
{
  if (is_void(arg)) arg = def;
  if (is_array(arg)) {
    if ((s = structof(arg)) == long) return 0;
    if (s == int || s == short || s == char) {
      arg = long(arg);
      return 0;
    }
  }
  return -1;
}

func _oifits_get_real_scalar(&arg, def)
{
  if (is_void(arg)) arg = def;
  if (is_scalar(arg)) {
    if ((s = structof(arg)) == double) return 0;
    if (s == float || s == long || s == int || s == short || s == char) {
      arg = double(arg);
      return 0;
    }
  }
  return -1;
}

func _oifits_get_real_vector(&arg, def)
{
  if (is_void(arg)) arg = def;
  if (is_vector(arg)) {
    if ((s = structof(arg)) == double) return 0;
    if (s == float || s == long || s == int || s == short || s == char) {
      arg = double(arg);
      return 0;
    }
  }
  return -1;
}

func _oifits_get_real_array(&arg, def)
{
  if (is_void(arg)) arg = def;
  if (is_array(arg)) {
    if ((s = structof(arg)) == double) return 0;
    if (s == float || s == long || s == int || s == short || s == char) {
      arg = double(arg);
      return 0;
    }
  }
  return -1;
}

func _oifits_copy_member(dst, src, key)
/**DOCUMENT _oifits_copy_member, dst, src, key;

     Copy value of member KEY from hash-table SRC into hash-table DST.

   SEE ALSO: h_get, h_set.
 */
{
  value = h_get(src, key); /* force a copy for array members */
  h_set, dst, key, value;
}

func _oifits_compare_tables(a, b, atol, rtol)
/** DOCUMENT _oifits_compare_tables(a, b, atol, rtol);

      Compare contents of hash tables A and B.  Return nothing if A and B have
      the same contents or an error message if they differ.  Hash keys
      starting with "__" are ignored.

      Keywords ATOL and RTOL can be used to specify the absolute and relative
      tolerances when comparing numerical values.

   SEE ALSO: oifits_merge, oifits_is_approx.
 */
{
  a_val = [];
  b_val = [];
  keys = h_keys(a);
  nkeys = numberof(keys);
  for (i = 1; i <= nkeys; ++i) {
    key = keys(i);
    if (strpart(key, 1:2) == "__") {
      continue;
    }
    if (! h_has(b, key)) {
      return swrite(format="missing key \"%s\"", key);
    }
    eq_nocopy, a_val, a(key);
    eq_nocopy, b_val, b(key);
    if ((id = identof(a_val)) != identof(b_val)) {
      return swrite(format="different types for key \"%s\"", key);
    }
    if (is_array(a_val)) {
      a_dims = dimsof(a_val);
      b_dims = dimsof(b_val);
      if (numberof(a_dims) != numberof(b_dims) || anyof(a_dims != b_dims)) {
        return swrite(format="different dimensions for key \"%s\"", key);
      }
      if ((id > Y_DOUBLE ? anyof(a_val != b_val) :
           ! oifits_is_approx(a_val, b_val, atol, rtol))) {
        return swrite(format="different value(s) for key \"%s\"", key);
      }
    } else {
      if (a_val != b_val) {
        return swrite(format="different contents for key \"%s\"", key);
      }
    }
  }
  keys = h_keys(b);
  nkeys = numberof(keys);
  for (i = 1; i <= nkeys; ++i) {
    key = keys(i);
    if (strpart(key, 1:2) != "__" && ! h_has(a, key)) {
      return swrite(format="missing key \"%s\"", key);
    }
  }
}

/*---------------------------------------------------------------------------*/
/* OI-FITS FORMAT DESCRIPTION TABLES
 *
 *   The format of the OI-FITS data block is described in what follows by
 *   strings like:
 *
 *     FLAGS KEYWORD FORMAT UNITS COMMENT
 *
 *   where:
 *
 *     FLAGS = 1st bit: optional value (else mandatory)
 *             2nd bit: column (else header keyword)
 *     KEYWORD = keyword for HDU header or column name for table (TTYPE)
 *     FORMAT = nL where n is an integer and L a letter
 *         for the HDU header: n should be 1
 *         for the table: a negative number means an array of rank equals to
 *                        abs(n) with all dimensions equal to NWAVE
 *     UNITS = default units
 *     COMMENT = description
 */

/*-------------------------------------------*/
/* OI_TARGET CLASS DEFINITION (1ST REVISION) */
/*-------------------------------------------*/

_OIFITS_CLASSDEF_TARGET_1 = \
["0 OI_REVN    1I -      revision number of the table definition",
 "2 TARGET_ID  1I -      index number",
 "2 TARGET    16A -      target name",
 "2 RAEP0      1D deg    RA at mean equinox",
 "2 DECEP0     1D deg    DEC at mean equinox",
 "2 EQUINOX    1E yr     equinox",
 "2 RA_ERR     1D deg    error in RA at mean equinox",
 "2 DEC_ERR    1D deg    error in DEC at mean equino",
 "2 SYSVEL     1D m/s    systemic radial velocity",
 "2 VELTYP     8A -      reference for radial velocity",
 "2 VELDEF     8A -      definition of radial velocity",
 "2 PMRA       1D deg/yr proper motion in RA",
 "2 PMDEC      1D deg/yr proper motion in DEC",
 "2 PMRA_ERR   1D deg/yr error of proper motion in RA",
 "2 PMDEC_ERR  1D deg/yr error of proper motion in DEC",
 "2 PARALLAX   1E deg    parallax",
 "2 PARA_ERR   1E deg    error in parallax",
 "2 SPECTYP   16A -      spectral type"];

/*-------------------------------------------*/
/* OI_TARGET CLASS DEFINITION (2ND REVISION) */
/*-------------------------------------------*/

_OIFITS_CLASSDEF_TARGET_2 = \
["0 OI_REVN    1I -      revision number of the table definition",
 "2 TARGET_ID  1I -      index number",
 "2 TARGET    16A -      target name",
 "2 RAEP0      1D deg    RA at mean equinox",
 "2 DECEP0     1D deg    DEC at mean equinox",
 "2 EQUINOX    1E yr     equinox",
 "2 RA_ERR     1D deg    error in RA at mean equinox",
 "2 DEC_ERR    1D deg    error in DEC at mean equino",
 "2 SYSVEL     1D m/s    systemic radial velocity",
 "2 VELTYP     8A -      reference for radial velocity",
 "2 VELDEF     8A -      definition of radial velocity",
 "2 PMRA       1D deg/yr proper motion in RA",
 "2 PMDEC      1D deg/yr proper motion in DEC",
 "2 PMRA_ERR   1D deg/yr error of proper motion in RA",
 "2 PMDEC_ERR  1D deg/yr error of proper motion in DEC",
 "2 PARALLAX   1E deg    parallax",
 "2 PARA_ERR   1E deg    error in parallax",
 "2 SPECTYP   16A -      spectral type",
 "3 CATEGORY   3A -      'CAL' or 'SCI'"];

/*------------------------------------------*/
/* OI_ARRAY CLASS DEFINITION (1ST REVISION) */
/*------------------------------------------*/

_OIFITS_CLASSDEF_ARRAY_1 = \
["0 OI_REVN    1I - revision number of the table definition",
 "0 ARRNAME    1A - array name for cross-referencing",
 "0 FRAME      1A - coordinate frame",
 "0 ARRAYX     1D m array center X-coordinate",
 "0 ARRAYY     1D m array center Y-coordinate",
 "0 ARRAYZ     1D m array center Z-coordinate",
 "2 TEL_NAME  16A - telescope name",
 "2 STA_NAME  16A - station name",
 "2 STA_INDEX  1I - station index",
 "2 DIAMETER   1E m element diameter",
 "2 STAXYZ     3D m station coordinates relative to array center"];

/*------------------------------------------*/
/* OI_ARRAY CLASS DEFINITION (2ND REVISION) */
/*------------------------------------------*/

_OIFITS_CLASSDEF_ARRAY_2 = \
["0 OI_REVN    1I -      revision number of the table definition",
 "0 ARRNAME    1A -      array name for cross-referencing",
 "0 FRAME      1A -      coordinate frame",
 "0 ARRAYX     1D m      array center X-coordinate",
 "0 ARRAYY     1D m      array center Y-coordinate",
 "0 ARRAYZ     1D m      array center Z-coordinate",
 "2 TEL_NAME  16A -      telescope name",
 "2 STA_NAME  16A -      station name",
 "2 STA_INDEX  1I -      station index",
 "2 DIAMETER   1E m      element diameter",
 "2 STAXYZ     3D m      station coordinates relative to array center",
 "2 FOV        1D arcsec photometric field of view",
 "2 FOVTYPE    6A -      model for FOV: 'FWHM' or 'RADIUS'"];

/*-----------------------------------------------*/
/* OI_WAVELENGTH CLASS DEFINITION (1ST REVISION) */
/*-----------------------------------------------*/

_OIFITS_CLASSDEF_WAVELENGTH_1 = \
["0 OI_REVN    1I - revision number of the table definition",
 "0 INSNAME    1A - name of detector for cross-referencing",
 "2 EFF_WAVE   1E m effective wavelength of channel",
 "2 EFF_BAND   1E m effective bandpass of channel"];

/*-----------------------------------------------*/
/* OI_WAVELENGTH CLASS DEFINITION (2ND REVISION) */
/*-----------------------------------------------*/

_OIFITS_CLASSDEF_WAVELENGTH_2 = _OIFITS_CLASSDEF_WAVELENGTH_1;

/*----------------------------------------*/
/* OI_VIS CLASS DEFINITION (1ST REVISION) */
/*----------------------------------------*/

_OIFITS_CLASSDEF_VIS_1 = \
["0 OI_REVN    1I -   revision number of the table definition",
 "0 DATE-OBS   1A -   UTC start date of observations",
 "1 ARRNAME    1A -   name of corresponding OI_ARRAY table",
 "0 INSNAME    1A -   name of corresponding OI_WAVELENGTH table",
 "2 TARGET_ID  1I -   target number as index into OI_TARGET table",
 "2 TIME       1D s   UTC time of observation",
 "2 MJD        1D day modified Julian Day",
 "2 INT_TIME   1D s   integration time",
 "2 VISAMP    -1D -   visibility amplitude",
 "2 VISAMPERR -1D -   error in visibility amplitude",
 "2 VISPHI    -1D deg visibility phase",
 "2 VISPHIERR -1D deg error in visibility phase",
 "2 UCOORD     1D m   U coordinate of the data",
 "2 VCOORD     1D m   V coordinate of the data",
 "2 STA_INDEX  2I -   station numbers contributing to the data",
 "2 FLAG      -1L -   flag"];

/*----------------------------------------*/
/* OI_VIS CLASS DEFINITION (2ND REVISION) */
/*----------------------------------------*/

_OIFITS_CLASSDEF_VIS_2 = \
["0 OI_REVN          1I -   revision number of the table definition",
 "0 DATE-OBS         1A -   UTC start date of observations",
 "0 ARRNAME          1A -   name of corresponding OI_ARRAY table",
 "0 INSNAME          1A -   name of corresponding OI_WAVELENGTH table",
 "1 CORRNAME         1A -   name of corresponding OI_CORR table",
 "1 AMPTYP           1A -   'absolute', 'differential', or 'correlated flux'",
 "1 PHITYP           1A -   'absolute', or 'differential'",
 "1 AMPORDER         1I -   polynomial fit order for differential chromatic amplitudes",
 "1 PHIORDER         1I -   polynomial fit order for differential chromatic phases",
 "2 TARGET_ID        1I -   target number as index into OI_TARGET table",
 "2 TIME             1D s   zero, for backward compatibility",
 "2 MJD              1D day modified Julian Day",
 "2 INT_TIME         1D s   integration time",
 "2 VISAMP          -1D -   visibility amplitude",
 "2 VISAMPERR       -1D -   error in visibility amplitude",
 "3 CORRINDX_VISAMP  1J -   index into correlation matrix for 1st VISAMP element",
 "2 VISPHI          -1D deg visibility phase",
 "2 VISPHIERR       -1D deg error in visibility phase",
 "3 CORRINDX_VISPHI  1J -   index into correlation matrix for 1st VISPHI element",
 "3 VISREFMAP       -2L -   matrix indicating, if true, which spectral channels were taken as reference for differential chromatic visibility computation",
 "3 RVIS            -1D -   real part of complex coherent flux",
 "3 RVISERR         -1D -   error on RVIS",
 "3 CORRINDX_RVIS    1J -   index into correlation matrix for 1st RVIS element",
 "3 IVIS            -1D -   imaginary part of complex coherent flux",
 "3 IVISERR         -1D -   error on IVIS",
 "3 CORRINDX_IVIS    1J -   index into correlation matrix for 1st IVIS element",
 "2 UCOORD           1D m   U coordinate of the data",
 "2 VCOORD           1D m   V coordinate of the data",
 "2 STA_INDEX        2I -   station numbers contributing to the data",
 "2 FLAG            -1L -   flag"];

/*-----------------------------------------*/
/* OI_VIS2 CLASS DEFINITION (1ST REVISION) */
/*-----------------------------------------*/

_OIFITS_CLASSDEF_VIS2_1 = \
["0 OI_REVN    1I -   revision number of the table definition",
 "0 DATE-OBS   1A -   UTC start date of observations",
 "1 ARRNAME    1A -   name of corresponding OI_ARRAY table",
 "0 INSNAME    1A -   name of corresponding OI_WAVELENGTH table",
 "2 TARGET_ID  1I -   target number as index into OI_TARGET table",
 "2 TIME       1D s   UTC time of observation",
 "2 MJD        1D day modified Julian Day",
 "2 INT_TIME   1D s   integration time",
 "2 VIS2DATA  -1D -   squared visibility",
 "2 VIS2ERR   -1D -   error in squared visibility",
 "2 UCOORD     1D m   U coordinate of the data",
 "2 VCOORD     1D m   V coordinate of the data",
 "2 STA_INDEX  2I -   station numbers contributing to the data",
 "2 FLAG      -1L -   flag"];

/*-----------------------------------------*/
/* OI_VIS2 CLASS DEFINITION (2ND REVISION) */
/*-----------------------------------------*/

_OIFITS_CLASSDEF_VIS2_2 = \
["0 OI_REVN            1I -   revision number of the table definition",
 "0 DATE-OBS           1A -   UTC start date of observations",
 "0 ARRNAME            1A -   name of corresponding OI_ARRAY table",
 "0 INSNAME            1A -   name of corresponding OI_WAVELENGTH table",
 "1 CORRNAME           1A -   name of corresponding OI_CORR table",
 "2 TARGET_ID          1I -   target number as index into OI_TARGET table",
 "2 TIME               1D s   zero, for backward compatibility",
 "2 MJD                1D day modified Julian Day",
 "2 INT_TIME           1D s   integration time",
 "2 VIS2DATA          -1D -   squared visibility",
 "2 VIS2ERR           -1D -   error in squared visibility",
 "3 CORRINDX_VIS2DATA  1J -   index into correlation matrix for 1st VIS2DATA element",
 "2 UCOORD             1D m   U coordinate of the data",
 "2 VCOORD             1D m   V coordinate of the data",
 "2 STA_INDEX          2I -   station numbers contributing to the data",
 "2 FLAG              -1L -   flag"];

/*---------------------------------------*/
/* OI_T3 CLASS DEFINITION (1ST REVISION) */
/*---------------------------------------*/

_OIFITS_CLASSDEF_T3_1 = \
["0 OI_REVN    1I -   revision number of the table definition",
 "0 DATE-OBS   1A -   UTC start date of observations",
 "1 ARRNAME    1A -   name of corresponding OI_ARRAY table",
 "0 INSNAME    1A -   name of corresponding OI_WAVELENGTH table",
 "2 TARGET_ID  1I -   target number as index into OI_TARGET table",
 "2 TIME       1D s   UTC time of observation",
 "2 MJD        1D day modified Julian Day",
 "2 INT_TIME   1D s   integration time",
 "2 T3AMP     -1D -   triple product amplitude",
 "2 T3AMPERR  -1D -   error in triple product amplitude",
 "2 T3PHI     -1D deg triple product phase",
 "2 T3PHIERR  -1D deg error in triple product phase",
 "2 U1COORD    1D m   U coordinate of baseline AB of the triangle",
 "2 V1COORD    1D m   V coordinate of baseline AB of the triangle",
 "2 U2COORD    1D m   U coordinate of baseline BC of the triangle",
 "2 V2COORD    1D m   V coordinate of baseline BC of the triangle",
 "2 STA_INDEX  3I -   station numbers contributing to the data",
 "2 FLAG      -1L -   flag"];

/*---------------------------------------*/
/* OI_T3 CLASS DEFINITION (2ND REVISION) */
/*---------------------------------------*/

_OIFITS_CLASSDEF_T3_2 = \
["0 OI_REVN         1I -   revision number of the table definition",
 "0 DATE-OBS        1A -   UTC start date of observations",
 "0 ARRNAME         1A -   name of corresponding OI_ARRAY table",
 "0 INSNAME         1A -   name of corresponding OI_WAVELENGTH table",
 "1 CORRNAME        1A -   name of corresponding OI_CORR table",
 "2 TARGET_ID       1I -   target number as index into OI_TARGET table",
 "2 TIME            1D s   zero, for backward compatibility",
 "2 MJD             1D day modified Julian Day",
 "2 INT_TIME        1D s   integration time",
 "2 T3AMP          -1D -   triple product amplitude",
 "2 T3AMPERR       -1D -   error in triple product amplitude",
 "3 CORRINDX_T3AMP  1J -   index into correlation matrix for 1st T3AMP element",
 "2 T3PHI          -1D deg triple product phase",
 "2 T3PHIERR       -1D deg error in triple product phase",
 "3 CORRINDX_T3PHI  1J -   index into correlation matrix for 1st T3PHI element",
 "2 U1COORD         1D m   U coordinate of baseline AB of the triangle",
 "2 V1COORD         1D m   V coordinate of baseline AB of the triangle",
 "2 U2COORD         1D m   U coordinate of baseline BC of the triangle",
 "2 V2COORD         1D m   V coordinate of baseline BC of the triangle",
 "2 STA_INDEX       3I -   station numbers contributing to the data",
 "2 FLAG           -1L -   flag"];

/*-----------------------------------------*/
/* OI_FLUX CLASS DEFINITION (1ST REVISION) */
/*-----------------------------------------*/

_OIFITS_CLASSDEF_FLUX_1 = \
["0 OI_REVN            1I -      revision number of the table definition",
 "0 DATE-OBS           1A -      UTC start date of observations",
 "0 INSNAME            1A -      name of corresponding OI_WAVELENGTH table",
 "1 ARRNAME            1A -      name of corresponding OI_ARRAY table",
 "1 CORRNAME           1A -      name of corresponding OI_CORR table",
 "1 FOV                1D arcsec area of sky over which flux is integrated",
 "1 FOVTYPE            1A -      model for FOV: 'FWHM' or 'RADIUS'",
 "0 CALSTAT            1A -      'C': spectrum is calibrated, 'U': uncalibrated",
 "2 TARGET_ID          1I -      target number as index into OI_TARGET table",
 "2 MJD                1D day    modified Julian Day",
 "2 INT_TIME           1D s      integration time",
 "2 FLUXDATA          -1D -      flux",
 "2 FLUXERR           -1D -      flux error",
 "3 CORRINDX_FLUXDATA  1J -      index into correlation matrix for 1st FLUXDATA element",
 "3 STA_INDEX          1I -      station number contributing to the data",
 "2 FLAG              -1L -      flag"];

/*-------------------------------------------*/
/* OI_INSPOL CLASS DEFINITION (1ST REVISION) */
/*-------------------------------------------*/

_OIFITS_CLASSDEF_INSPOL_1 = \
["0 OI_REVN    1I -   revision number of the table definition",
 "0 DATE-OBS   1A -   UTC start date of observations",
 "0 NPOL       1I -   number of polarisation types in this table",
 "0 ARRNAME    1A -   name of corresponding OI_ARRAY table",
 "0 ORIENT     1A -   orientation of the Jones matrix",
 "0 MODEL      1A -   method for estimating the Jones matrix",
 "2 TARGET_ID  1I -   target number as index into OI_TARGET table",
 "2 INSNAME    1A -   name of corresponding OI_WAVELENGTH table",
 "2 MJD_OBS    1D day modified Julian Day, start of time lapse",
 "2 MJD_END    1D day modified Julian Day, end of time lapse",
 "2 JXX       -1C -   complex Jones matrix component along X axis",
 "2 JYY       -1C -   complex Jones matrix component along Y axis",
 "2 JXY       -1C -   complex Jones matrix component between X and Y axis",
 "2 JYX       -1C -   complex Jones matrix component between Y and X axis",
 "3 STA_INDEX  1I -   station number for the Jones matrix"];

/*-----------------------------------------*/
/* OI_CORR CLASS DEFINITION (1ST REVISION) */
/*-----------------------------------------*/

_OIFITS_CLASSDEF_CORR_1 = \
["0 OI_REVN   1I - revision number of the table definition",
 "0 CORRNAME  1A - name of correlated data set",
 "0 NDATA     1I - number of correlated data",
 "2 IINDX     1J - first index of correlation matrix element",
 "2 JINDX     1J - second index of correlation matrix element",
 "2 CORR     -1D - matrix element at (IINDX,JINDX)"];

/*---------------------------------------------------------------------------*/
/* INITIALIZATION OF OI-FITS TABLES AND CONSTANTS */

/* Some constants. */
local OIFITS_PI, OIFITS_MICRON;
local OIFITS_DEGREE, OIFITS_ARCSECOND, OIFITS_MILLIARCSECOND;
/* DOCUMENT OIFITS_PI             = 3.1415.....
         or OIFITS_MICRON         = micron to meter conversion factor
         or OIFITS_DEGREE         = degree to radian conversion factor
         or OIFITS_ARCSECOND      = arcsecond to radian conversion factor
         or OIFITS_MILLIARCSECOND = milliarcsecond to radian conversion factor

   SEE ALSO: oifits_get_type.
 */
OIFITS_PI = 3.141592653589793238462643383279503;
OIFITS_DEGREE = OIFITS_PI/180;
OIFITS_ARCSECOND = OIFITS_DEGREE/3600;
OIFITS_MILLIARCSECOND = 1e-3*OIFITS_ARCSECOND;
OIFITS_MICRON = 1e-6;

local OIFITS_TYPE_TARGET, OIFITS_TYPE_WAVELENGTH, OIFITS_TYPE_ARRAY;
local OIFITS_TYPE_VIS, OIFITS_TYPE_VIS2, OIFITS_TYPE_T3, OIFITS_TYPE_FLUX;
local OIFITS_TYPE_INSPOL, OIFITS_TYPE_CORR;
func oifits_get_type(db) { return db.__type; }
/* DOCUMENT oifits_get_type(db)
     Returns OI-FITS type identifier for datablock DB, one of:
       OIFITS_TYPE_TARGET     - for an OI-FITS data block with the list
                                of targets;
       OIFITS_TYPE_WAVELENGTH - for an instrumental  OI-FITS data block;
       OIFITS_TYPE_ARRAY      - for an OI-FITS data block describing the
                                array of telescopes;
       OIFITS_TYPE_VIS        - for an OI-FITS data block with measured
                                complex visibilities;
       OIFITS_TYPE_VIS2       - for an OI-FITS data block with measured
                                squared visibilities;
       OIFITS_TYPE_T3         - for an OI-FITS data block with measured
                                triple products (bispectrum).
       OIFITS_TYPE_FLUX       - for an OI-FITS data block with measured
                                target(s) spectrum.
       OIFITS_TYPE_INSPOL     - for an OI-FITS instrumental polarisation.
       OIFITS_TYPE_CORR       - for an OI-FITS correlation matrix.

   SEE ALSO: oifits_new. */

struct _oifits_classdef {
  string member, keyword, letter, units, comment;
  long multiplier, ctype;
  int optional;
};
local _oifits_classdef;
local _OIFITS_CLASSDEF_FORMAT;
local _OIFITS_CTYPE_LOGICAL, _OIFITS_CTYPE_INTEGER;
local _OIFITS_CTYPE_REAL, _OIFITS_CTYPE_STRING;
func _oifits_classdef_name(class, revn)
{
  if (is_integer(class)) class = _OIFITS_CLASS_NAME_TABLE(class);
  return swrite(format=_OIFITS_CLASSDEF_FORMAT, class, revn);
}
func _oifits_classdef_spec(class, revn)
{
  name = _oifits_classdef_name(class, revn);
  if (! symbol_exists(name)) {
    error, "invalid data-block class or revision number";
  }
  return symbol_def(name);
}
func _oifits_classdef_header(class, revn)
{ return *_oifits_classdef_spec(class, revn)(1); }
func _oifits_classdef_column(class, revn)
{ return *_oifits_classdef_spec(class, revn)(2); }
/* DOCUMENT _oifits_classdef_name(class, revn)
         or _oifits_classdef_column(class, revn)
         or _oifits_classdef_header(class, revn)

     _oifits_classdef_name returns the name of the global variable where
     is stored the class definition for CLASS datablock for revision number
     REVN of the OIFITS format.

     _oifits_classdef_column / _oifits_classdef_header return definition
     table for columns / header keywords of CLASS datablock for revision
     number REVN of the OIFITS format.

   SEE ALSO: */
_OIFITS_CLASSDEF_FORMAT = "_OIFITS_CLASSDEF_%s_%d";
_OIFITS_CTYPE_LOGICAL = 1; /* for format letter 'L' */
_OIFITS_CTYPE_INTEGER = 2; /* for format letters 'I' or 'J' */
_OIFITS_CTYPE_REAL    = 3; /* for format letters 'D' or 'E' */
_OIFITS_CTYPE_STRING  = 4; /* for format letter 'A' */
_OIFITS_CTYPE_COMPLEX = 5; /* for format letter 'C' */

local _OIFITS_REVN_MAX;
local _OIFITS_REVN_DEFAULT;
/* DOCUMENT _OIFITS_REVN_MAX      maximum OI-FITS revision number
            _OIFITS_REVN_DEFAULT  default OI-FITS revision number
   SEE ALSO: */
_OIFITS_REVN_MAX = 2;
_OIFITS_REVN_DEFAULT = 2;

/* Hash table for fast identification of OI-FITS data blocks. */
local _OIFITS_DATABLOCK_CLASS;
local _OIFITS_TYPE_TABLE, _OIFITS_CLASS_NAME_TABLE;
func _oifits_init
{
  extern OIFITS_TYPE_TARGET, OIFITS_TYPE_WAVELENGTH, OIFITS_TYPE_ARRAY;
  extern OIFITS_TYPE_VIS, OIFITS_TYPE_VIS2, OIFITS_TYPE_T3;
  extern OIFITS_TYPE_FLUX, OIFITS_TYPE_INSPOL, OIFITS_TYPE_CORR;
  extern _OIFITS_TYPE_TABLE, _OIFITS_CLASS_NAME_TABLE;
  extern _OIFITS_DATABLOCK_CLASS;

  _OIFITS_DATABLOCK_CLASS = ["TARGET", "WAVELENGTH", "ARRAY", "VIS", "VIS2",
                             "T3", "FLUX", "INSPOL", "CORR"];

  _OIFITS_TYPE_TABLE = h_new();
  _OIFITS_CLASS_NAME_TABLE = array(string, numberof(_OIFITS_DATABLOCK_CLASS));
  for (type = numberof(_OIFITS_DATABLOCK_CLASS); type >= 1; --type) {
    class = _OIFITS_DATABLOCK_CLASS(type);
    h_set, _OIFITS_TYPE_TABLE, "OI_"+class, type;
    _OIFITS_CLASS_NAME_TABLE(type) = class;
    symbol_set, "OIFITS_TYPE_"+class, type;
  }

  /* Character conversion table. */
  convert = array(char, 256);
  convert('_') = '_';
  convert('-') = '_';
  convert('0':'9') = indgen('0':'9');
  convert('a':'z') = indgen('a':'z');
  convert('A':'Z') = indgen('a':'z');

  /* Parse class definition tables. */
  local table;
  letter2code = h_new(l=_OIFITS_CTYPE_LOGICAL,
                      L=_OIFITS_CTYPE_LOGICAL,
                      i=_OIFITS_CTYPE_INTEGER,
                      I=_OIFITS_CTYPE_INTEGER,
                      j=_OIFITS_CTYPE_INTEGER,
                      J=_OIFITS_CTYPE_INTEGER,
                      e=_OIFITS_CTYPE_REAL,
                      E=_OIFITS_CTYPE_REAL,
                      d=_OIFITS_CTYPE_REAL,
                      D=_OIFITS_CTYPE_REAL,
                      a=_OIFITS_CTYPE_STRING,
                      A=_OIFITS_CTYPE_STRING,
                      c=_OIFITS_CTYPE_COMPLEX,
                      C=_OIFITS_CTYPE_COMPLEX); /* fast decoder for
                                                  CTYPE letter */
  format = "%d %s %d%s %s %[^\n]"; /* to decode a single definition line */
  flags = long();
  keyword = string();
  letter = string();
  multiplier = long();
  units = string();
  comment = string();
  for (revn = 1; revn <= _OIFITS_REVN_MAX; ++revn) {
    for (type = numberof(_OIFITS_DATABLOCK_CLASS); type >= 1; --type) {
      class = _OIFITS_DATABLOCK_CLASS(type);
      tablename = _oifits_classdef_name(class, revn);
      if (! symbol_exists(tablename)) {
        continue;
      }
      eq_nocopy, table, symbol_def(tablename);
      if (is_void(table)) continue;
      number = numberof(table);
      cdef = array(_oifits_classdef, number);
      header = array(long, number);
      mdef = h_new(); /* to check for multiple definitions of members */
      kdef = h_new(); /* to check for multiple definitions of keywords */
      for (j = 1; j <= number; ++j) {
        /* Parse header description table. */
        if (sread(table(j), format=format, flags, keyword, multiplier,
                  letter, units, comment) < 5
            || is_void((ctype = letter2code(letter)))
            || strlen(letter) != 1) {
          error, swrite(format="syntax error in line %d of table %s",
                        j, tablename);
        }
        member = strchar(convert(strchar(keyword)(1:-1)));
        if (member == "oi_revn") member = "revn";
        if (strpart(member, -5:0) == "_units") {
          error, swrite(format="illegal member name '%s' at line %d of table %s",
                        member, j, tablename);
        }
        if (h_has(mdef, member)) {
          error, swrite(format="multiple definitions of member '%s' at line %d of table %s",
                        member, j, tablename);
        } else {
          h_set, mdef, member, 1;
        }
        if (h_has(kdef, keyword)) {
          error, swrite(format="multiple definitions of keyword '%s' at line %d of table %s",
                        keyword, j, tablename);
        } else {
          h_set, kdef, keyword, 1;
        }
        in_header = ((flags & 2) == 0);
        if ((in_header && multiplier != 1) ||
            (! in_header && multiplier == 0)) {
          error, swrite(format="illegal multiplier line %d of table %s",
                        j, tablename);
        }
        header(j) = in_header;
        cdef(j) =_oifits_classdef(member = member,
                                  keyword = keyword,
                                  letter = letter,
                                  units = units,
                                  comment = comment,
                                  multiplier = multiplier,
                                  ctype = ctype,
                                  optional = ((flags & 1) != 0));
      }
      symbol_set, tablename, [&cdef(where(header)),
                              &cdef(where(! header))];
    }
  }
}

_oifits_init; /* must be last statement */
_oifits_init = []; /* avoids calling it again */
