# grib2sqlite

Interpolate grib files to sqlite tables for the harp verification tool.
Can be called via command line or used as a module.

## Command line

    python3 -m grib2sqlite <options> <grib file(s)>

usage: grib2sqlite [-h] [-p param_file] [-s station_file] [-m model_name]
                   [-t output_template] [-o output_path] [-d debug_level]
                   /path/to/file(s) [/path/to/file(s) ...]

### positional arguments:

    /path/to/file(s)
    Wild cards are allowed.

### options:

*  -h, --help          show help message and exit
*  -p <param_file>       parameter list file (json)
*  -s <station_file>     station list file (csv)
*  -m <model_name>       base model name used in SQLite (deterministic fields get a _det suffix, ensemble members get _mbrNNN)
*  -t <output_template>  template used for SQLite files
*  -o <output_path>      path used for SQLite files
*  -d <debug_level>      verbosity level (0...3) default: 1
*  -e, --ensemble        treat input as ensemble data (members in _mbrNNN columns)
*  --det-member <N>      with -e, member number used for records without a perturbationNumber (default: 0, i.e. _mbr000)
*  --cleanup            delete existing rows for the written keys before upserting (default: merge/upsert only)

Option `-d0` surpresses all messages except errors, `-d1` adds warnings, and so on up to `-d3` that will output a large amount of debugging information.

## Ensemble mode

With `-e`, records that have a `perturbationNumber` are written to `<model>_mbrNNN` columns. Records without one (deterministic fields) are written to member 0 by default. Use `--det-member <N>` to choose another member number, e.g. when each member forecast is supplied as a separate deterministic GRIB file:

    python3 -m grib2sqlite -e --det-member 7 -m IFSENS member007.grib2
    python3 -m grib2sqlite -e --det-member 8 -m IFSENS member008.grib2

Both runs must use the same model name, output template/path and parameter list. Rows are merged by upsert, so the second run adds its member column without erasing the first. Do **not** combine this with `--cleanup`: it deletes the whole row for a key before writing, which would remove member columns written by earlier runs.

Note that the productDefinitionTemplateNumber and typeOfStatisticalProcessing keys of ensemble data differ from deterministic data, so an ensemble-specific parameter list is usually needed.

## Station list

The station list shoud be given as a .csv file. As an example, there is a default list in the data sub-directory containing a global set of SYNOP stations.

The .csv should have (at least) the following columns:

* **SID**: Station ID (e.g. WMO code)
* **lat**, **lon**: latitude and longitude
* **elev**: elevation above sea level (optional, currently not used)
* **name**

## Parameter list

The parameters are prescribed in a .json file. There is a default file that shows the main structures available. The json file should contain a (single) list of parameter descriptors explained below.

```
"SID","lat","lon","elev","name"
1001,70.9331,-8.6667,9.4,"JAN MAYEN"
1002,80.0592,16.25,8,"VERLEGENHUKEN"
```

### Simple parameters

Basic variables that can be read from a single GRIB record, can be described by a dict with the following three elements

* **harp_param**: the "harp" name for the variable
* **method**: interpolation method (bilin or nearest)
* **grib_id**: a list of key values that uniquely define the GRIB record.

Examples:

```
[
{ "harp_param":"T2m",
  "method":"bilin",
  "grib_id":{ "shortName":"2t", "productDefinitionTemplateNumber":"0"}
},
{ "harp_param":"MAXT2m",
  "method":"bilin",
  "grib_id":{ "shortName":"2t", "productDefinitionTemplateNumber":"8", "typeOfStatisticalProcessing":"max"}
},
]
```

Note that in the first example the template number is added because Max/min T2m would have the same shortName, but different productDefinitionTemplateNumber (8). Without this specification the description would not be unique and the extracted value would depend on the order that the fields were found in the file.

### Multiple levels

The key values can also be lists. In this case, the fields will all be extracted. Typically, this would be for levels (pressure, height, model level...).

Examples:

```
{ "harp_param":"RH",
  "method":"bilin",
  "grib_id":{"shortName":"r", "productDefinitionTemplateNumber":"0",
             "typeOfLevel":"isobaricInhPa",
             "level":[ 50, 100, 250, 300, 600, 700, 800, 925, 950]}
},
```

### Combined fields

Some parameters are not available as a single GRIB record, but must be constructed from various fields. For instance wind speed from (u,v) components or total precipitation from sub-types.
For this, **grib_id** can contain the different fields in a list (**NOTE:** in this case the keys themselves form a list, not the values!). There are also two additional entries to the dict:

- **common**: the GRIB keys that the various fields should have in common, like pressure level etc.
- **function**: the method to combine the various fields. Currently, the possibilities are 
    - **vector_norm** (wind speed),
    - **vector_angle** (wind direction),
    - **sum** (precipitation sub-types).  

Future additions may include elevation-correction for T2m, dew point temperature (Magnus formula), ...

**Limitation:** combined fields assume that each GRIB file contains a single forecast date and lead time. The component fields are cached (one slot per parameter and level) and only combined after the whole file has been read, using the forecast date and lead time of the last matching record in the file. If a file contains several lead times, only one combined value per parameter is written, and it may mix components from different lead times or be stored under the wrong lead time. No warning is given. Single (direct) parameters are not affected. Until this is fixed, split such files into one file per lead time before processing (e.g. with `grib_copy -w step=6 in.grib2 out.grib2`, or `grib_copy in.grib2 'out_[step].grib2'`).

Examples:

```
{ "harp_param":"Pcp",
  "method":"nearest",
  "grib_id":[ {"parameterNumber":"65"},
              {"parameterNumber":"66"},
              {"parameterNumber":"75"}
  ],
  "common":{ "productDefinitionTemplateNumber":"8",
             "typeOfStatisticalProcessing":"accum",
             "parameterCategory":"1"
  },
  "function":"sum"
},
```



