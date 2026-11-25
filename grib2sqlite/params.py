"""Parameter lists, matching, caching, and field combination."""

import json
from copy import deepcopy

import numpy as np

from .config import logger
from .grib import get_keylist, get_proj4


def read_param_list(param_file):
    """Read a parameter file (json).

    Args:
        param_file: name of a json file

    Returns:
        A dict containing the parsed parameter list
    """
    try:
        with open(param_file) as pf:
            parameter_list = json.load(pf)
            pf.close()
    except OSError:
        logger.error("Can not read parameter file %s", param_file)
        raise

    return parameter_list


def param_match(gid, parameter_list):
    """Check whether a grib record is in the list of required parameters.

    TODO: can we re-organise the code to avoid getting the same key multiple times?
      But I suspect the impact is minimal

    Args:
        gid: grib handle
        parameter_list: list of parameter descriptors

    Returns:
        Parameter descriptor (from the list) that matches the current grib handle,
        or None.
    """
    for param in parameter_list:
        # param['grib_id'] is a dictionary of keys|values that describe the parameter
        # so you need to check all of these to know if the grib record matches
        plist = param["grib_id"]
        if isinstance(plist, dict):
            plist = [plist]

        for par in plist:
            keylist = list(par.keys())
            ginfo = get_keylist(gid, keylist, "string")
            # we make a deepcopy, because we may have to modify some key values...
            if match_keys(par, ginfo):
                pmatch = deepcopy(param)
                # make sure we have a single grib_id in the pmatch!
                pmatch["grib_id"] = deepcopy(par)
                pmatch["units"] = get_keylist(gid, ["parameterUnits"], "string")[
                    "parameterUnits"
                ]
                # make sure level information is explicitly passed!
                # e.g. for "2t" it may not be in the grib_id list.
                level_type = get_keylist(gid, ["typeOfLevel"], "string")["typeOfLevel"]

                pmatch["typeOfLevel"] = level_type
                pmatch["level"] = int(get_keylist(gid, ["level"], "long")["level"])
                if level_type == "isobaricInhPa":
                    pmatch["level_name"] = "p"
                elif level_type == "heightAboveGround":
                    pmatch["level_name"] = "z"
                elif level_type == "hybridLevel":
                    pmatch["level_name"] = "ml"
                elif level_type == "surface":
                    pmatch["level_name"] = None  #  'level' #  "sfc"
                else:
                    pmatch["level_name"] = None  # "level" #  "xx"

                return pmatch

    # no match was found
    return None


def combine_fields(param, station_list):
    """Combine multiple decoded fields into one parameter.

    We need lat/lon and proj4 values when calculating wind speed & direction.
    Because we should correct for the local rotation of the axes.

    Args:
        param: full parameter descriptor
        station_list: a pandas table with at least lon & lat columns

    Returns:
        a new data matrix that combines the input fields according to parameter descriptor
    """
    nfields = len(param["data"])
    parname = param["harp_param"]
    if "level" in param:
        # NOTE: adding level to the name may sometimes lead to strange results
        #       e.g. S10m10, but it's just for debug messages, so we don't care.
        parname += str(param["level"])

    if any(dd is None for dd in param["data"]):
        # there are missing fields
        # don't make this a "warning", it would show too often...
        logger.info("COMBINE : missing data fields for parameter %s.", parname)
        return None

    logger.debug("COMBINING %i fields for %s", nfields, parname)

    npoints = len(param["data"][0])
    result = np.zeros(npoints)

    if param["function"] == "vector_norm":
        logger.debug("WIND SPEED")
        for ff in range(nfields):
            result += param["data"][ff] * param["data"][ff]

        lat = np.array(station_list["lat"].tolist())
        lon = np.array(station_list["lon"].tolist())
        angle, _mapfactor = rotate_wind(lon, lat, param["proj4"])

        return np.sqrt(result)  # * mapfactor ???

    if param["function"] == "sum":
        logger.debug("SUM")
        for ff in range(nfields):
            result += param["data"][ff]
        return result

    if param["function"] == "vector_angle":
        logger.debug("WIND DIRECTION")
        # NOTE:  make sure we get the wind components in the right order!
        #        They appear in the same order as in config.
        #        So we assume the first entry is "u" (or equivalent).
        param["units"] = "deg"

        # 1. check there are exactly 2 components
        if len(param["data"]) != 2:
            logger.error("ERROR: angle always needs exactly two components.")

        u = param["data"][0]
        v = param["data"][1]
        # 2. Angle is not well defined if both components are too small (no wind).
        # NOTE: numpy knows how to handle u==0. ( np.arctan2(+/-1, 0) == +/- np.pi / 2. )
        # BUT if both u, v -> 0, result depends significantly on sign of u (0 vs 180).
        # FIXME: what should happen when wind speed is (close to) zero?
        #        NaN ? Something else?
        wdir = np.mod(180 + np.arctan2(u, v) * 180.0 / np.pi, 360)
        if param["gridinfo"]["uvRelativeToGrid"] == 1:
            lat = np.array(station_list["lat"].tolist())
            lon = np.array(station_list["lon"].tolist())
            logger.debug("Correcting wind angle.")
            angle, _mapfactor = rotate_wind(lon, lat, param["proj4"])
            # keep the result in [0, 360)
            wdir = np.mod(wdir + angle, 360)

        return wdir

    return None


def rotate_wind(lon, lat, p4):
    """Rotate wind from projection grid to geographical axes for correct wind direction.

    Args:
        lon: longitude (numpy vectors or single value)
        lat: latitude (numpy vectors or single value)
        p4: proj4 definition as a list (not single string)

    Returns:
        angle: correction angle in deg (vector) corresponding to every lat/lon location.
               It is added to the grid-relative wind direction.
    """
    if p4["proj"] == "lcc":
        rad = np.pi / 180.0
        lat1 = p4["lat_1"]
        lat2 = p4.get("lat_2", lat1)
        refcos = np.cos(lat1 * rad)
        refsin = np.sin(lat1 * rad)
        # Cone constant: sin(lat_1) for a tangent cone, the secant formula
        # otherwise (Snyder, eq. 15-3). Negative in the southern hemisphere.
        if np.isclose(lat1, lat2):
            cone = refsin
        else:
            cone = np.log(np.cos(lat1 * rad) / np.cos(lat2 * rad)) / np.log(
                np.tan(np.pi / 4 + lat2 * rad / 2) / np.tan(np.pi / 4 + lat1 * rad / 2)
            )
        # Longitude difference in [-180, 180): GRIB gives LoV in [0, 360),
        # station longitudes are usually in [-180, 180].
        dlon = np.mod(np.asarray(lon) - p4["lon_0"] + 180.0, 360.0) - 180.0
        angle = cone * dlon
        # NOTE: mapfactor still assumes a tangent cone (it is not used).
        mapfactor = np.power(refcos / np.cos(lat * rad), 1 - refcos) * np.power(
            (1 + refsin) / (1 + np.sin(lat * rad)), refsin
        )
    elif p4["proj"] == "latlong":
        angle = np.zeros(len(lat))
        mapfactor = np.ones(len(lat))
    else:
        logger.error("Unimplemented wind rotation for projection %s.", p4["proj"])
        angle = np.zeros(len(lat))
        mapfactor = np.ones(len(lat))

    return angle, mapfactor


def parse_parameter_list(param_list):
    """Parse a (json) structure of required parameters.

    Handle single vs combined fields, model levels etc.

    Args:
        param_list: the list of parameters (from json file)

    Returns:
        param_sgl_list: the list of fields that should be decoded
        param_cmb_list: details for combined fields
    """
    param_sgl_list = []
    param_cmb_list = []

    for param in param_list:
        if isinstance(param["grib_id"], dict):
            # a single field is decoded for this parameter
            if "level" in param["grib_id"] and isinstance(
                param["grib_id"]["level"], list
            ):
                # expand to multiple fields
                for lev in param["grib_id"]["level"]:
                    pid = deepcopy(param)
                    pid["grib_id"]["level"] = str(lev)
                    param_sgl_list.append(pid)
            else:
                pid = deepcopy(param)
                # single field parameters are not cached for now
                # so "data" entry is not needed
                param_sgl_list.append(pid)

        else:
            # the parameter requires multiple fields
            # so we will have to cache the data and calculate later
            # for e.g. wind fields we must also cache grid and projection settings
            nfields = len(param["grib_id"])
            if "common" not in param:
                # the most simple combine case: no common keys
                pid = deepcopy(param)
                pid["data"] = [None] * nfields
                pid["gridinfo"] = None
                pid["proj4"] = None
                param_cmb_list.append(pid)

            elif "level" in param["common"] and isinstance(
                param["common"]["level"], list
            ):
                # there is a common level key, so multiple entries
                for lev in param["common"]["level"]:
                    pid = deepcopy(param)
                    pid["common"]["level"] = str(lev)
                    for f in range(len(pid["grib_id"])):
                        pid["grib_id"][f].update(pid["common"])
                    pid.pop("common")
                    pid["data"] = [None] * nfields
                    pid["gridinfo"] = None
                    pid["proj4"] = None
                    pid["level"] = lev
                    param_cmb_list.append(pid)

            else:
                # there are common keys, but no level expansion
                pid = deepcopy(param)
                for f in range(len(pid["grib_id"])):
                    pid["grib_id"][f].update(pid["common"])
                pid.pop("common")
                pid["data"] = [None] * nfields
                pid["gridinfo"] = None
                pid["proj4"] = None
                param_cmb_list.append(pid)
    return param_sgl_list, param_cmb_list


def match_keys(p1, p2):
    """Match two sets of grib key values.

    Note that a key missing in p1 may be None in p2.

    Args:
      p1: a dictionary with GRIB2 key values
      p2: a dictionary with GRIB2 key values

    Returns:
      True if all (not None) values are equal. False otherwise.
    """
    for k in list(set(list(p1.keys()) + list(p2.keys()))):
        if k not in p1:
            if p2[k] is not None:
                return False
        elif k not in p2:
            if p1[k] is not None:
                return False
        elif p1[k] != p2[k]:
            return False
    return True


def cache_field(
    param, data, param_cmb_list, gid, perturbation_number=None, n_members=None
):
    """Check if a decoded field needs to be cached for later parameters.

    Args:
      param: parameter description of the current data
      data: decode (and interpolated) data
      param_cmb_list: a list (may be MODIFIED!)
      gid: GRIB handle (needed to cache the projection and grid details for every field)
      perturbation_number: ensemble member number, or None for deterministic data
      n_members: ensemble size reported by GRIB, or None for deterministic data

    Returns:
      A count of the number of matching parameters.
      Usually 0 or 1, but wind components may be used for direction and speed,
      so may be matched twice.
    """
    # NOTE: when storing a deterministic run as an ensemble member,
    #       n_members may be "1" while det_member may be > 1.
    #       So often n_members reflect the "expected" number.
    count = 0
    gridinfo_list = ["uvRelativeToGrid"]
    for cmb in param_cmb_list:
        nfields = len(cmb["grib_id"])
        for ff in range(nfields):
            if match_keys(cmb["grib_id"][ff], param["grib_id"]):
                logger.debug("Caching output for %s", cmb["harp_param"])
                if perturbation_number is None:
                    cmb["data"][ff] = np.array(data)
                else:
                    # store each member's component data separately under
                    # cmb["member_data"][perturbation_number][field_index]
                    if "member_data" not in cmb:
                        cmb["member_data"] = {}
                    if perturbation_number not in cmb["member_data"]:
                        cmb["member_data"][perturbation_number] = [None] * nfields
                    cmb["member_data"][perturbation_number][ff] = np.array(data)
                # we assume the unit is the same for all constituents and result
                # this is NOT the case for e.g. wind direction
                # probably should be fixed in the final combination function
                cmb["units"] = param["units"]
                cmb["level"] = param["level"]
                cmb["level_name"] = param["level_name"]
                if n_members is not None:
                    cmb["n_members"] = n_members
                if cmb["gridinfo"] is None:
                    cmb["gridinfo"] = get_keylist(gid, gridinfo_list, "long")
                    cmb["proj4"] = get_proj4(gid)
                count += 1
                continue
    if count > 0:
        logger.debug("Found %i matching combined parameters.", count)
    return count
