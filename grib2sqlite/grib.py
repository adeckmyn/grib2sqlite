"""GRIB metadata and field decoding helpers."""

from datetime import datetime, timezone

import eccodes

from .config import logger


def get_date_info(gid):
    """Forecast time, lead time, accumulation time etc. from GRIB record.

    Args:
        gid: grib handle

    Returns:
        Forecast date/time and lead time as read from the grib record
    """
    keys = [
        "editionNumber",
        "dataDate",
        "hour",
        "minute",
        "second",
    ]
    info = get_keylist(gid, keys, "long")
    if info["editionNumber"] == 1:
        logger.info("GRIB-1 files are experimental.")
        keys2 = [
            "stepUnits",
            "endStep",
        ]
        info2 = get_keylist(gid, keys2, "long")
        info = info | info2

    else:
        # NOTE: using stepUnits as number is bugged (Nov 2024) for sub-hourly
        #       certainly in grib-2, probably also grib-1
        #       but avoid grib-1 for sub-hourly data!
        keys2 = [
            "indicatorOfUnitOfTimeRange",
            "forecastTime",
            "productDefinitionTemplateNumber",
        ]
        info2 = get_keylist(gid, keys2, "long")
        info = info | info2

        # Statistically processed fields (PDT 8, but also ensemble PDT 11, 12,
        # ...) have a time interval. The keys are absent for instantaneous
        # templates, in which case get_keylist returns None.
        keys3 = [
            "indicatorOfUnitForTimeRange",
            "lengthOfTimeRange",
        ]
        info3 = get_keylist(gid, keys3, "long")
        info = info | info3
    return date_from_gribinfo(info)


# Seconds per unit, following GRIB2 code table 4.4 (also used by eccodes for
# the "stepUnits" key in GRIB-1).
_TIME_UNIT_SECONDS = {
    0: 60.0,  # minute
    1: 3600.0,  # hour
    2: 86400.0,  # day
    10: 3 * 3600.0,  # 3 hours
    11: 6 * 3600.0,  # 6 hours
    12: 12 * 3600.0,  # 12 hours
    13: 1.0,  # second
}


def _time_unit_seconds(tunit, keyname):
    """Convert a GRIB code table 4.4 time unit to seconds.

    Raises:
        ValueError: for missing or unsupported units (e.g. month, year).
    """
    try:
        return _TIME_UNIT_SECONDS[tunit]
    except KeyError:
        raise ValueError(f"Unsupported {keyname}: {tunit}") from None


def date_from_gribinfo(info):
    """Interprete list of date information as given by eccodes.

    Args:
        info: a dict with all necessary grib key values

    Returns:
        forecast date (datetime object) and lead time (in seconds!)

    Raises:
        ValueError: if the time unit(s) can not be interpreted.
    """
    fcdate = datetime.strptime(
        "{:8}T{:02}{:02}{:02}".format(
            info["dataDate"], info["hour"], info["minute"], info["second"]
        ),
        "%Y%m%dT%H%M%S",
    ).replace(tzinfo=timezone.utc)

    if info["editionNumber"] == 1:
        tunit = info["stepUnits"]
        leadtime = info["endStep"]
        unit_key = "stepUnits"
    else:
        tunit = info["indicatorOfUnitOfTimeRange"]
        leadtime = info["forecastTime"]
        unit_key = "indicatorOfUnitOfTimeRange"

    leadtime *= _time_unit_seconds(tunit, unit_key)

    # At this point, leadtime may be just the START of accumulation (or
    # min/max/mean) time. Any template with a time interval (PDT 8, 11, 12...)
    # carries lengthOfTimeRange: add it to get the end of the interval.
    if info["editionNumber"] == 2 and info.get("lengthOfTimeRange") is not None:
        tunit2 = info.get("indicatorOfUnitForTimeRange")  # NOTE: "For" vs "Of"
        leadtime += info["lengthOfTimeRange"] * _time_unit_seconds(
            tunit2, "indicatorOfUnitForTimeRange"
        )

    return (fcdate, leadtime)


def get_proj4(gid):
    """Read all projection details and return proj4 string.

    Args:
        gid: GRIB handle

    Returns:
        projection function corresponding to proj4 string
    """
    # NOTE: for now, we assume standard ACCORD Earth shape etc.
    # FIXME: for GRIB-2, we should retrieve the actual earth shape!
    proj4 = {"R": 6371229}
    gridtype = get_keylist(gid, ["gridType"], "string")["gridType"]
    if gridtype in ["lambert", "lambert_lam"]:
        pkeys = ["Latin1InDegrees", "Latin2InDegrees", "LoVInDegrees"]
        p4 = get_keylist(gid, pkeys, "double")

        proj4["proj"] = "lcc"
        proj4["lon_0"] = p4["LoVInDegrees"]
        proj4["lat_1"] = p4["Latin1InDegrees"]
        proj4["lat_2"] = p4["Latin2InDegrees"]
        # TODO: Northern of Southern hemisphere? projectionCentreFlag (?)

    elif gridtype == "polar_stereographic":
        pkeys = ["LoVInDegrees"]
        p4 = get_keylist(gid, pkeys, "double")
        proj4["proj"] = "stere"
        proj4["lon_0"] = p4["LoVInDegrees"]
        proj4["lat_0"] = 90.0  # FIXME: assuming Northern hemisphere!
    elif gridtype == "regular_ll":
        # CHECK
        pkeys = [
            "longitudeOfFirstGridPointInDegrees",
            "longitudeOfLastGridPointInDegrees",
        ]
        p4 = get_keylist(gid, pkeys, "double")
        proj4["proj"] = "latlong"
        proj4["lon_wrap"] = round(sum(p4.values()) / 2.0)

        # use lon_wrap (center longitude) to make sure all values fall in same interval
        # if min_lon = -180: lon_wrap=0 (default)
        # if min_lon = 0: lon_wrap = 180 !!!
        # but you must also be able to handle e.g. [-20,270]

    elif gridtype == "rotated_ll":
        pkeys = [
            "angleOfRotationInDegrees",
            "latitudeOfSouthernPoleInDegrees",
            "longitudeOfSouthernPoleInDegrees",
        ]
        sp_lat = pkeys["latitudeOfSouthernPoleInDegrees"]
        sp_lon = pkeys["longitudeOfSouthernPoleInDegrees"]
        sp_angle = pkeys["angleOfRotationInDegrees"]
        if sp_angle != 0:
            # TODO: I don't know how to handle this. Need examples!
            logger.error("Rotated LatLon with SPangle not supported yet.")
        proj4["proj"] = "ob_tran"
        proj4["o_proj"] = "latlong"
        proj4["o_lat_p"] = -sp_lat
        proj4["o_lon_p"] = 0.0
        proj4["lon_0"] = sp_lon
    #  [rotated] mercator?
    else:
        return None
    return proj4


def get_gridinfo(gid):
    """Read all grid details and return all necessary data.

    Args:
        gid: GRIB handle

    Returns:
        a dictionary with the main grid characteristics
    """
    gridtype = get_keylist(gid, ["gridType"], "string")["gridType"]
    gkeys = [
        "Nx",
        "Ny",
        "uvRelativeToGrid",
        "iScansPositively",
        "jScansPositively",
        "latitudeOfFirstGridPointInDegrees",
        "longitudeOfFirstGridPointInDegrees",
        "latitudeOfLastGridPointInDegrees",
        "longitudeOfLastGridPointInDegrees",
        "DxInMetres",
        "DyInMetres",
        "iDirectionIncrementInDegrees",
        "jDirectionIncrementInDegrees",
    ]
    gg = get_keylist(gid, gkeys, "double")

    result = {
        "nlon": int(gg["Nx"]),
        "nlat": int(gg["Ny"]),
        "dx": gg["DxInMetres"],
        "dy": gg["DyInMetres"],
        "lon0": gg["longitudeOfFirstGridPointInDegrees"],
        "lat0": gg["latitudeOfFirstGridPointInDegrees"],
        "wrap_x": False,
    }
    if gridtype in ["regular_ll", "rotated_ll"]:
        result["dx"] = gg["iDirectionIncrementInDegrees"]
        result["dy"] = gg["jDirectionIncrementInDegrees"]
        logger.debug("%i vs %i", abs(360 - result["dx"] * result["nlon"]), result["dx"])
        if abs(360 - result["dx"] * result["nlon"]) < result["dx"]:
            result["wrap_x"] = True

    if not gg["iScansPositively"]:
        result["lon0"] = gg["longitudeOfLastGridPointInDegrees"]
    if not gg["jScansPositively"]:
        result["lat0"] = gg["latitudeOfLastGridPointInDegrees"]
    return result


def proj4_to_string(proj4):
    """Transfrom a proj4 from dictionary to string.

    Args:
        proj4: a dictionary

    Returns:
        a single string
    """
    result = " ".join([f"+{p}={proj4[p]}" for p in proj4])
    return result


def get_keylist(gid, keylist, ktype="string"):
    """Get list of grib keys. Return "" or None for keys that are not found.

    Args:
        gid: GRIB handle
        keylist: list of key names
        ktype: required output type ("string", "double" or "long")

    Returns:
        A dictionary with all key values. Keys that are not present in the grib record
        return None.
    """
    if ktype == "string":
        func = eccodes.codes_get_string
        miss = ""
    elif ktype == "long":
        func = eccodes.codes_get_long
        miss = None
    elif ktype == "double":
        func = eccodes.codes_get_double
        miss = None

    ginfo = {}
    for kk in keylist:
        try:
            val = func(gid, kk)
        except eccodes.KeyValueNotFoundError:
            val = miss
        ginfo[kk] = val
    return ginfo


def get_ensemble_info(gid):
    """Read ensemble member information from a GRIB record.

    Args:
        gid: GRIB handle

    Returns:
        A tuple (perturbation_number, n_members) where both are None if the
        record does not belong to an ensemble.
    """
    info = get_keylist(
        gid,
        ["perturbationNumber", "numberOfForecastsInEnsemble"],
        "long",
    )
    pert = info["perturbationNumber"]
    n_members = info["numberOfForecastsInEnsemble"]
    # eccodes returns None for missing keys; treat 0 members as non-ensemble too
    if pert is None or n_members is None or n_members == 0:
        return None, None
    return int(pert), int(n_members)


def get_grid_values(gid):
    """Decode a GRIB2 field.

    Args:
      gid: A grib handle

    Returns:
      Numpy array of decoded field.
    """
    gkeys = [
        "Nx",
        "Ny",
        "iScansNegatively",
        "jScansPositively",
        "jPointsAreConsecutive",
        "alternativeRowScanning",
    ]
    ginfo = get_keylist(gid, gkeys, "long")
    nx = ginfo["Nx"]
    ny = ginfo["Ny"]

    # data given by column (C) or by row (Fortran)?
    order = "C" if ginfo["jPointsAreConsecutive"] else "F"
    data = eccodes.codes_get_values(gid).reshape(nx, ny, order=order)
    if ginfo["iScansNegatively"]:
        # logger.warning("Untested data ordering iScansNegatively=1")
        data[range(nx), :] = data[range(nx)[::-1], :]
    if not ginfo["jScansPositively"]:
        # logger.warning("Untested data ordering jScansPositively=0")
        data[:, range(ny)] = data[:, range(ny)[::-1]]

    return data
