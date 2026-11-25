"""GRIB-to-SQLite extraction pipeline."""

from copy import deepcopy

import eccodes
import numpy as np

from .config import (
    default_parameter_list,
    default_sqlite_template,
    default_station_list,
    Grib2SqliteError,
    logger,
)
from .grib import (
    get_date_info,
    get_ensemble_info,
    get_grid_signature,
    get_keylist,
    same_grid,
)
from .params import (
    cache_field,
    combine_fields,
    param_match,
    parse_parameter_list,
    read_param_list,
)
from .spatial import (
    interp_from_weights,
    points_restrict,
    read_station_list,
    train_weights,
)
from .sqlite import create_table, sqlite_name, write_to_sqlite


class _EnsembleAccumulator:
    """Collect interpolated member data for one (param, fcdate, leadtime) slot.

    The accumulator waits until all members have been seen before flushing
    to SQLite.  Member numbers are 1-based (as reported by perturbationNumber).

    Attributes:
        n_members  : total number of members expected
        members    : dict mapping perturbation_number -> data_vector
        param      : parameter descriptor (from the first member seen)
        fcdate     : forecast date
        leadtime   : lead time in seconds
    """

    def __init__(self, n_members, param, fcdate, leadtime):
        self.n_members = n_members
        self.param = deepcopy(param)
        self.fcdate = fcdate
        self.leadtime = leadtime
        self.members = {}  # {perturbation_number: data_vector}

    def add(self, perturbation_number, data_vector):
        self.members[perturbation_number] = data_vector

    @property
    def is_complete(self):
        return len(self.members) == self.n_members

    def missing_members(self):
        if self.is_complete:
            return []
        if not self.members:
            return []
        # Members may be numbered 1..N or 0..(N-1). Use the observed maximum
        # as well as the minimum to infer the numbering: N-1 means 0-based,
        # N means 1-based.
        # If unclear, default to 0-based (like in ENS).
        # This is indicative only: the exact expected set is not known here.
        mx = max(self.members)
        if mx == self.n_members:
            start = 1
        elif 0 in self.members or mx == self.n_members - 1:
            # NOTE: if in reality it is 1-based and member N is missing
            # this will report "0" as missing in stead. Can't be helped.
            start = 0
        else:
            start = 0
        expected = set(range(start, start + self.n_members))
        return sorted(expected - set(self.members))


def _ens_key(param, fcdate, leadtime):
    """Hashable key for the ensemble accumulator dict."""
    level = param.get("level")
    return (param["harp_param"], fcdate, leadtime, level)


def _flush_ensemble(
    accum,
    station_list,
    model_name,
    sqlite_template,
    cleanup=False,
    model_elevation=None,
):
    """Write all accumulated ensemble members for one (param, fcdate, leadtime) to SQLite.

    Each member is stored in a separate column named <model_name>_mbr001, …

    Args:
        accum          : a completed (or partial) _EnsembleAccumulator
        station_list   : pandas table with station info
        model_name     : base model name (no suffix)
        sqlite_template: template string for the output SQLite path
        cleanup        : if True, delete existing rows for the key before upserting
        model_elevation: model orography at the stations, or None
    """
    param = accum.param
    fcdate = accum.fcdate
    leadtime = accum.leadtime
    n_members = accum.n_members

    sqlite_file = sqlite_name(param, fcdate, model_name, sqlite_template)
    # Build one combined pandas table with all member columns in a single pass.
    table = create_table(
        accum.members,
        station_list,
        param,
        fcdate,
        leadtime,
        model_name,
        model_elevation=model_elevation,
    )
    write_to_sqlite(table, sqlite_file, param, model_name, n_members, cleanup=cleanup)

    logger.debug(
        "SQLITE: flushed %i members for %s fcdate=%s leadtime=%s",
        len(accum.members),
        param["harp_param"],
        fcdate,
        leadtime,
    )


# Standard gravity, to convert surface geopotential (m2 s-2) to metres.
GRAVITY = 9.80665
# shortNames accepted as model orography (at typeOfLevel "surface").
OROGRAPHY_NAMES = ("z", "orog", "h", "gh")


def read_model_elevation(clim_file, station_list, weights=None):
    """Read the model orography from a (clim) GRIB file and interpolate it.

    The first record with one of the shortNames in OROGRAPHY_NAMES at
    typeOfLevel "surface" is used. Geopotential (m2 s-2) is converted to m.
    If no weights are given, the station list is restricted to the domain and
    the interpolation weights are trained on this record. So the clim file
    must be on the same grid as the data (see same_grid).

    Args:
        clim_file: name of a GRIB file containing the model orography
        station_list: pandas table with station info
        weights: interpolation weights (if None, they are calculated)

    Returns:
        station_list (possibly restricted), weights, model elevation at the
        stations (numpy array, m) and the grid signature of the clim record.

    Raises:
        Grib2SqliteError: if no usable orography record is found, its units are
                          not recognised or no station lies inside the domain.
    """
    with open(clim_file, "rb") as gfile:
        while True:
            gid = eccodes.codes_grib_new_from_file(gfile)
            if gid is None:
                break
            try:
                info = get_keylist(
                    gid, ["shortName", "typeOfLevel", "parameterUnits"], "string"
                )
                if (
                    info["shortName"] in OROGRAPHY_NAMES
                    and info["typeOfLevel"] == "surface"
                ):
                    return _interp_model_elevation(
                        gid, info, clim_file, station_list, weights
                    )
            finally:
                eccodes.codes_release(gid)
    raise Grib2SqliteError(
        f"No model orography ({'/'.join(OROGRAPHY_NAMES)} at surface) "
        f"found in {clim_file}."
    )


def _interp_model_elevation(gid, info, clim_file, station_list, weights):
    """Interpolate an orography record to the stations (see read_model_elevation)."""
    units = info["parameterUnits"]
    # eccodes spells units differently between versions ("m**2 s**-2" vs "m2 s-2")
    units_norm = units.replace("**", "").replace("^", "").replace(" ", "")
    if units_norm == "m2s-2":
        scale = 1.0 / GRAVITY
    elif units_norm in ("m", "gpm"):
        scale = 1.0
    else:
        raise Grib2SqliteError(
            f"Unexpected units {units!r} for orography {info['shortName']} "
            f"in {clim_file}."
        )

    if weights is None:
        nstation_orig = station_list.shape[0]
        station_list = points_restrict(gid, station_list)
        if station_list.shape[0] == 0:
            raise Grib2SqliteError(f"No stations inside the domain of {clim_file}.")
        logger.info(
            "SQLITE: selected %i stations inside domain from %i.",
            station_list.shape[0],
            nstation_orig,
        )
        logger.debug("SQLITE: training interpolation weights on clim file.")
        weights = train_weights(station_list, gid, lsm=False)

    model_elevation = np.asarray(interp_from_weights(gid, weights, "bilin")) * scale
    logger.info(
        "SQLITE: model elevation from %s (%s, %s).",
        clim_file,
        info["shortName"],
        units,
    )
    return station_list, weights, model_elevation, get_grid_signature(gid)


def parse_grib_file(
    infile,
    param_list=default_parameter_list,
    station_list=default_station_list,
    sqlite_template=default_sqlite_template,
    weights=None,
    model_name="TEST",
    ensemble=False,
    det_member=0,
    cleanup=False,
    clim_file=None,
):
    """Read a GRIB2 file and extract all required data points to SQLite.

    Args:
      infile: input grib2 file
      param_list: list of parameters or name of a json file
      station_list: pandas table with all stations or name of a csv file
      sqlite_template: template for sqlite output files
      weights: interpolation weights (if None, they are calculated)
      model_name: model name (string) used in the SQLite file name and data columns
      ensemble: if True, treat the file as an ensemble GRIB file.
                Member data are stored in columns named <model_name>_mbr001 etc.
                Note that the productDefinitionTemplateNumber and typeOfStatisticalProcessing
                keys differ from deterministic, therefore the param_list must be eps-specific.
      det_member: member number to use when ensemble=True and a record has no
                  perturbationNumber; written as <model_name>_mbrNNN. Default 0.
      cleanup: if True, delete all existing rows for the written keys before
               upserting. By default rows are merged (upsert only).
      clim_file: GRIB file with the model orography (on the same grid as
                 infile), written as model_elevation for T2m. If None,
                 model_elevation is NULL.
    Returns:
        Total number of GRIB records and number of matching parameters found.

    Raises:
        Grib2SqliteError: if the file can not be processed (no stations inside
                          the domain, problems with the clim file, ...).
    """
    if det_member < 0:
        raise ValueError("det_member must be a non-negative integer")

    if isinstance(station_list, str):
        station_list = read_station_list(station_list)
    if isinstance(param_list, str):
        param_list = read_param_list(param_list)

    # Model orography is read first, so it is available for all T2m fields.
    model_elevation = None
    clim_grid = None
    if clim_file is not None:
        station_list, weights, model_elevation, clim_grid = read_model_elevation(
            clim_file, station_list, weights
        )

    # split into "combined" and "direct" parameters
    param_sgl_list, param_cmb_list = parse_parameter_list(param_list)

    logger.info(
        "SQLITE: expecting %i single and %i combined parameters.",
        len(param_sgl_list),
        len(param_cmb_list),
    )

    # For ensemble mode: accumulate member data before writing.
    # Key: _ens_key(param, fcdate, leadtime)  ->  _EnsembleAccumulator
    ens_accumulators = {}

    fcdate = None
    leadtime = None
    gt = 0
    gi = 0
    gd = 0
    gc = 0

    with open(infile, "rb") as gfile:
        while True:
            # loop over all grib records in the file
            gid = eccodes.codes_grib_new_from_file(gfile)
            if gid is None:
                # We've reached the last grib record in the file
                break
            gt += 1
            param = param_match(gid, param_sgl_list)
            if param is not None:
                direct = True
                gd += 1
                logger.debug("SQLITE: found parameter %s", param["harp_param"])
            else:
                # the record is not in our "direct" parameter list
                param = param_match(gid, param_cmb_list)
                if param is None:
                    eccodes.codes_release(gid)
                    continue
                # the record is in the "combined" list
                direct = False
                gc += 1
                logger.debug(
                    "SQLITE: found combinded parameter %s:%s",
                    param["harp_param"],
                    param["grib_id"],
                )

            # we have a matching parameter
            gi += 1
            # NOTE: actually, we only need to get fcdate & lead time once
            #       and we could even consider those to be known
            try:
                fcdate, leadtime = get_date_info(gid)
            except ValueError as err:
                logger.error(
                    "SQLITE: skipping record %i (%s): %s",
                    gt,
                    param["harp_param"],
                    err,
                )
                eccodes.codes_release(gid)
                continue
            logger.debug(
                "SQLITE: gt=%i gi=%i fcdate=%s, leadtime=%i, direct=%s",
                gt,
                gi,
                fcdate.isoformat(),
                leadtime,
                direct,
            )

            if clim_grid is not None:
                # The orography and weights come from the clim file: check
                # (once) that the data are on the same grid.
                if not same_grid(clim_grid, get_grid_signature(gid)):
                    eccodes.codes_release(gid)
                    raise Grib2SqliteError(
                        f"The grid of {infile} differs from the grid of the "
                        f"clim file {clim_file}."
                    )
                clim_grid = None

            if weights is None:
                # We assume that station list and weights are the same for all GRIB fields
                # so we only "train" once
                # First reduce the station list to points inside the domain
                nstation_orig = station_list.shape[0]
                station_list = points_restrict(gid, station_list)
                if station_list.shape[0] == 0:
                    # In this case, we can not extract any points
                    eccodes.codes_release(gid)
                    raise Grib2SqliteError(
                        f"No stations inside the model domain of {infile}."
                    )

                logger.info(
                    "SQLITE: selected %i stations inside domain from %i.",
                    station_list.shape[0],
                    nstation_orig,
                )
                # create a list of interpolation weights
                logger.debug("SQLITE: training interpolation weights.")
                weights = train_weights(station_list, gid, lsm=False)

            # by default, we do bilinear interpolation
            method = param.get("method", "bilin")
            # interpolate to station points
            data_vector = interp_from_weights(gid, weights, method)

            if ensemble:
                # NOTE: ensemble and deterministic fields occur in the same file
                #       we could get inconsistent behaviour because
                #       n_members would change randomly between 1 and N
                #       But this could be fixed with good parameter lists
                #       If template is used, ens and det can not be
                #       the same parameter, so it should be OK.
                perturbation_number, n_members = get_ensemble_info(gid)

                if perturbation_number is None:
                    logger.debug(
                        "SQLITE: ensemble=True but record has no perturbationNumber "
                        "for parameter %s – setting to mbr%03d.",
                        param["harp_param"],
                        det_member,
                    )
                    perturbation_number = det_member
                    n_members = 1
                    #eccodes.codes_release(gid)
                    #continue
            else:
                perturbation_number = None
                n_members = None

            # Cache data for combined parameters (wind speed/direction etc.)
            # Each member's components are stored separately so that
            # combine_fields can be called independently per member.
            cache_field(
                param, data_vector, param_cmb_list, gid, perturbation_number, n_members
            )

            if direct:
                if ensemble:
                    key = _ens_key(param, fcdate, leadtime)
                    if key not in ens_accumulators:
                        ens_accumulators[key] = _EnsembleAccumulator(
                            n_members, param, fcdate, leadtime
                        )
                    ens_accumulators[key].add(perturbation_number, data_vector)

                    if ens_accumulators[key].is_complete:
                        accum = ens_accumulators.pop(key)
                        _flush_ensemble(
                            accum,
                            station_list,
                            model_name,
                            sqlite_template,
                            cleanup=cleanup,
                            model_elevation=model_elevation,
                        )
                else:
                    # this is a "direct" field, create a table and write to SQLite
                    sqlite_file = sqlite_name(
                        param, fcdate, model_name, sqlite_template
                    )
                    data = create_table(
                        data_vector,
                        station_list,
                        param,
                        fcdate,
                        leadtime,
                        model_name,
                        model_elevation=model_elevation,
                    )
                    write_to_sqlite(
                        data,
                        sqlite_file,
                        param,
                        model_name,
                        None,
                        cleanup=cleanup,
                    )

            eccodes.codes_release(gid)

    # at this point, all grib records have been parsed
    # Flush any accumulators that never reached n_members (missing members at
    # end of file).  This is normal for split ensemble files where the control
    # member (e.g. member 51) lives in a separate file.
    for key, accum in ens_accumulators.items():
        missing = accum.missing_members()
        logger.info(
            "SQLITE: partial ensemble flush for %s fcdate=%s leadtime=%s – "
            "members not in this file: %s.",
            accum.param["harp_param"],
            accum.fcdate,
            accum.leadtime,
            missing,
        )
        _flush_ensemble(
            accum,
            station_list,
            model_name,
            sqlite_template,
            cleanup=cleanup,
            model_elevation=model_elevation,
        )

    # Check all combined parameters and write to SQLite
    logger.debug("SQLITE: checking cached combined fields")
    for param in param_cmb_list:
        if ensemble:
            # In ensemble mode each member's components were stored separately
            # by cache_field.  Combine and accumulate per member, then
            # flush exactly like a direct ensemble parameter.
            member_data = param.get("member_data", {})
            if not member_data:
                logger.info(
                    "COMBINE: no ensemble member data cached for %s.",
                    param["harp_param"],
                )
                continue

            # Prefer the ensemble size reported by GRIB; fall back to the
            # highest perturbation number when that key was not available.
            n_members = param.get("n_members") or max(member_data.keys())

            combined = {}
            for mbr_num, component_data in member_data.items():
                # Temporarily swap in this member's component data so that
                # combine_fields (which reads param["data"]) works unchanged.
                # TODO: we should modify combine_fields to accept an ensemble data argument.

                param["data"] = component_data
                data_vector = combine_fields(param, station_list)
                if data_vector is not None:
                    combined[mbr_num] = data_vector
                else:
                    logger.info(
                        "COMBINE: missing component for %s member %i – skipping member.",
                        param["harp_param"],
                        mbr_num,
                    )

            if combined:
                # Create the accumulator only now: it takes a copy of param,
                # and combine_fields may have modified it (e.g. units="deg"
                # for wind direction).
                cmb_accum = _EnsembleAccumulator(n_members, param, fcdate, leadtime)
                for mbr_num, data_vector in combined.items():
                    cmb_accum.add(mbr_num, data_vector)
                missing = cmb_accum.missing_members()
                if missing:
                    logger.info(
                        "SQLITE: partial combined ensemble flush for %s – "
                        "members not in this file: %s.",
                        param["harp_param"],
                        missing,
                    )
                _flush_ensemble(
                    cmb_accum,
                    station_list,
                    model_name,
                    sqlite_template,
                    cleanup=cleanup,
                    model_elevation=model_elevation,
                )
        else:
            # Deterministic: single combine_fields call, write directly.
            data_vector = combine_fields(param, station_list)
            # only write to SQLite if ALL components were found!
            if data_vector is not None:
                logger.debug("SQLITE: writing combined field")
                sqlite_file = sqlite_name(param, fcdate, model_name, sqlite_template)
                # Pass the base model name: create_table()/mbr_column() append
                # the "_det" suffix once for deterministic fields.
                data = create_table(
                    data_vector,
                    station_list,
                    param,
                    fcdate,
                    leadtime,
                    model_name,
                    model_elevation=model_elevation,
                )
                write_to_sqlite(data, sqlite_file, param, model_name, cleanup=cleanup)

    logger.info(
        "SQLITE: Total %i records. %i matching of which %i direct and %i combined.",
        gt,
        gi,
        gd,
        gc,
    )
    # Return: total count and # of matching param
    return gt, gi
