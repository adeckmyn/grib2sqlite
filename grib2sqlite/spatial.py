"""Grid geometry, station selection, and interpolation helpers."""

import numpy as np
import pandas
from pyproj import Proj

from .config import logger
from .grib import get_grid_values, get_gridinfo, get_proj4


def read_station_list(station_file):
    """Read a station file (csv).

    Args:
        station_file: name of a csv file

    Returns:
        A pandas table containing the parsed station list
    """
    try:
        station_list = pandas.read_csv(station_file, skipinitialspace=True)
    except OSError:
        logger.error("Can not read station file %s", station_file)
        raise
    return station_list


def get_grid_limits(gid):
    """Find bounding box (lat/lon) of a grid."""
    lons, lats = get_grid_boundary(gid)

    minlat = np.floor(np.min(lats)) - 1
    minlon = np.floor(np.min(lons)) - 1
    maxlat = np.ceil(np.max(lats)) + 1
    maxlon = np.ceil(np.max(lons)) + 1

    minlat = np.max([minlat, -90])
    minlon = np.max([minlon, -180])
    maxlat = np.min([maxlat, 90])
    maxlon = np.min([maxlon, 180])

    return minlon, maxlon, minlat, maxlat


def points_restrict(gid, plist):
    """Restrict the station list to points inside the current domain.

    NOTE: * eccodes returns distances in kilometer
          * store the list for use in next run? run this function seperately?

    Args:
        gid: grib handle
        plist: list of stations

    Returns:
        reduced station list that contains only stations inside the domain
    """
    # 1. Get the bounding box lat/lon values
    #    that is a fast way to eliminate most outside points

    minlon, maxlon, minlat, maxlat = get_grid_limits(gid)
    gridinfo = get_gridinfo(gid)
    nlon = gridinfo["nlon"]
    nlat = gridinfo["nlat"]

    # reduce the table to the bounding box
    # Make a copy! The original retains old row numbers and becomes hard to manage.
    # FIXME: for a global grid, this will delete points that are
    #        "between" the max and min longitude
    #        e.g. if the grid is [0,...,359.95], what about 359.99?
    #        the 0 meridian is pretty important, so we need to fix this
    # NOTE: minlon, maxlon are always according to [-180,180],
    #       even if the grid itself is [0,360] !!!
    if gridinfo["wrap_x"]:
        p1 = plist[(plist["lat"] >= minlat) & (plist["lat"] <= maxlat)].copy()
    else:
        p1 = plist[
            (plist["lat"] >= minlat)
            & (plist["lat"] <= maxlat)
            & (plist["lon"] >= minlon)
            & (plist["lon"] <= maxlon)
        ].copy()

    # 2. Now use grid index (requires projection)
    #    to decide which stations are inside the grid
    #    NOTE: we could skip step 1 and project all stations...

    lon = p1["lon"].tolist()
    lat = p1["lat"].tolist()

    i, j = get_gridindex(lon, lat, gid)
    if gridinfo["wrap_x"]:
        # if x wraps around the globe, don't restrict at all
        # BUT: make sure to take this into account when interpolating!
        # NOTE: Amundsen-Scott SP base has j==0 !
        # With current code, i, j == 0 are OK, but == nlat|nlon-1 needs work
        p2 = p1[(j >= 0) & (j < nlat - 1)].copy()
    else:
        p2 = p1[(i >= 0) & (i < nlon - 1) & (j >= 0) & (j < nlat - 1)].copy()
    return p2


def get_grid_points(gid):
    """Get all lat/lon co-ordinates of the grid points.

    Args:
        gid: GRIB handle
    Returns:
        Numpy arrays with all lat/lon values.
    """
    gridinfo = get_gridinfo(gid)
    p4 = get_proj4(gid)
    proj = Proj(p4)
    nlon = int(gridinfo["nlon"])
    nlat = int(gridinfo["nlat"])
    dx = gridinfo["dx"]
    dy = gridinfo["dy"]

    x0, y0 = proj(gridinfo["lon0"], gridinfo["lat0"])

    xxx = np.empty(nlon)
    yyy = np.empty(nlat)
    for i in range(nlon):
        xxx[i] = x0 + (float(i) * dx)
    for j in range(nlat):
        yyy[j] = y0 + (float(j) * dy)

    x_v, y_v = np.meshgrid(xxx, yyy)
    lons, lats = proj(x_v, y_v, inverse=True)
    # NOTE: the inverse projection of a "wrapped" latlong
    #       turns out to be [-180,180]
    return lons, lats


def get_grid_boundary(gid):
    """Get lat/lon co-ordinates of the grid boundary points.

    Args:
        gid: GRIB handle
    Returns:
        Numpy arrays with lat/lon values.
    """
    gridinfo = get_gridinfo(gid)
    p4 = get_proj4(gid)
    proj = Proj(p4)
    nlon = int(gridinfo["nlon"])
    nlat = int(gridinfo["nlat"])
    dx = gridinfo["dx"]
    dy = gridinfo["dy"]
    # get SW corner
    x0, y0 = proj(gridinfo["lon0"], gridinfo["lat0"])
    xxx = np.fromiter((x0 + i * dx for i in range(nlon)), float)
    yyy = np.fromiter((y0 + i * dy for i in range(nlat)), float)
    x_v = np.concatenate(
        (
            xxx,
            np.full(nlat - 2, xxx[nlon - 1]),
            xxx[::-1],
            np.full(nlat - 2, xxx[0]),
        )
    )
    y_v = np.concatenate(
        (
            np.full(nlon - 1, yyy[0]),
            yyy,
            np.full(nlon - 2, yyy[nlat - 1]),
            yyy[:0:-1],
        )
    )

    lons, lats = proj(x_v, y_v, inverse=True)
    # NOTE: the inverse projection of a "wrapped" latlong
    #       turns out to be [-180,180]
    return lons, lats


def get_gridindex(lon, lat, gid):
    """Convert lat/lon values to (zero-offset) grid indices.

    The SW corner has index (0,0). Points are projected on the grid
    and a (non-integer) index is calculated.

    Args:
        lon: single value or list of longitude values
        lat: single value of list of latitude values
        gid: grib handle

    Returns:
        two vectors with (non-integer) index values.
    """
    gridinfo = get_gridinfo(gid)
    p4 = get_proj4(gid)
    proj = Proj(p4)
    x0, y0 = proj(gridinfo["lon0"], gridinfo["lat0"])
    dx = gridinfo["dx"]
    dy = gridinfo["dy"]

    x, y = proj(np.array(lon), np.array(lat))
    i = (x - x0) / dx
    j = (y - y0) / dy
    return i, j


def train_weights(station_list, gid, lsm=False):
    """Train weights for bilinear and nearest neighbour interpolation.

    We train two kinds of interpolation at once.
    NOTE: Land/Sea Mask is not yet implemented.

    Args:
        station_list: pandas table
        gid: grib handle
        lsm: use Land/Sea mask (ignored for now)

    Returns:
        interpolation weights (nearest neighbour & bilinear)
    """
    # TODO: land/sea mask for T2m...
    if lsm:
        logger.warning("SQLITE: ignoring land/sea mask!")
    # usually, station_list is a pandas table
    # but why not also allow dicts with 'lat' and 'lon' lists
    if isinstance(station_list, dict):
        lat = np.array(station_list["lat"])
        lon = np.array(station_list["lon"])
    else:
        lat = np.array(station_list["lat"].tolist())
        lon = np.array(station_list["lon"].tolist())

    nstations = len(lon)
    i, j = get_gridindex(lon, lat, gid)

    gridinfo = get_gridinfo(gid)
    # eccodes python interface only gives us distances
    # so we need to do some math for bilinear weights
    # assuming the 4 closest points are exactly what we need (OK if dx=dy)
    # NOTE: this assumes dx == dy !!!
    # otherwise, you have to check whether 2nd point is along X or Y axis from first
    # probably easy, just look at the index
    ic = np.round(i).astype("i4")
    jc = np.round(j).astype("i4")
    # FIXME: in rare cases, when i/j are exactly integer, floor==ceil
    # this happens e.g. with Amundsen-Scott South Pole base, j==0.
    # In such cases we should interpolate between only 2 points.
    i0 = np.floor(i).astype("i4")
    i1 = i0 + 1  # np.ceil(i)
    # In the very exceptional case that i0==nlon-1 (*exactly* on the outside border)
    # we may not use i0+1, but since the weights will be zero anyway, we can change to
    # any other value...
    # if i0 == nlon - 1:
    #    if i != i0:
    #        # The points were not correctly constrained to domain!
    j0 = np.floor(j).astype("i4")
    j1 = j0 + 1  # np.ceil(j)
    di = i - i0
    dj = j - j0
    w1 = (1 - di) * (1 - dj)
    w2 = (1 - di) * dj
    w3 = di * (1 - dj)
    w4 = di * dj
    if gridinfo["wrap_x"]:
        logger.info("The domain wraps around the globe.")
        ic[ic == gridinfo["nlon"]] = 0
        ic[ic == -1] = gridinfo["nlon"] - 1
        i1[i1 == gridinfo["nlon"]] = 0
        i0[i0 == -1] = gridinfo["nlon"] - 1

    # Weights per method as arrays of shape (nstations, npoints):
    # grid indices i, j and weights w, such that
    #     value[p] = sum_k data[i[p, k], j[p, k]] * w[p, k]
    weights = {
        "bilin": {
            "i": np.stack([i0, i0, i1, i1], axis=1),
            "j": np.stack([j0, j1, j0, j1], axis=1),
            "w": np.stack([w1, w2, w3, w4], axis=1),
        },
        "nearest": {
            "i": ic[:, np.newaxis],
            "j": jc[:, np.newaxis],
            "w": np.ones((nstations, 1)),
        },
    }
    return weights


def _weight_arrays(method_weights):
    """Return the (i, j, w) arrays of one interpolation method.

    Also accepts the old format (a list per station of [[i, j], w] pairs),
    for weights created outside this package.
    """
    if isinstance(method_weights, dict):
        return method_weights["i"], method_weights["j"], method_weights["w"]
    idx = np.array([[pt[0] for pt in st] for st in method_weights], dtype="i4")
    w = np.array([[pt[1] for pt in st] for st in method_weights], dtype=float)
    return idx[..., 0], idx[..., 1], w


def interp_from_weights(gid, weights, method):
    """Interpolate a GRIB field to points using the given weights.

    Args:
        gid: grib handle
        weights: interpolation weights (see train_weights)
        method: interpolation method (bilin, nearest)

    Returns:
        numpy array of interpolated values (one per station)
    """
    # NOTE: this assumes all records in a file use the same grid representation
    data = get_grid_values(gid)
    i, j, w = _weight_arrays(weights[method])
    return (data[i, j] * w).sum(axis=1)
