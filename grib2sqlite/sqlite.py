"""SQLite naming, schema, and write helpers."""

import os
import re
import sqlite3
import threading
from contextlib import closing

from .config import logger


# Per-SQLite-file locks used by parallel workers (threads) to serialise the
# file-creation step and avoid "duplicate column name" races.
# Keys must be normalised.

# NOTE: the per-file locks are only valid for threads inside this process.
#       we use WAL for safe readign durign writes, but if two processes write
#       at the same time and the file system does not coprrectly report locking,
#       data corruption is still possible.
_sqlite_file_locks: dict = {}
_sqlite_file_locks_lock = threading.Lock()


def _get_file_lock(path: str) -> threading.Lock:
    """Return a per-file threading.Lock, creating it on first use.

    All workers that share the same process (e.g. ThreadPoolExecutor) will
    receive the same lock object for the same absolute path, serialising the
    SQLite file-creation and schema-patching steps without blocking unrelated
    files.

    Args:
        path: absolute or relative path to the SQLite file (normalised here)

    Returns:
        A threading.Lock unique to this resolved path.
    """
    key = os.path.realpath(os.path.abspath(os.path.normpath(path)))
    with _sqlite_file_locks_lock:
        if key not in _sqlite_file_locks:
            _sqlite_file_locks[key] = threading.Lock()
        return _sqlite_file_locks[key]

def _configure_connection(con):
    """Apply the standard pragmas to a freshly opened SQLite connection.

    busy_timeout is set before journal_mode so the mode change and DDL can
    also wait for a lock held by another process. WAL may be unsupported on
    network filesystems; warn and continue rather than break single-process
    use, but make the degraded state visible.
    """
    con.execute("PRAGMA busy_timeout = 30000")
    journal_mode = con.execute("PRAGMA journal_mode=WAL").fetchone()[0]
    if str(journal_mode).lower() != "wal":
        logger.warning(
            "SQLITE: WAL mode unavailable (journal_mode=%s); concurrent "
            "writers may fail with SQLITE_BUSY.",
            journal_mode,
        )

def _validate_path_component(value, label):
    """Reject path separators and traversal in template substitutions.

    Args:
        value: value that will be inserted into the SQLite file path
        label: name used in the error message

    Returns:
        The validated value.
    """
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a non-empty string")
    if "\x00" in value:
        raise ValueError(f"{label} contains a NUL character")
    if value in (".", "..") or "/" in value or "\\" in value:
        raise ValueError(
            f"{label} must not contain path separators or be '.'/'..': {value!r}"
        )
    return value


def sqlite_name(param, fcdate, model_name, sqlite_template):
    """Create the full name of the SQLite file from template and date.

    Args:
        param: parameter descriptor
        fcdate: forecast date
        model_name: model name (string) used in the SQLite file name and data columns
        sqlite_template: template for SQLite file name

    Returns:
        full name of SQLite file
    """
    # The model name and parameter name are inserted into the file path: make
    # sure they cannot escape the output directory.
    if "{PP}" in sqlite_template:
        _validate_path_component(param["harp_param"], "harp_param")
    if "{MODEL}" in sqlite_template:
        _validate_path_component(model_name, "model_name")
    result = sqlite_template
    result = result.replace("{PP}", param["harp_param"])
    result = result.replace("{YYYY}", fcdate.strftime("%Y"))
    result = result.replace("{MM}", fcdate.strftime("%m"))
    result = result.replace("{DD}", fcdate.strftime("%d"))
    result = result.replace("{HH}", fcdate.strftime("%H"))
    result = result.replace("{MODEL}", model_name)
    return result


def ensure_table_columns(con, member_columns):
    """Add any missing _mbrNNN columns to an existing FC table.

    This is needed when the SQLite file was created by an earlier run
    with a different (smaller) member set.

    Args:
        con            : open sqlite3 connection
        member_columns : iterable of column names
    """
    # PRAGMA table_info returns rows of (cid, name, type, notnull, dflt_value, pk)

    with con, closing(con.cursor()) as cur:
        existing = {row[1] for row in cur.execute("PRAGMA table_info(FC)").fetchall()}
        for col in sorted(member_columns):
            if col not in existing:
                cur.execute(
                    f"ALTER TABLE FC ADD COLUMN {_quote_identifier(col)} DOUBLE"
                )


def mbr_column(model_name, perturbation_number):
    """Build the SQLite column name for a given ensemble member.

    Args:
        model_name: base model name (without any suffix)
        perturbation_number: integer member index (1-based)

    Returns:
        column name string, e.g. "AROME_mbr003"
    """
    if perturbation_number is None:
        return f"{model_name}_det"
    else:
        return f"{model_name}_mbr{perturbation_number:03d}"


def list_mbr_columns(table, model_name):
    """List all model data columns in a pandas table.

    These are the columns produced by mbr_column(): <model_name>_det and
    <model_name>_mbrNNN. An exact pattern is used (not a substring match), so
    that e.g. a model name "t" does not select "lat" or "units".

    Args:
        table: pandas DataFrame
        model_name: base model name (without any suffix)

    Returns:
        list of column names
    """
    pattern = re.compile(re.escape(model_name) + r"_(det|mbr\d+)")
    return [col for col in table.columns if pattern.fullmatch(str(col))]


def create_table(data, station_list, param, fcdate, leadtime, model_name):
    """Put a data vector (interpolated field) to pandas table.

    Args:
      data: numpy vector or a dict of such vectors
      station_list: pandas table with station list
      param: parameter descriptor
      fcdate: forecast date (datetime object)
      leadtime: lead time (int)
      model_name: model name (string)

    Returns:
      pandas table
    """
    station_keys = ["SID", "lat", "lon"]
    table = station_list[station_keys].copy()

    # prepare SQLITE output
    fcd = int(fcdate.timestamp())
    vad = int(fcdate.timestamp() + leadtime)

    table["fcst_dttm"] = fcd
    # NOTE: currently FCTABLE still expects leadtime in hours!
    table["lead_time"] = leadtime / 3600.0
    table["valid_dttm"] = vad
    table["parameter"] = param["harp_param"]
    table["units"] = param["units"]
    logger.debug("SQLITE: writing parameter %s", param["harp_param"])
    # model_elevation only relevant for T2m
    # NOTE: we can not yet read model elevation from clim file!
    # otherwise remove column
    if "T2m" in param["harp_param"]:
        if "model_elevation" not in table.columns.to_numpy().tolist():
            table["model_elevation"] = 0
        if "elev" not in table.columns.to_numpy().tolist():
            table["elev"] = 0
    else:
        if "model_elevation" in table.columns.to_numpy().tolist():
            table = table.drop("model_elevation", axis=1)
        if "elev" in table.columns.to_numpy().tolist():
            table = table.drop("elev", axis=1)

    if param["level"] is not None and param["level_name"] is not None:
        table[param["level_name"]] = int(param["level"])
    # the actual data column(s)
    # note that we want to order the members in the final table!
    if isinstance(data, dict):
        for mbr, data_vector in sorted(data.items()):
            col_name = mbr_column(model_name, mbr)
            table[col_name] = data_vector
    else:
        col_name = mbr_column(model_name, None)
        table[col_name] = data

    return table


def _quote_identifier(name):
    """Quote an SQLite identifier (column or table name)."""
    name = str(name)
    if "\x00" in name:
        raise ValueError("SQL identifier contains a NUL character")
    return '"' + name.replace('"', '""') + '"'


def _upsert_conflict_columns(con, table):
    """Find a unique FC index that can be used as an upsert conflict target.

    Args:
        con: database connection
        table: pandas table about to be written

    Returns:
        List of unique-index column names, or None if no suitable index exists.
    """
    table_columns = set(table.columns)
    for row in con.execute("PRAGMA index_list('FC')").fetchall():
        index_name = row[1]
        is_unique = row[2]
        if not is_unique:
            continue
        index_columns = [
            info[2]
            for info in con.execute(
                "SELECT * FROM pragma_index_info(?)", (index_name,)
            ).fetchall()
        ]
        if index_columns and set(index_columns) <= table_columns:
            return index_columns
    return None


def _upsert_table(con, table, conflict_columns):
    """Insert rows, updating existing rows on a unique-index conflict.

    Only columns present in the table are updated, so member columns that are
    not part of this write keep their existing values.

    Args:
        con: database connection
        table: pandas table to write
        conflict_columns: columns of the FC unique index used as conflict target
    """
    columns = list(table.columns)
    column_list = ", ".join(_quote_identifier(col) for col in columns)
    placeholders = ", ".join("?" for _ in columns)
    conflict = ", ".join(_quote_identifier(col) for col in conflict_columns)
    key_columns = set(conflict_columns)
    update_columns = [col for col in columns if col not in key_columns]

    if update_columns:
        updates = ", ".join(
            f"{_quote_identifier(col)} = excluded.{_quote_identifier(col)}"
            for col in update_columns
        )
        action = f"DO UPDATE SET {updates}"
    else:
        action = "DO NOTHING"

    sql = (
        f"INSERT INTO FC ({column_list}) VALUES ({placeholders}) "
        f"ON CONFLICT ({conflict}) {action}"
    )
    con.executemany(sql, table.itertuples(index=False, name=None))


def write_to_sqlite(
    table, sqlite_file, param, model_name, n_members=None, cleanup=False
):
    """Write a data table to SQLite.

    Rows are upserted on the FC unique index. Columns not present in the table
    (for example ensemble members from another file) keep their existing
    values. With cleanup=True the existing rows for the written keys are
    deleted first, restoring the old "replace the whole key" behaviour.

    Args:
        table: a data table
        sqlite_file: file name
        param: parameter descriptor
        model_name: model name used in file path and data table
        n_members: number of ensemble members
        cleanup: delete existing rows for the key before upserting
    """
    logger.debug("Writing to %s", sqlite_file)
    fcdate = float(table.iloc[0]["fcst_dttm"])
    leadtime = float(table.iloc[0]["lead_time"])
    mbr_columns = sorted(list_mbr_columns(table, model_name))
    # Serialise table creation and schema patching across parallel workers.
    # Without this lock, two threads can both find no FC table (or both add
    # the same missing member column), producing "duplicate column name" errors.
    file_lock = _get_file_lock(sqlite_file)
    with file_lock:
        sqlite_path = os.path.dirname(sqlite_file)
        if sqlite_path and not os.path.isdir(sqlite_path):
            logger.info("SQLITE: Creating directory %s.", sqlite_path)
            os.makedirs(sqlite_path, exist_ok=True)
        with closing(sqlite3.connect(sqlite_file)) as con:
            _configure_connection(con)

            # Check for the FC table, not for the file: a run that was
            # interrupted (or failed) may have left a file without FC table.
            if not _fc_table_exists(con):
                logger.info("SQLITE: Creating FC table in %s.", sqlite_file)
                # NOTE: to create the table, we should know exactly which members are present
                #       so we don't create empty columns.
                _create_fc_table(con, param, model_name, n_members, mbr_columns)
            # Patch any member columns that may be missing (e.g. table created
            # by an earlier partial flush, a different run or another process).
            ensure_table_columns(con, mbr_columns)

            conflict_columns = _upsert_conflict_columns(con, table)
            if conflict_columns is None:
                raise RuntimeError(
                    f"SQLITE: no unique FC index in {sqlite_file} usable for upsert."
                )

            logger.debug("leadtime (h): %i", leadtime)
            # Cleanup and upsert share one transaction, so a failed write
            # cannot lose rows that were already deleted.
            with con:
                if cleanup:
                    db_cleanup(param, fcdate, leadtime, con, commit=False)
                _upsert_table(con, table, conflict_columns)


def db_cleanup(param, fcd, leadtime, con, commit=True):
    """Delete any prior version of the date (it's a primary key).

    Args:
        param: parameter description
        fcd: forecast date (string or int)
        leadtime: lead time (int)
        con: database connection
        commit: commit the deletion (set False when part of a larger transaction)
    """
    logger.debug("Cleanup: ldt=%i", leadtime)
    cleanup = "DELETE from FC WHERE fcst_dttm=? AND lead_time=?"
    # NOTE: variables must be cast to int (or float),
    # because "numpy.int64" will not work in SQL
    # FIXME: better look whether the table has a "level" column!
    with closing(con.cursor()) as cur:
        cn1 = cur.execute("select name from PRAGMA_TABLE_INFO('FC')").fetchall()
        colnames = [x[0] for x in cn1]

        if param["level_name"] in colnames:
            lname = param["level_name"]
            parlev = int(param["level"])
            cleanup = cleanup + f" AND {_quote_identifier(lname)}=?"
            cur.execute(
                cleanup,
                (
                    float(fcd),
                    float(leadtime),
                    int(parlev),
                ),
            )
        else:
            cur.execute(
                cleanup,
                (
                    float(fcd),
                    float(leadtime),
                ),
            )
    if commit:
        con.commit()


def db_create(sqlite_file, param, model_name, n_members=None, member_columns=None):
    """Create new SQLite file with FC table.

    The table is designed to be compatible with harp, so uses the same
    "unique index" approach.
    The table contains one column per member: <model_name>_mbr001

    NOTE: no longer called internally: write_to_sqlite() creates the table
          itself via _create_fc_table(). Kept for external code that imports
          it (it is re-exported by the grib2sqlite.grib2sqlite facade).

    Args:
      sqlite_file: SQLite file name
      param: parameter descriptor
      model_name: model name used in SQLite
      n_members: number of ensemble members (default=None for deterministic)
      member_columns: list of ensemble members available in the data
    """
    with closing(sqlite3.connect(sqlite_file)) as con:
        _configure_connection(con)
        _create_fc_table(con, param, model_name, n_members, member_columns)


def _fc_table_exists(con):
    """Return True if the database already contains an FC table."""
    row = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'FC'"
    ).fetchone()
    return row is not None


def _create_fc_table(con, param, model_name, n_members=None, member_columns=None):
    """Create the FC table and its unique index on an open connection.

    Both statements run in one explicit transaction, so the table never exists
    without its unique index (the sqlite3 module would otherwise commit each
    DDL statement separately). Both use IF NOT EXISTS, so this is a no-op for
    an existing table.

    Args:
      con: open sqlite3 connection (not inside a transaction)
      param: parameter descriptor
      model_name: model name used in SQLite
      n_members: number of ensemble members (default=None for deterministic)
      member_columns: list of ensemble members available in the data
    """
    primary_keys, other_keys = fctable_definition(param, model_name, n_members, member_columns)
    all_keys = {**primary_keys, **other_keys}

    fc_def = (
        "CREATE table if not exists FC ("
        + ",".join(
            f"{_quote_identifier(name)} {sql_type}"
            for name, sql_type in all_keys.items()
        )
        + ")"
    )
    index_name = "index_" + "_".join(primary_keys.keys())
    pk_def = (
        "CREATE unique INDEX IF NOT EXISTS "
        + _quote_identifier(index_name)
        + " ON FC("
        + ",".join(_quote_identifier(name) for name in primary_keys)
        + ")"
    )

    # IMMEDIATE takes the write lock up front (waiting up to busy_timeout if
    # another process holds it).
    con.execute("BEGIN IMMEDIATE")
    try:
        con.execute(fc_def)
        con.execute(pk_def)
    except BaseException:
        con.rollback()
        raise
    con.commit()


def fctable_definition(param, model_name, n_members=None, member_columns=None):
    """Create the SQL command for FCtable definition.

    Args:
        param: the parameter descriptor
        model_name: model name to be used
        n_members: number of ensemble members (None for deterministic)
        member_columns: the ensemble members present in the data table.

    Returns:
        S(primary_keys, other_keys) – both ordered dicts of {col_name: sql_type}
    """
    primary_keys = {"fcst_dttm": "DOUBLE", "lead_time": "DOUBLE", "SID": "INT"}
    if param["level_name"] in ["z", "p", "ml"]:
        primary_keys[param["level_name"]] = "INT"

    other_keys = {
        "lat": "DOUBLE",
        "lon": "DOUBLE",
        "valid_dttm": "INT",
        "parameter": "TEXT",
        "units": "TEXT",
    }
    # when called from the standard pipeline, 
    # the table should always have data columns for all members
    if member_columns: # EPS members (<>_mbrxxx) or DET (<>_det)
        for col in sorted(member_columns):
            other_keys[col] = "DOUBLE"
    else:
        # This shouldn't happen, but maybe if the table definition is called directly from outside.
        logger.warning("no model data columns. Going by n_members=%s", n_members)
        if n_members is None: # Deterministic
            other_keys[mbr_column(model_name, None)] = "DOUBLE"
        else:
            for mbr in range(0, n_members):
                other_keys[mbr_column(model_name, mbr)] = "DOUBLE"
    # for T2m correction
    if "T2m" in param["harp_param"]:
        other_keys["elev"] = "DOUBLE"
        other_keys["model_elevation"] = "DOUBLE"

    return primary_keys, other_keys
