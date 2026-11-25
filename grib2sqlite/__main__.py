"""Command line interface: python3 -m grib2sqlite <options> <grib file(s)>.

Exit status:
    0  all files processed
    1  at least one file could not be processed (or an unexpected error)
    2  usage error: invalid arguments, input files not found, unreadable
       parameter or station list
"""

import os
import sqlite3
import sys
from argparse import ArgumentParser
from glob import glob

import eccodes

from .config import Grib2SqliteError, basedir, logger
from .params import read_param_list
from .pipeline import parse_grib_file
from .spatial import read_station_list

EXIT_OK = 0
EXIT_FAILED = 1
EXIT_USAGE = 2

# Errors that concern one input file: reported, then the next file is processed.
# Anything else is considered a bug and ends the run with a traceback.
FILE_ERRORS = (Grib2SqliteError, OSError, sqlite3.Error, eccodes.CodesInternalError)


def build_parser():
    """Create the command line argument parser."""
    parser = ArgumentParser(prog="extractsqlite")
    # positional argument:
    # we allow multiple GRIB files, but some might even be globs
    parser.add_argument("gribfile", metavar="/path/to/file(s)", nargs="+")
    parser.add_argument(
        "-p",
        metavar="param_file",
        help="parameter list file (json)",
        default=basedir + "/data/param_list_default.json",
    )
    parser.add_argument(
        "-s",
        metavar="station_file",
        help="station list file (csv)",
        default=basedir + "/data/station_list_default.csv",
    )
    parser.add_argument(
        "-m",
        metavar="model_name",
        help="base model name used in SQLite (deterministic fields get a _det suffix)",
        default="TEST",
    )
    parser.add_argument(
        "-t",
        metavar="output_template",
        help="template used for SQLite files",
        default="{MODEL}/{YYYY}/{MM}/FCTABLE_{PP}_{YYYY}{MM}.sqlite",
    )
    parser.add_argument(
        "-o",
        metavar="output_path",
        help="path used for SQLite files",
        default=os.getcwd() + "/FCTABLE",
    )
    parser.add_argument(
        "-d",
        metavar="debug_level",
        type=int,
        help="verbosity level (0...3) default: 1",
        default=1,
    )
    parser.add_argument(
        "-e",
        "--ensemble",
        action="store_true",
        help="treat input as ensemble GRIB file (members stored in _mbr001, _mbr002, ... columns)",
    )
    parser.add_argument(
        "--det-member",
        metavar="N",
        type=int,
        default=0,
        help="member number for records without perturbationNumber when -e is used "
        "(default: 0 -> <model>_mbr000)",
    )
    parser.add_argument(
        "-c",
        "--clim",
        metavar="clim_file",
        default=None,
        help="GRIB file with the model orography (surface z, orog, h or gh) on the "
        "same grid as the data; written as model_elevation for T2m",
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="delete existing rows for the written keys before upserting "
        "(default: merge/upsert only)",
    )
    return parser


def expand_input_files(patterns):
    """Expand file names / glob patterns to a list of existing files.

    Args:
        patterns: list of file names or glob patterns

    Returns:
        (files, missing): the files found (in order, without duplicates) and
        the patterns that did not match any file.
    """
    files = []
    missing = []
    for pattern in patterns:
        matches = sorted(glob(pattern))
        logger.debug(f"{pattern} -> {matches}")
        found = [f for f in matches if os.path.isfile(f)]
        for f in matches:
            if not os.path.isfile(f):
                logger.warning(f"Not a file, ignored: {f}")
        if not found:
            missing.append(pattern)
        for f in found:
            if f not in files:
                files.append(f)
    return files, missing


def main(argv=None):
    """Run the command line interface.

    Args:
        argv: list of arguments (default: sys.argv[1:])

    Returns:
        exit status (see module docstring)
    """
    parser = build_parser()
    args = parser.parse_args(argv)  # exits with status 2 on invalid arguments
    if args.det_member < 0:
        parser.error("--det-member must be a non-negative integer")

    if args.d <= 0:
        logger.setLevel("ERROR")
    elif args.d == 1:
        logger.setLevel("WARNING")
    elif args.d == 2:
        logger.setLevel("INFO")
    else:
        logger.setLevel("DEBUG")

    sqlite_template = args.o + "/" + args.t

    logger.debug(f"Binary path    : {basedir} ")
    logger.info(f"GRIB file      : {args.gribfile}")
    logger.info(f"Parameter list : {args.p}")
    logger.info(f"Station list   : {args.s}")
    logger.info(f"Model name     : {args.m}")
    logger.info(f"Clim file      : {args.clim}")
    logger.info(f"Output files   : {sqlite_template}")

    # Check all input before processing anything.
    file_list, missing = expand_input_files(args.gribfile)
    if missing:
        logger.error("No input files found for: %s", ", ".join(missing))
        return EXIT_USAGE
    # Read the parameter and station lists once (not again for every file).
    try:
        param_list = read_param_list(args.p)
        station_list = read_station_list(args.s)
    except (OSError, ValueError) as err:
        logger.error("Can not read parameter or station list: %s", err)
        return EXIT_USAGE

    failed = []
    for fi, file_name in enumerate(file_list, start=1):
        logger.info(f"{fi:3d}:  Reading {file_name}")
        try:
            gt, gi = parse_grib_file(
                infile=file_name,
                param_list=param_list,
                station_list=station_list,
                sqlite_template=sqlite_template,
                model_name=args.m,
                weights=None,
                ensemble=args.ensemble,
                det_member=args.det_member,
                cleanup=args.cleanup,
                clim_file=args.clim,
            )
        except FILE_ERRORS as err:
            logger.error("%s: %s", file_name, err)
            failed.append(file_name)
            continue
        if gi == 0:
            logger.warning("%s: no matching GRIB records (of %i).", file_name, gt)

    n_files = len(file_list)
    if failed:
        logger.error(
            "Processed %i files: %i OK, %i failed (%s).",
            n_files,
            n_files - len(failed),
            len(failed),
            ", ".join(failed),
        )
        return EXIT_FAILED
    logger.info("Processed %i files: all OK.", n_files)
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
