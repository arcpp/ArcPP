#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os
import csv


def merge_files(file_in=""):
    """
    Filters ArcPP_results_PSMs.csv for glycopeptides and uses SugarPy's
    glycopeptide fragmentor to search for glycopeptide-specific fragment ions.

    Usage:
        python <script_name.py> <path_to_ArcPP_PSMs_file>
    """
    uc = ursgal.UController()
    uc.params["prefix"] = "Glyco_everywhere"
    uc.params["csv_filter_rules"] = [
        ["Modifications", "contains", "Hex"],
        # ['Sequence','contains_glycosite', 'N[^P][ST]']
    ]
    Glyco_filtered = uc.execute_misc_engine(
        input_file=file_in,
        engine="filter_csv",
    )

    uc.params["prefix"] = ""

    PRIDE_folders = [
        "PXD021874",
        "PXD007061",  # should not include requirement for 1 oxonoium ion
        "PXD013046",  # should not include requirement for 1 oxonoium ion
        "PXD011012",
        "PXD006877",  # should not include requirement for 1 oxonoium ion
        "PXD011218",  # should not include requirement for 1 oxonoium ion
        "PXD009116",  # should not include requirement for 1 oxonoium ion
        "PXD011056",
        "PXD000202",  # should not include requirement for 1 oxonoium ion
        "PXD011015",
        "PXD011050",
        "PXD014974",  # should not include requirement for 1 oxonoium ion
        "PXD010824",
    ]

    all_gfrag_ion_unfiltered = []
    for pride_id in PRIDE_folders[:]:
        mzML_files = []
        for mzml in glob.glob(os.path.join("{0}".format(pride_id), "*.idx.gz")):
            mzML_files.append(mzml)

        uc.params["csv_filter_rules"] = [
            ["Dataset", "contains", pride_id],
        ]
        uc.params["prefix"] = pride_id
        filtered_glyco = uc.execute_misc_engine(
            input_file=Glyco_filtered,
            engine="filter_csv",
        )
        uc.params["prefix"] = ""

        uc.params.update(
            {
                "mzml_input_files": mzML_files,
                "frag_mass_tolerance": 50,
                "frag_mass_tolerance_unit": "ppm",
                "min_oxonium_ions": 1,
                "min_y_ions": 2,
            }
        )
        glyco_frag_ions_results = uc.execute_misc_engine(
            input_file=filtered_glyco,
            engine="glycopeptide_fragmentor",
            force=False,
        )

        all_gfrag_ion_unfiltered.append(glyco_frag_ions_results)

    merged_results = uc.execute_misc_engine(
        input_file=all_gfrag_ion_unfiltered,
        engine="merge_csvs",
    )


if __name__ == "__main__":
    print(sys.argv[1:])
    merge_files(file_in=sys.argv[1:])
