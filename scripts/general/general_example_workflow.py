#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os

samples = {
    "Sample1": [
        "File1.raw",
        "File2.raw",
        "File3.raw",
    ],
    "Sample2": [
        "File4.raw",
        "File5.raw",
        "File6.raw",
    ],
}


def main(folder=None, target_decoy_database=None):
    """
    Workflow for the analysis of a dataset with muliple MS raw files per sample.
    Usage:
        python do_it_all_folder_wide_general_workflow.py <folder_with_raw_files> <path_to_database>
    """
    # define folder with raw_files as sys.argv[1]
    raw_files = []
    for raw in glob.glob(os.path.join(folder, "*.raw")):
        raw_files.append(os.path.basename(raw))
    sample_files = []
    for sample in samples.keys():
        for spec_file in samples[sample]:
            sample_files.append(spec_file)
    for raw in raw_files:
        if raw not in sample_files:
            print("raw file in folder but NOT in samples dict: {}".format(raw))
            exit()
    for sample in sample_files:
        if sample not in raw_files:
            print("Sample in samples dict but raw file NOT in folder: {}".format(sample))
            exit()

    mass_spectrometer = "QExactive+"
    search_engines = [
        "xtandem_vengeance",
        "msfragger_2_3",
        "msgfplus_v2019_07_03",
    ]

    validation_engine = "percolator_3_4_0"

    params = {
        "database": target_decoy_database,
        "enzyme": "trypsin",
        "precursor_mass_tolerance_minus": 10,
        "precursor_mass_tolerance_plus": 10,
        "precursor_mass_tolerance_unit": "ppm",
        "frag_mass_tolerance": 10,
        "frag_mass_tolerance_unit": "ppm",
        "max_missed_cleavages": 2,
        "-xmx": "32g",
        "modifications": [
            "C,fix,any,Carbamidomethyl",
            "M,opt,any,Oxidation",
            "*,opt,Prot-N-term,Acetyl",
        ],
        "rounded_mass_decimals": 2,
        "peptide_mapper_class_version": "UPeptideMapper_v4",
        "use_pyqms_for_mz_calculation": True,
        "percolator_post_processing": "mix-max",
        "psm_defining_colnames": [
            "Spectrum Title",
            "Sequence",
            "Modifications",
            "Charge",
            "Is decoy",
        ],
    }

    uc = ursgal.UController(profile=mass_spectrometer, params=params)

    all_result_files = []
    for n, sample in enumerate(samples.keys()):
        validated_result_files = []
        combined_pep_result_files = []
        for search_engine in search_engines:
            results = []
            for spec_file in samples[sample]:
                basename = spec_file
                dirname = os.path.join(folder)
                spec_file_path = os.path.join(dirname, basename)

                mzml_file = uc.convert(
                    input_file=spec_file_path,
                    engine="thermo_raw_file_parser_1_1_2",
                )

                mgf_file = uc.convert(
                    input_file=mzml_file,
                    engine="mzml2mgf_2_0_0",
                )

                search_result = uc.search_mgf(
                    input_file=mgf_file,
                    engine=search_engine,
                )

                converted_result = uc.convert(
                    input_file=search_result,
                    guess_engine=True,
                )

                mapped_results = uc.execute_misc_engine(
                    input_file=converted_result,
                    engine="upeptide_mapper",
                )

                unified_search_results = uc.execute_misc_engine(
                    input_file=mapped_results, engine="unify_csv"
                )

                results.append(unified_search_results)

            uc.params["prefix"] = sample
            results_one_engine = uc.execute_misc_engine(
                input_file=results,
                engine="merge_csvs",
                # merge_duplicates=True,
            )
            uc.params["prefix"] = ""

            validated_csv = uc.validate(
                input_file=results_one_engine,
                engine=validation_engine,
            )

            validated_result_files.append(validated_csv)

        combined_results = uc.combine_search_results(
            input_files=validated_result_files,
            engine="combine_pep_1_0_0",
        )

        uc.params["csv_filter_rules"] = [
            ["combined PEP", "lte", 0.01],
            ["Conflicting uparam", "contains_not", "enzyme"],
        ]
        filtered_combined_results = uc.execute_misc_engine(
            input_file=combined_results,
            engine="filter_csv",
        )
        all_result_files.append(filtered_combined_results)

    results_all_files = uc.execute_misc_engine(
        input_file=all_result_files,
        engine="merge_csvs",
        merge_duplicates=True,
    )

    uc.params.update(
        {
            "validation_score_field": "combined PEP",
            "bigger_scores_better": False,
            "num_compared_psms": 10,
            "accept_conflicting_psms": False,
            "threshold_is_log10": True,
            "score_diff_threshold": 1,
            "psm_defining_colnames": [
                "Spectrum Title",
                "Sequence",
            ],
        }
    )

    sanitized_combined_results = uc.execute_misc_engine(
        input_file=results_all_files,
        engine="sanitize_csv",
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
        exit()
    main(
        folder=sys.argv[1],
        target_decoy_database=sys.argv[2],
    )
