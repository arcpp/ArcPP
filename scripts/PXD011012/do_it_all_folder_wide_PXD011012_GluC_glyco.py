#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os

offsets = {
    "Rep1_0-7": {
        "StS_Hvo_iTRAQ_0+7_GluC_Big12_3h_09102016.mzML": -2,
    },
    "Rep1_1-2-0": {
        "StS_Hvo_iTRAQ_1+2+0_GluC_Big12_3h_09102016.mzML": 0,
    },
    "Rep1_1-2": {
        "StS_Hvo_iTRAQ_1+2_GluC_Big12_3h_09102016.mzML": -2,
    },
    "Rep1_1-3-0": {
        "StS_Hvo_iTRAQ_1+3+0_GluC_Big12_3h_09102016.mzML": 0,
    },
    "Rep1_1-3-5": {
        "StS_Hvo_iTRAQ_1+3+5_GluC_Big12_3h_09102016.mzML": 2,
    },
    "Rep1_1-3": {
        "StS_Hvo_iTRAQ_1+3_GluC_Big12_3h_09102016.mzML": -2,
    },
    "Rep1_2-4-0": {
        "StS_Hvo_iTRAQ_2+4+0_GluC_Big12_3h_09102016.mzML": 4,
    },
    "Rep1_2-4-0-II": {
        "StS_Hvo_iTRAQ_2+4+0-II_GluC_Big12_3h_09102016.mzML": 2,
    },
    "Rep1_2-4-6": {
        "StS_Hvo_iTRAQ_2+4+6_GluC_Big12_3h_09102016.mzML": -2,
    },
    "Rep1_2-4": {
        "StS_Hvo_iTRAQ_2+4_GluC_Big12_3h_09102016.mzML": -2,
    },
    "Rep1_3-5-0": {
        "StS_Hvo_iTRAQ_3+5+0_GluC_Big12_3h_09102016.mzML": 0,
    },
    "Rep1_4-6": {
        "StS_Hvo_iTRAQ_4+6_GluC_Big12_3h_09102016.mzML": 0,
    },
    "Rep1_5-6": {
        "StS_Hvo_iTRAQ_5+6+0_GluC_Big12_3h_09102016.mzML": 4,
    },
    "Rep3_0-7-8": {
        "StS_ITRAQ-GluC-0-7-8_18082017.mzML": -4,
    },
    "Rep2_0-7-8": {
        "StS_ITRAQ-GluC-0-7-8_30052017_170625002335.mzML": -2,
    },
    "Rep3_1": {
        "StS_ITRAQ-GluC-1_18082017.mzML": 0,
    },
    "Rep3_1-2-0": {
        "StS_ITRAQ-GluC-1-2-0_18082017.mzML": -4,
    },
    "Rep2_1-2-0": {
        "StS_ITRAQ-GluC-1-2-0_30052017_170625101559.mzML": -2,
    },
    "Rep3_1-3-5-0": {
        "StS_ITRAQ-GluC-1-3-5-0_18082017.mzML": -2,
    },
    "Rep2_1-3-5-0": {
        "StS_ITRAQ-GluC-1-3-5-0_30052017_170625165053.mzML": 0,
    },
    "Rep3_1-3-5": {
        "StS_ITRAQ-GluC-1-3-5_18082017.mzML": -4,
    },
    "Rep2_1-3-5": {
        "StS_ITRAQ-GluC-1-3-5_30052017_170625034101.mzML": 2,
    },
    "Rep3_2-4-6-0": {
        "StS_ITRAQ-GluC-2-4-6-0_18082017.mzML": 0,
    },
    "Rep2_2-4-6-0": {
        "StS_ITRAQ-GluC-2-4-6-0_30052017_170625200827.mzML": 0,
    },
    "Rep3_2-4-6": {
        "StS_ITRAQ-GluC-2-4-6_18082017.mzML": -2,
    },
    "Rep2_2-4-6": {
        "StS_ITRAQ-GluC-2-4-6_30052017_170624210606.mzML": 0,
    },
    "Rep3_3-4-0": {
        "StS_ITRAQ-GluC-3-4-0_18082017.mzML": 2,
    },
    "Rep2_3-4-0": {
        "StS_ITRAQ-GluC-3-4-0_30052017_170625133327.mzML": -2,
    },
    "Rep3_4": {
        "StS_ITRAQ-GluC-4_18082017.mzML": 0,
    },
    "Rep3_5-6-0": {
        "StS_ITRAQ-GluC-5-6-0_18082017.mzML": -4,
    },
    "Rep2_5-6-0": {
        "StS_ITRAQ-GluC-5-6-0_30052017_170625065828.mzML": -2,
    },
}


def main(folder=None, enzyme=None, target_decoy_database=None):
    """
    Workflow for the analysis a dataset with one run per sample.
    Usage:
        python <script_name.py> <folder_with_mzML> <enzyme> <path_to_database>
    """
    # define folder with mzML_files as sys.argv[1]
    mzML_files = []
    for mzml in glob.glob(os.path.join(folder, "*.mzML")):
        mzML_files.append(os.path.basename(mzml))
    offset_files = []
    for sample in offsets.keys():
        for spec_file in offsets[sample].keys():
            offset_files.append(spec_file)
    for mzml in mzML_files:
        if mzml not in offset_files:
            print("mzML file in folder but NOT in offset dict: {}".format(mzml))
            exit()
    for sample in offset_files:
        if sample not in mzML_files:
            print("Sample in offset dict but mzML file NOT in folder: {}".format(sample))
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
        "enzyme": enzyme,
        "precursor_mass_tolerance_minus": 10,
        "precursor_mass_tolerance_plus": 10,
        "frag_mass_tolerance": 10,
        "frag_mass_tolerance_unit": "ppm",
        "rounded_mass_decimals": 2,
        "-xmx": "32g",
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
        "max_missed_cleavages": 3,
    }

    # glycans are defined as variable modifications
    # Hex and Hex(1)HexA(1) (=1427) are existing unimod modifications
    Hvo_Glyco = [
        "",
        "N,opt,any,Hex",
        "N,opt,any,1427",
        "N,opt,any,Hex(1)HexA(2),C18H26O17",
        "N,opt,any,Hex(1)HexA(3),C24H34O23",
        "N,opt,any,Hex(1)HexA(2)MeHexA(1)Hex(1),C31H46O28",
        "N,opt,any,Hex(1)HexA(2)MeHexA(1),C25H36O23",
        "N,opt,any,SO3Hex(1),C6H10O8S1",
        "N,opt,any,SO3Hex(1)Hex(1),C12H20O13S1",
        "N,opt,any,SO3Hex(1)Hex(2),C18H30O18S1",
        "N,opt,any,SO3Hex(1)Hex(2)dHex(1),C24H40O22S1",
    ]

    uc = ursgal.UController(profile=mass_spectrometer, params=params)

    combined_pep_result_files = []
    for n, sample in enumerate(sorted(offsets.keys(), reverse=True)):
        validated_result_files = []
        for search_engine in search_engines:
            engine_results_validated = []
            for n, mod in enumerate(Hvo_Glyco):
                results_one_mod = []
                for spec_file in sorted(offsets[sample].keys()):
                    basename = spec_file
                    dirname = os.path.join(folder)
                    offset = offsets[sample][basename]
                    spec_file_path = os.path.join(dirname, basename)
                    if offset == "skip":
                        continue
                    uc.params["machine_offset_in_ppm"] = offset
                    mgf_file = uc.convert(
                        input_file=spec_file_path,
                        engine="mzml2mgf_2_0_0",
                    )

                    if n == 0:
                        uc.params["modifications"] = [
                            "C,fix,any,Carbamidomethyl",
                            "M,opt,any,Oxidation",
                            "*,fix,N-term,iTRAQ4plex",
                            "K,opt,any,iTRAQ4plex",
                            "Y,opt,any,iTRAQ4plex",
                        ]
                    else:
                        uc.params["modifications"] = [
                            "C,fix,any,Carbamidomethyl",
                            "M,opt,any,Oxidation",
                            "*,fix,N-term,iTRAQ4plex",
                            "K,opt,any,iTRAQ4plex",
                            "Y,opt,any,iTRAQ4plex",
                            "S,opt,any,Hex(2)",
                            "T,opt,any,Hex(2)",
                        ]
                        uc.params["modifications"].append(mod)
                        uc.params["prefix"] = mod.split(",")[3]

                    search_result = uc.search_mgf(
                        input_file=mgf_file,
                        engine=search_engine,
                    )
                    uc.params["prefix"] = ""

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

                    results_one_mod.append(unified_search_results)

                uc.params["prefix"] = sample
                merged_1engine_1mod_1sample = uc.execute_misc_engine(
                    input_file=results_one_mod,
                    engine="merge_csvs",
                    # merge_duplicates=True,
                )
                uc.params["prefix"] = ""
                # engine_results_unvalidated.append(merged_1engine_1mod_1sample)

                validated_csv = uc.validate(
                    input_file=merged_1engine_1mod_1sample,
                    engine=validation_engine,
                )
                engine_results_validated.append(validated_csv)

            merged_1engine_all_mods_validated = uc.execute_misc_engine(
                input_file=engine_results_validated,
                engine="merge_csvs",
                merge_duplicates=False,
            )
            validated_result_files.append(merged_1engine_all_mods_validated)

        uc.params["prefix"] = sample
        combined_pep_validated = uc.combine_search_results(
            input_files=validated_result_files,
            engine="combine_pep_1_0_0",
        )
        uc.params["prefix"] = ""
        uc.params["csv_filter_rules"] = [
            # ["Is decoy", "equals", "false"],
            ["combined PEP", "lte", 0.01],
            ["Conflicting uparam", "contains_not", "enzyme"],
        ]
        filtered_validated_results = uc.execute_misc_engine(
            input_file=combined_pep_validated,
            engine="filter_csv",
        )
        combined_pep_result_files.append(filtered_validated_results)

        # uc.params['peptide_forest_initial_engine'] = 'msfragger_2_3'
        # uc.params['peptide_forest_file_params'] = {}
        # uc.params['prefix'] = 'peptide_forest_' + sample
        # peptide_forest_validated =  uc.validate(
        #     input_file=unvalidated_result_files,
        #     engine='peptide_forest',
        # )
        # uc.params['csv_filter_rules'] = [
        #     ['Is decoy', 'equals', 'false'],
        #     ['q-value_RF-reg','lte', 0.01],
        #     ['Conflicting uparam', 'contains_not', 'enzyme'],
        # ]
        # filtered_peptide_forest = uc.execute_misc_engine(
        #     input_file = peptide_forest_validated,
        #     engine='filter_csv',
        # )
        # peptide_forest_result_files.append(filtered_peptide_forest)

    uc.params["prefix"] = ""
    results_all_combined_pep = uc.execute_misc_engine(
        input_file=combined_pep_result_files,
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
        input_file=results_all_combined_pep,
        engine="sanitize_csv",
    )

    uc.params["prefix"] = "Glyco_everywhere"
    uc.params["csv_filter_rules"] = [
        ["Modifications", "contains", "Hex"],
        # ['Sequence','contains_glycosite', 'N[^P][ST]']
    ]
    Glyco_filtered = uc.execute_misc_engine(
        input_file=sanitized_combined_results,
        engine="filter_csv",
    )

    uc.params["prefix"] = "Glyco_glycosite"
    uc.params["csv_filter_rules"] = [
        ["Modifications", "contains", "Hex"],
        ["Sequence", "contains_glycosite", "N[^P][ST]"],
    ]
    Glyco_filtered = uc.execute_misc_engine(
        input_file=sanitized_combined_results,
        engine="filter_csv",
    )
    uc.params["prefix"] = ""


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
        exit()
    main(
        folder=sys.argv[1],
        enzyme=sys.argv[2],
        target_decoy_database=sys.argv[3],
    )
