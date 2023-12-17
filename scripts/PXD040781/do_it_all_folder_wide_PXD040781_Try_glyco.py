#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os
import pprint


offsets = {
    "H53_RM_EL_rep1": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_1a_330ng_60min_Wat25cmBEH_1.mzML": 0.65,
    },
    "H53_RM_EL_rep2": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_1b_330ng_60min_Wat25cmBEH_1.mzML": 1.93,
    },
    "H53_RM_EL_rep3": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_1c_330ng_60min_Wat25cmBEH_1.mzML": 1.61,
    },
    "H53_RM_LL_rep1": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_2a_330ng_60min_Wat25cmBEH_1.mzML": 1.28,
    },
    "H53_RM_LL_rep2": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_2b_330ng_60min_Wat25cmBEH_1.mzML": 1.85,
    },
    "H53_RM_LL_rep3": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_2c_330ng_60min_Wat25cmBEH_1.mzML": 1.98,
    },
    "H53_CM_EL_rep1": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_3a_330ng_60min_Wat25cmBEH_1.mzML": 1.67,
    },
    "H53_CM_EL_rep2": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_3b_330ng_60min_Wat25cmBEH_1.mzML": 1.51,
    },
    "H53_CM_EL_rep3": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_3c_330ng_60min_Wat25cmBEH_1.mzML": 1.04,
    },
    "2174_RM_EL_rep1": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_4a_330ng_60min_Wat25cmBEH_1.mzML": 1.38,
    },
    "2174_RM_EL_rep2": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_4b_330ng_60min_Wat25cmBEH_1.mzML": 1.05,
    },
    "2174_RM_EL_rep3": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_4c_330ng_60min_Wat25cmBEH_1.mzML": 1.17,
    },
    "2174_RM_LL_rep1": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_5a_330ng_60min_Wat25cmBEH_1.mzML": 1.21,
    },
    "2174_RM_LL_rep2": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_5b_330ng_60min_Wat25cmBEH_1.mzML": 1.18,
    },
    "2174_RM_LL_rep3": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_5c_330ng_60min_Wat25cmBEH_1.mzML": 1.96,
    },
    "2176_RM_EL_rep1": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_6a_330ng_60min_Wat25cmBEH_1.mzML": 1.27,
    },
    "2176_RM_EL_rep2": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_6b_330ng_60min_Wat25cmBEH_1.mzML": 1.61,
    },
    "2176_RM_EL_rep3": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_6c_330ng_60min_Wat25cmBEH_1.mzML": 1.34,
    },
    "2176_RM_LL_rep1": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_7a_330ng_60min_Wat25cmBEH_1.mzML": 1.82,
    },
    "2176_RM_LL_rep2": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_7b_330ng_60min_Wat25cmBEH_1.mzML": 1.04,
    },
    "2176_RM_LL_rep3": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_7c_330ng_60min_Wat25cmBEH_1.mzML": 1.86,
    },
    "2176_CM_EL_rep1": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_8a_330ng_60min_Wat25cmBEH_1.mzML": 0.88,
    },
    "2176_CM_EL_rep2": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_8b_330ng_60min_Wat25cmBEH_1.mzML": 0.95,
    },
    "2176_CM_EL_rep3": {
        "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_8c_330ng_60min_Wat25cmBEH_1.mzML": 1.67
    },
}


def main(folder=None, enzyme=None, target_decoy_database=None):
    """
    Workflow for the glycoproteomic analysis of a dataset.
    Usage:
        python <script_name.py> <folder_with_mzML> <enzyme> <path_to_database>
    """
    # define folder with mzML_files as sys.argv[1]
    raw_files = []
    for raw in glob.glob(os.path.join(folder, "*.mzML")):
        raw_files.append(os.path.basename(raw))
    offset_files = []
    for sample in offsets.keys():
        for spec_file in offsets[sample].keys():
            offset_files.append(spec_file)
    for raw in raw_files:
        if raw not in offset_files:
            print("raw file in folder but NOT in offset dict: {}".format(raw))
            pprint.pprint(offset_files)
            exit()

    mass_spectrometer = "QExactive+"
    search_engines = [
        "xtandem_vengeance",
        "msfragger_3_0",
        "msgfplus_v2019_07_03",
    ]

    validation_engine = "percolator_3_4_0"

    params = {
        "database": target_decoy_database,
        "enzyme": enzyme,
        "precursor_mass_tolerance_minus": 10,
        "precursor_mass_tolerance_plus": 10,
        "frag_mass_tolerance": 15,
        "frag_mass_tolerance_unit": "ppm",
        "rounded_mass_decimals": 2,
        "-xmx": "48g",
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
        "max_missed_cleavages": 2,
        "calibrate_mass": "no_calibration",
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
    for n, sample in enumerate(sorted(offsets.keys(), reverse=False)[:]):
        validated_result_files = []
        for search_engine in search_engines:
            engine_results_validated = []
            for n, mod in enumerate(Hvo_Glyco):
                results_one_mod = []
                for spec_file in sorted(offsets[sample].keys()):
                    uc.params["prefix"] = ""
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
                            "*,opt,Prot-N-term,Acetyl",
                        ]
                    else:
                        uc.params["modifications"] = [
                            "C,fix,any,Carbamidomethyl",
                            "M,opt,any,Oxidation",
                            "*,opt,Prot-N-term,Acetyl",
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


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
        exit()
    main(
        folder=sys.argv[1],
        enzyme=sys.argv[2],
        target_decoy_database=sys.argv[3],
    )
