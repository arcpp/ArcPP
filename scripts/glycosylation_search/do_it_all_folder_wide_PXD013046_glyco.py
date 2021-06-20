#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os

offsets = {
    "ABI_Mem_exp_rep1": {
        "10_131129073758.mzML": 2,
        "10.mzML": 0,
    },
    "ABI_Mem_exp_rep2": {
        "10A_131129150241.mzML": 0,
        "10A.mzML": 2,
    },
    "ABI_Mem_exp_rep3": {
        "10B_131129222708.mzML": 4,
        "10B.mzML": 0,
    },
    "ABI_Mem_exp_rep4": {
        "10C_131130055123.mzML": 2,
        "10C.mzML": 4,
    },
    "ABI_Cyt_exp_rep1": {
        "9_131125033003.mzML": 0,
        "9.mzML": -4,
    },
    "ABI_Cyt_exp_rep2": {
        "9A_131125110110.mzML": 0,
        "9A.mzML": -2,
    },
    "ABI_Cyt_exp_rep3": {
        "9B_131125223541.mzML": -2,
        "9B.mzML": -2,
    },
    "ABI_Cyt_exp_rep4": {
        "9C_131126060708.mzML": 0,
        "9C.mzML": -2,
    },
    "ABI_Mem_stat_rep1": {
        "12_131130131543.mzML": 0,
        "12.mzML": -2,
    },
    "ABI_Mem_stat_rep2": {
        "12A_131130204006.mzML": 0,
        "12A.mzML": 0,
    },
    "ABI_Mem_stat_rep3": {
        "12B_131201040426.mzML": 0,
        "12B.mzML": 2,
    },
    "ABI_Mem_stat_rep4": {
        "12C_131201112849.mzML": 2,
        "12C.mzML": 0,
    },
    "ABI_Cyt_stat_rep1": {
        "11_131126133808.mzML": -2,
        "11.mzML": -4,
    },
    "ABI_Cyt_stat_rep2": {
        "11A_131126210917.mzML": -2,
        "11A.mzML": 0,
    },
    "ABI_Cyt_stat_rep3": {
        "11B_131127044005.mzML": -0,
        "11B.mzML": -2,
    },
    "ABI_Cyt_stat_rep4": {
        "11C_131127121052.mzML": -2,
        "11C.mzML": -4,
    },
    "HVLON3_Mem_exp_rep1": {
        "6_top20.mzML": 2,
        "6_botttom20.mzML": 2,
    },
    "HVLON3_Mem_exp_rep2": {
        "6A_botttom20.mzML": 2,
        "6A_top20.mzML": 0,
    },
    "HVLON3_Mem_exp_rep3": {
        "6B_131101122740.mzML": 0,
        "6B.mzML": 0,
    },
    "HVLON3_Mem_exp_rep4": {
        "6C_131128164906.mzML": 2,
        "6C.mzML": 0,
    },
    "HVLON3_Cyt_exp_rep1": {
        "5_131112222612.mzML": 0,
        "5.mzML": 2,
    },
    "HVLON3_Cyt_exp_rep2": {
        "5A_131113055702.mzML": 4,
        "5A.mzML": 2,
    },
    "HVLON3_Cyt_exp_rep3": {
        "5B_131113132754.mzML": 0,
        "5B.mzML": 0,
    },
    "HVLON3_Cyt_exp_rep4": {
        "5C_131113205840.mzML": 0,
        "5C.mzML": 0,
    },
    "HVLON3_Mem_stat_rep1": {
        "8.mzML": 2,
        "8_131101201658.mzML": 2,
    },
    "HVLON3_Mem_stat_rep2": {
        "8A_131102040617.mzML": 0,
        "8A.mzML": 2,
    },
    "HVLON3_Mem_stat_rep3": {
        "8B_131102115537.mzML": 0,
        "8B.mzML": 2,
    },
    "HVLON3_Mem_stat_rep4": {
        "8C_131129001334.mzML": 0,
        "8C.mzML": -2,
    },
    "HVLON3_Cyt_stat_rep1": {
        "7.mzML": 0,
        "7_131122192428.mzML": -2,
    },
    "HVLON3_Cyt_stat_rep2": {
        "7A_131123025552.mzML": -2,
        "7A.mzML": 0,
    },
    "HVLON3_Cyt_stat_rep32": {
        "7B_131125150408.mzML": -2,
        "7B.mzML": 0,
    },
    "HVLON3_Cyt_stat_rep4": {
        "7C_131124195907.mzML": -2,
        "7C.mzML": -2,
    },
    "WT_Mem_exp_rep1": {
        "2_top20.mzML": 2,
        "2_botttom20.mzML": 0,
    },
    "WT_Mem_exp_rep2": {
        "2A_botttom20.mzML": 0,
        "2A_top20_140615200858.mzML": 0,
    },
    "WT_Mem_exp_rep3": {
        "2B_131030002035.mzML": 2,
        "2B.mzML": 2,
    },
    "WT_Mem_exp_rep4": {
        "2C_131128020026.mzML": 2,
        "2C.mzML": 0,
    },
    "WT_Cyt_exp_rep1": {
        "1.mzML": 0,
        "1_131111162318.mzML": 4,
    },
    "WT_Cyt_exp_rep2": {
        "1A_131111235403.mzML": 2,
        "1A.mzML": 2,
    },
    "WT_Cyt_exp_rep3": {
        "1B_131112072445.mzML": -2,
        "1B.mzML": -2,
    },
    "WT_Cyt_exp_rep4": {
        "1C_131112145528.mzML": 0,
        "1C.mzML": 4,
    },
    "WT_Mem_stat_rep1": {
        "4_top20.mzML": 4,
        "4_botttom20.mzML": 0,
    },
    "WT_Mem_stat_rep2": {
        "4A_botttom20.mzML": 0,
        "4A_top20.mzML": 2,
    },
    "WT_Mem_stat_rep3": {
        "4B_131030234822.mzML": -8,
        "4B.mzML": 0,
    },
    "WT_Mem_stat_rep4": {
        "4C_131128092443.mzML": 2,
        "4C.mzML": 4,
    },
    "WT_Cyt_stat_rep1": {
        "3.mzML": -4,
        "3_131114042927.mzML": -2,
    },
    "WT_Cyt_stat_rep2": {
        "3A_131114120015.mzML": -4,
        "3A.mzML": 0,
    },
    "WT_Cyt_stat_rep3": {
        "3B_131114193105.mzML": -4,
        "3B.mzML": -4,
    },
    "WT_Cyt_stat_rep4": {
        "3C_131115030157.mzML": 0,
        "3C.mzML": -2,
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
        "precursor_mass_tolerance_minus": 8,
        "precursor_mass_tolerance_plus": 8,
        "frag_mass_tolerance": 0.4,
        "frag_mass_tolerance_unit": "da",
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
        "max_missed_cleavages": 2,
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
                            # "C,fix,any,Carbamidomethyl",
                            "M,opt,any,Oxidation",
                            "*,opt,Prot-N-term,Acetyl",
                        ]
                    else:
                        uc.params["modifications"] = [
                            # "C,fix,any,Carbamidomethyl",
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
