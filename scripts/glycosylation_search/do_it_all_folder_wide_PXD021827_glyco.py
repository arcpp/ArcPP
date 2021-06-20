#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os


offsets = {
    "WT_Rep1": {
        "OF_20190919_OEF_01.mzML": -0.3,
    },
    "ksgA_Rep1": {
        "OF_20190919_OEF_03.mzML": -0.2,
    },
    "ksgA-reint_Rep1": {
        "OF_20190919_OEF_05.mzML": -0.2,
    },
    "WT_Rep2": {
        "OF_20190920_OEF_03.mzML": -0.1,
    },
    "WT_Rep3": {
        "OF_20190920_OEF_04.mzML": -0.1,
    },
    "ksgA_Rep2": {
        "OF_20190920_OEF_05.mzML": -0.2,
    },
    "ksgA_Rep3": {
        "OF_20190921_OEF_01.mzML": -0.1,
    },
    "ksgA-reint_Rep2": {
        "OF_20190921_OEF_02.mzML": 0.1,
    },
    "ksgA-reint_Rep3": {
        "OF_20190921_OEF_03.mzML": -0.1,
    },
    "ksgA-reintE48A_Rep1": {
        "OF_20190921_OEF_04.mzML": 0.0,
    },
    "ksgA-reintE48A_Rep3": {
        "OF_20190923_OEF_01.mzML": 0.0,
    },
    "ksgA-reintE48A_Rep2": {
        "OF_20190923_OEF_03.mzML": -0.1,
    },
}


def main(folder=None, enzyme=None, target_decoy_database=None):
    """"""
    # define folder with raw_files as sys.argv[1]
    mzml_files = []
    for mzml in glob.glob(os.path.join(folder, "*.mzML")):
        mzml_files.append(os.path.basename(mzml))
    offset_files = []
    for sample in offsets.keys():
        for spec_file in offsets[sample].keys():
            offset_files.append(spec_file)
    for mzml in mzml_files:
        if mzml not in offset_files:
            print("mzml file in folder but NOT in offset dict: {}".format(mzml))
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
        "precursor_mass_tolerance_minus": 8,
        "precursor_mass_tolerance_plus": 8,
        "frag_mass_tolerance": 0.4,
        "frag_mass_tolerance_unit": "da",
        "rounded_mass_decimals": 2,
        "-xmx": "24g",
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
        "calibrate_mass": True,
    }

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
    for n, sample in enumerate(sorted(offsets.keys(), reverse=True)[:]):
        print(n)
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

                    if n == 0:
                        prefix = ""
                    else:
                        prefix = "{0}_".format(
                            mod.split(",")[3],
                        )

                    unified_search_results = os.path.join(
                        folder,
                        search_engine,
                        "{0}{1}_{2}_pmap_unified.csv".format(
                            prefix, spec_file.replace(".mzML", ""), search_engine
                        ),
                    )
                    if os.path.exists(unified_search_results):
                        pass
                    else:
                        uc.params["machine_offset_in_ppm"] = offset
                        # mzml_file = uc.convert(
                        #     input_file=spec_file_path,
                        #     engine="thermo_raw_file_parser_1_1_2",
                        # )
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

                continue
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

            continue
            merged_1engine_all_mods_validated = uc.execute_misc_engine(
                input_file=engine_results_validated,
                engine="merge_csvs",
                merge_duplicates=False,
            )
            validated_result_files.append(merged_1engine_all_mods_validated)

        continue
        uc.params["prefix"] = sample
        combined_pep_validated = uc.combine_search_results(
            input_files=validated_result_files,
            engine="combine_pep_1_0_0",
        )
        uc.params["prefix"] = ""
        uc.params["csv_filter_rules"] = [
            ["Is decoy", "equals", "false"],
            ["combined PEP", "lte", 0.01],
            ["Conflicting uparam", "contains_not", "enzyme"],
        ]
        filtered_validated_results = uc.execute_misc_engine(
            input_file=combined_pep_validated,
            engine="filter_csv",
        )
        combined_pep_result_files.append(filtered_validated_results)

    exit()  # uc.params["prefix"] = ""
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
