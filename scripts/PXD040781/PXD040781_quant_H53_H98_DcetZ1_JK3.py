#!/usr/bin/env python
import pathlib
import sys
import ursgal
import os
import glob


def main(mzml_folder, merged_result):
    params = {
        "isotopic_distribution_tolerance": 5,
        "normalize_intensities": True,
        "integrate_peak_areas": False,
        "only_precursor_charge": False,
        "match_between_runs": True,
        "match_between_runs_RT_window": 0.5,  # maybe 1?
        "require_msms_id": False,
        "bayesian_fold_change": True,
        "bayesian_fold_change_control_condition": "JK3_L",  # H = Healthy, A = ALS
        "fold_change_cutoff": 0.1,
        "markov_chain_iterations": 3000,
        "markov_chain_burn_in_iterations": 1000,
        "use_shared_peptides": True,
        "random_seed": 200,
    }
    uc = ursgal.UController(
        verbose=True,
        params=params,
        profile="QExactive+",
    )

    experiment_setup = {
        "1": {
            "FileName": "210224_StS_Hfx_H53_E1_60min",
            "Condition": "H53_E",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "2": {
            "FileName": "210224_StS_Hfx_H53_E2_60min",
            "Condition": "H53_E",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "3": {
            "FileName": "210224_StS_Hfx_H53_L1_60min",
            "Condition": "H53_L",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "4": {
            "FileName": "210224_StS_Hfx_H53_L2_60min",
            "Condition": "H53_L",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "5": {
            "FileName": "210224_StS_Hfx_cetZ_E1_60min",
            "Condition": "cetZ_E",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "6": {
            "FileName": "210224_StS_Hfx_cetZ_E2_60min",
            "Condition": "cetZ_E",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "7": {
            "FileName": "210224_StS_Hfx_cetZ_L1_60min",
            "Condition": "cetZ_L",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "8": {
            "FileName": "210224_StS_Hfx_cetZ_L2_60min",
            "Condition": "cetZ_L",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "9": {
            "FileName": "210224_StS_Hfx_H98_E1_60min",
            "Condition": "H98_E",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "10": {
            "FileName": "210224_StS_Hfx_H98_E2_60min",
            "Condition": "H98_E",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "11": {
            "FileName": "210224_StS_Hfx_H98_L1_60min",
            "Condition": "H98_L",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "12": {
            "FileName": "210224_StS_Hfx_H98_L2_60min",
            "Condition": "H98_L",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "13": {
            "FileName": "210224_StS_Hfx_JK3_E1_60min",
            "Condition": "JK3_E",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "14": {
            "FileName": "210224_StS_Hfx_JK3_E2_60min",
            "Condition": "JK3_E",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "15": {
            "FileName": "210224_StS_Hfx_JK3_L1_60min",
            "Condition": "JK3_L",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "16": {
            "FileName": "210224_StS_Hfx_JK3_L2_60min",
            "Condition": "JK3_L",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
    }

    uc.params["experiment_setup"] = experiment_setup
    uc.params["quantification_evidences"] = merged_result

    mzml_files = []
    for mzml in glob.glob(os.path.join(mzml_folder, "*.mzML")):
        mzml_files.append(mzml)

    quantified_peaks = uc.quantify(
        input_file=mzml_files, multi=True, engine="flash_lfq_1_1_1"
    )


if __name__ == "__main__":
    main(
        sys.argv[1],
        sys.argv[2],
    )
    pass
