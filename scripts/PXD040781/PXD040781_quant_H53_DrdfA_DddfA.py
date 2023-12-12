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
        "bayesian_fold_change_control_condition": "H53_EL",  # H = Healthy, A = ALS
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
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_1a_330ng_60min_Wat25cmBEH_1",
            "Condition": "H53_EL",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "2": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_1b_330ng_60min_Wat25cmBEH_1",
            "Condition": "H53_EL",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "3": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_1c_330ng_60min_Wat25cmBEH_1",
            "Condition": "H53_EL",
            "Biorep": 3,
            "Fraction": 1,
            "Techrep": 1,
        },
        "4": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_2a_330ng_60min_Wat25cmBEH_1",
            "Condition": "H53_LL",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "5": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_2b_330ng_60min_Wat25cmBEH_1",
            "Condition": "H53_LL",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "6": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_H53_2c_330ng_60min_Wat25cmBEH_1",
            "Condition": "H53_LL",
            "Biorep": 3,
            "Fraction": 1,
            "Techrep": 1,
        },
        "7": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_4a_330ng_60min_Wat25cmBEH_1",
            "Condition": "2174_EL",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "8": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_4b_330ng_60min_Wat25cmBEH_1",
            "Condition": "2174_EL",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "9": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_4c_330ng_60min_Wat25cmBEH_1",
            "Condition": "2174_EL",
            "Biorep": 3,
            "Fraction": 1,
            "Techrep": 1,
        },
        "10": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_5a_330ng_60min_Wat25cmBEH_1",
            "Condition": "2174_LL",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "11": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_5b_330ng_60min_Wat25cmBEH_1",
            "Condition": "2174_LL",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "12": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2174_5c_330ng_60min_Wat25cmBEH_1",
            "Condition": "2174_LL",
            "Biorep": 3,
            "Fraction": 1,
            "Techrep": 1,
        },
        "13": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_6a_330ng_60min_Wat25cmBEH_1",
            "Condition": "2176_EL",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "14": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_6b_330ng_60min_Wat25cmBEH_1",
            "Condition": "2176_EL",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "15": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_6c_330ng_60min_Wat25cmBEH_1",
            "Condition": "2176_EL",
            "Biorep": 3,
            "Fraction": 1,
            "Techrep": 1,
        },
        "16": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_7a_330ng_60min_Wat25cmBEH_1",
            "Condition": "2176_LL",
            "Biorep": 1,
            "Fraction": 1,
            "Techrep": 1,
        },
        "17": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_7b_330ng_60min_Wat25cmBEH_1",
            "Condition": "2176_LL",
            "Biorep": 2,
            "Fraction": 1,
            "Techrep": 1,
        },
        "18": {
            "FileName": "20220622_jcl_SSchulze_Pohlschroder_UPenn_Haloferax_d2176_7c_330ng_60min_Wat25cmBEH_1",
            "Condition": "2176_LL",
            "Biorep": 3,
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
