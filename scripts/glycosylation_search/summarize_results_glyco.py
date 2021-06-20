#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import os
import csv
import pickle
import copy
import pprint
import glob
import itertools
import collections
import math
import numpy as np

# lookup for dataset, lab and location of the file
datasets = {
    "PXD006877": {
        "lab": "Maupin-Furlow",
        "num_spectra": 2498860,
        "num_files": 56,
        "instrument": "Q Exacive Plus",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "Lana-SILAC-_-____pmap_unified_merged_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD007061": {
        "lab": "De_Castro",
        "num_spectra": 6847689,
        "num_files": 118,
        "instrument": "LTQ-Orbitrap Elite",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "labeling",
        "result_file": "Rep___-top20___pmap_unified_merged_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD009116": {
        "lab": "De_Castro",
        "num_spectra": 1322477,
        "num_files": 24,
        "instrument": "LTQ-Orbitrap Elite",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "WT__ed_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD010824": {
        "lab": "Pohlschroder",
        "num_spectra": 131989,
        "num_files": 8,
        "instrument": "Q Exacive Plus",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "Hvo_0405_____pmap_unified_merged_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD011012": {
        "lab": "Pohlschroder",
        "num_spectra": 3446807,
        "num_files": 62,
        "instrument": "Q Exacive Plus",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "Merged_gluc_try_glyco_sanitized.csv",
    },
    "PXD011015": {
        "lab": "Pohlschroder",
        "num_spectra": 473288,
        "num_files": 28,
        "instrument": "Q Exacive Plus",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "ed_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD011050": {
        "lab": "Pohlschroder",
        "num_spectra": 452855,
        "num_files": 16,
        "instrument": "Q Exacive Plus",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "Merged_gluc_try_glyco_sanitized.csv",
    },
    "PXD011056": {
        "lab": "Marchfelder",
        "num_spectra": 1063443,
        "num_files": 32,
        "instrument": "TripleTOF 5600+",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "ZJ_____pmap_unified_merged_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD011218": {
        "lab": "De_Castro",
        "num_spectra": 1926117,
        "num_files": 48,
        "instrument": "LTQ-Orbitrap Elite",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "ed_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD013046": {
        "lab": "De_Castro",
        "num_spectra": 4291777,
        "num_files": 96,
        "instrument": "LTQ-Orbitrap Elite",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "pmap_unified_merged_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD014974": {
        "lab": "Ferreira-Cerca",
        "num_spectra": 301610,
        "num_files": 4,
        "instrument": "Q Exacive HF",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "POL11409X002__A_rep____pmap_unified_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD014974": {
        "lab": "Ferreira-Cerca",
        "num_spectra": 1792418,
        "num_files": 12,
        "instrument": "Orbitrap Fusion",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "pmap_unified_percolator_3_4_0_validated_merged_combine_pep_1_0_0_accepted_merged_sanitized.csv",
    },
    "PXD021874": {
        "lab": "Pohlschroder",
        "num_spectra": 1660539,
        "num_files": 62,
        "instrument": "Q Exacive HF",
        "identified_proteins": dict(),
        "identified_peptides": dict(),
        "folders": "",
        "result_file": "Merged_gluc_try_glyco_sanitized.csv",
    },
}

# parameters
SEQ_Q_VALUE_THRESHOLD = 0.01
PROT_Q_VALUE_THRESHOLD = 0.005

# lookup mapping MS file name to sample
ms_filename2sample = {}


def main(dirpath, skip_old=False, num_specs=1):
    uc = ursgal.UController()
    uc.params.update(
        {
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
    pkl_name = os.path.join(dirpath, "datasets_result.pkl")
    fdr_pkl_name = os.path.join(dirpath, "fdr_result.pkl")
    old_exists = False
    if os.path.exists(pkl_name) and skip_old is True:
        # load results from previous analysis
        # will only add datasets that are not part of it already
        print(">>>>>>>> loading pkl <<<<<<<<<<<")
        results_dict = pickle.load(open(pkl_name, "rb"))
        fdr_dict = pickle.load(open(fdr_pkl_name, "rb"))
        old_exists = True
    else:
        # collect proteins and peptides from result csv,
        # store in dict with all important data
        results_dict = {
            "all": {
                "num_spectra": 0,
                "instrument": set(),
                "lab": set(),
                # protein_groups, proteins and peptides are dicts that contain sets for each level of confidence
                "protein_groups": {
                    "all": set(),
                    "safe_psm": set(),
                    "safe_seq": set(),
                    "safe_seq_num_spec": set(),
                    "safe_seq_num_spec_0005": set(),
                },
                "proteins": {
                    "all": set(),
                    "safe_psm": set(),
                    "safe_seq": set(),
                    "safe_seq_num_spec": set(),
                    "safe_seq_num_spec_0005": set(),
                },
                "peptides": {"all": set(), "safe": set(), "safe_num_specs": set()},
                "spectra": {"all": set()},
                # protein_dict in contrast is a nested dict with protein/protein_group --> peptide sequence --> spectral information
                # (containing lists of 'spec_title', 'bayes_pep', modifications', 'charge', 'psm_q_value', 'start_stop')
                "protein_dict": {},
            }
        }
        fdr_dict = {
            "peptides_seq_level": {},
            "peptides_psm_level": {},
            "peptides_seq_level_2specs": {},
            "glycopeptides_psm_level": {},
            "glycopeptides_seq_level": {},
            "glycopeptides_seq_level_2specs": {},
            "proteins_seq_level": {},
            "proteins_psm_level": {},
            "proteins_seq_level_2specs": {},
        }

    result_file_list = []
    for PRIDE_ID in datasets.keys():
        if skip_old is True and old_exists is True and PRIDE_ID in results_dict:
            continue
        print("reading:", PRIDE_ID)
        instrument = datasets[PRIDE_ID]["instrument"]
        results_dict["all"]["instrument"].add(instrument)
        lab = datasets[PRIDE_ID]["lab"]
        results_dict["all"]["lab"].add(lab)
        results_dict["all"]["num_spectra"] += datasets[PRIDE_ID]["num_spectra"]
        if PRIDE_ID not in results_dict.keys():
            results_dict[PRIDE_ID] = {
                "num_spectra": datasets[PRIDE_ID]["num_spectra"],
                "instrument": instrument,
                "lab": lab,
                "protein_groups": {
                    "all": set(),
                    "safe_psm": set(),
                    "safe_seq": set(),
                    "safe_seq_num_spec": set(),
                    "safe_seq_num_spec_0005": set(),
                },
                "proteins": {
                    "all": set(),
                    "safe_psm": set(),
                    "safe_seq": set(),
                    "safe_seq_num_spec": set(),
                    "safe_seq_num_spec_0005": set(),
                },
                "peptides": {"all": set(), "safe": set(), "safe_num_specs": set()},
                "spectra": {"all": set()},
                "protein_dict": {},
            }

        # collect proteins, peptides and corresponding spectrum_titles
        if datasets[PRIDE_ID]["folders"] != "":
            PRIDE_folder = os.path.join(PRIDE_ID, datasets[PRIDE_ID]["folders"])
        else:
            PRIDE_folder = PRIDE_ID
        merged_file = os.path.join(
            PRIDE_folder,
            datasets[PRIDE_ID]["result_file"],
        )
        result_file_list.append(merged_file)
        protein_ids = set()
        protein_groups = set()
        with open(merged_file, "r") as in_file:
            result_csv = csv.DictReader(in_file)
            for line_dict in result_csv:
                seq = line_dict["Sequence"]  # + line_dict['Modifications']
                mod = line_dict["Modifications"]
                mods = []
                # In contrast to the original ArcPP analysis, modifications are taken into account
                # except for optional modifications that depend on the sample preparation.
                # In the following, commented out sections indicate the use of "seq"
                # that has now been changed to "seq_mod".
                for m in line_dict["Modifications"].split(";"):
                    if "iTRAQ4plex" in m or "Label:" in m or "Oxidation" in m:
                        continue
                    mods.append(m)
                charge = line_dict["Charge"]
                # seq_mod = '{0}#{1}'.format(seq, mod)
                seq_mod = "{0}#{1}".format(seq, ";".join(mods))
                seq_length = len(seq)
                spec_title = line_dict["Spectrum Title"]
                sample = spec_title.split(".")[0]
                is_decoy = line_dict["Is decoy"]
                prot = line_dict["Protein ID"]
                start = line_dict["Sequence Start"]
                stop = line_dict["Sequence Stop"]
                pre = line_dict["Sequence Pre AA"]
                post = line_dict["Sequence Post AA"]
                psm_q_value = float(line_dict["combined PEP"])
                bayes_pep = float(line_dict["Bayes PEP"])
                if psm_q_value <= 0.01:
                    if seq_length not in fdr_dict["peptides_psm_level"].keys():
                        fdr_dict["peptides_psm_level"][seq_length] = {}
                    if seq not in fdr_dict["peptides_psm_level"][seq_length].keys():
                        fdr_dict["peptides_psm_level"][seq_length][seq] = (
                            psm_q_value,
                            is_decoy,
                        )
                    elif (
                        psm_q_value < fdr_dict["peptides_psm_level"][seq_length][seq][0]
                    ):
                        fdr_dict["peptides_psm_level"][seq_length][seq] = (
                            psm_q_value,
                            is_decoy,
                        )
                    if (
                        "Hex" in mod
                        and seq_length not in fdr_dict["glycopeptides_psm_level"].keys()
                    ):
                        fdr_dict["glycopeptides_psm_level"][seq_length] = {}
                    if (
                        "Hex" in mod
                        and seq_mod
                        not in fdr_dict["glycopeptides_psm_level"][seq_length].keys()
                    ):
                        fdr_dict["glycopeptides_psm_level"][seq_length][seq_mod] = (
                            psm_q_value,
                            is_decoy,
                        )
                    elif (
                        "Hex" in mod
                        and psm_q_value
                        < fdr_dict["glycopeptides_psm_level"][seq_length][seq_mod][0]
                    ):
                        fdr_dict["glycopeptides_psm_level"][seq_length][seq_mod] = (
                            psm_q_value,
                            is_decoy,
                        )
                else:
                    print(
                        "Results should be filtered by combined PEP <= 1% (but should contain targets and decoys)"
                    )
                    sys.exit(1)

                # differentiate between protein groups and proteins
                # and remove contaminants
                if len(prot.split("<|>")) > 1:
                    contaminants = True
                    for p in prot.split("<|>"):
                        prot_id = p.split(" ")[0]
                        if "HVO" not in prot_id:
                            continue
                        else:
                            contaminants = False
                    # contaminants = False
                    if contaminants is False and is_decoy == "false":
                        results_dict[PRIDE_ID]["protein_groups"]["all"].add(
                            line_dict["Protein ID"]
                        )
                        results_dict[PRIDE_ID]["peptides"]["all"].add(seq_mod)
                        # results_dict[PRIDE_ID]['peptides']['all'].add(seq)
                        results_dict[PRIDE_ID]["spectra"]["all"].add(spec_title)
                else:
                    contaminants = False
                    prot_id = prot.split(" ")[0]
                    if "HVO" not in prot_id:
                        contaminants = True
                    if contaminants is False and is_decoy == "false":
                        results_dict[PRIDE_ID]["proteins"]["all"].add(
                            line_dict["Protein ID"]
                        )
                        results_dict[PRIDE_ID]["peptides"]["all"].add(seq_mod)
                        # results_dict[PRIDE_ID]['peptides']['all'].add(seq)
                        results_dict[PRIDE_ID]["spectra"]["all"].add(spec_title)

                # add info to protein_dict
                if prot not in results_dict[PRIDE_ID]["protein_dict"].keys():
                    results_dict[PRIDE_ID]["protein_dict"][prot] = {}
                # if seq not in results_dict[PRIDE_ID]['protein_dict'][prot].keys():
                #     results_dict[PRIDE_ID]['protein_dict'][prot][seq] = {
                if seq_mod not in results_dict[PRIDE_ID]["protein_dict"][prot].keys():
                    results_dict[PRIDE_ID]["protein_dict"][prot][seq_mod] = {
                        "spec_title": [],
                        "bayes_pep": [],
                        "modifications": [],
                        "charge": [],
                        "psm_q_value": [],
                        "start_stop": (start, stop, pre, post),
                    }
                results_dict[PRIDE_ID]["protein_dict"][prot][seq_mod][
                    "spec_title"
                ].append(spec_title)
                results_dict[PRIDE_ID]["protein_dict"][prot][seq_mod][
                    "bayes_pep"
                ].append(bayes_pep)
                results_dict[PRIDE_ID]["protein_dict"][prot][seq_mod][
                    "psm_q_value"
                ].append(psm_q_value)
                results_dict[PRIDE_ID]["protein_dict"][prot][seq_mod][
                    "modifications"
                ].append(mod)
                results_dict[PRIDE_ID]["protein_dict"][prot][seq_mod]["charge"].append(
                    charge
                )
                # results_dict[PRIDE_ID]['protein_dict'][prot][seq]['spec_title'].append(spec_title)
                # results_dict[PRIDE_ID]['protein_dict'][prot][seq]['bayes_pep'].append(bayes_pep)
                # results_dict[PRIDE_ID]['protein_dict'][prot][seq]['psm_q_value'].append(psm_q_value)
                # results_dict[PRIDE_ID]['protein_dict'][prot][seq]['modifications'].append(mod)
                # results_dict[PRIDE_ID]['protein_dict'][prot][seq]['charge'].append(charge)

        # merge identifications from each dataset into "all"
        for level in ["protein_groups", "proteins", "peptides", "spectra"]:
            results_dict["all"][level]["all"] |= results_dict[PRIDE_ID][level]["all"]
        for prot in results_dict[PRIDE_ID]["protein_dict"].keys():
            if prot not in results_dict["all"]["protein_dict"].keys():
                results_dict["all"]["protein_dict"][prot] = {"datasets": set()}
            results_dict["all"]["protein_dict"][prot]["datasets"].add(PRIDE_ID)
            for seq in results_dict[PRIDE_ID]["protein_dict"][prot].keys():
                start_stop = results_dict[PRIDE_ID]["protein_dict"][prot][seq][
                    "start_stop"
                ]
                if seq not in results_dict["all"]["protein_dict"][prot].keys():
                    results_dict["all"]["protein_dict"][prot][seq] = {
                        "spec_title": [],
                        "bayes_pep": [],
                        "modifications": [],
                        "charge": [],
                        "psm_q_value": [],
                        "start_stop": start_stop,
                    }
                for k, v in results_dict[PRIDE_ID]["protein_dict"][prot][seq].items():
                    if k == "start_stop":
                        continue
                    results_dict["all"]["protein_dict"][prot][seq][k].extend(v)

    # Calculate q-values
    # peptides first, then proteins
    for PRIDE_ID in results_dict.keys():
        # generate input dict for q_value calculation function
        seq_q_value_dict = {}
        for prot in results_dict[PRIDE_ID]["protein_dict"].keys():
            for seq_mod in results_dict[PRIDE_ID]["protein_dict"][prot].keys():
                if seq_mod == "datasets":
                    continue
                seq_length = len(seq_mod.split("#")[0])
                min_bayes_pep = min(
                    results_dict[PRIDE_ID]["protein_dict"][prot][seq_mod]["bayes_pep"]
                )
                if seq_length not in seq_q_value_dict.keys():
                    seq_q_value_dict[seq_length] = {}
                if "decoy_" in prot:
                    is_decoy = True
                else:
                    is_decoy = False
                seq_q_value_dict[seq_length][seq_mod] = {
                    "Bayes PEP": min_bayes_pep,
                    "Is decoy": is_decoy,
                }
            # for seq in results_dict[PRIDE_ID]['protein_dict'][prot].keys():
            #     if seq == 'datasets':
            #         continue
            #     seq_length = len(seq)
            #     min_bayes_pep = min(
            #         results_dict[PRIDE_ID]['protein_dict'][prot][seq]['bayes_pep']
            #     )
            #     if seq_length not in seq_q_value_dict.keys():
            #         seq_q_value_dict[seq_length] = {}
            #     if 'decoy_' in prot:
            #         is_decoy = True
            #     else:
            #         is_decoy = False
            #     seq_q_value_dict[seq_length][seq] = {
            #         'Bayes PEP' : min_bayes_pep,
            #         'Is decoy' : is_decoy,
            #     }

        print("calculating q-values on peptide level")
        seq_calc_q_value_dict = calculate_q_value_by_group(
            seq_q_value_dict, sliding=False
        )

        # read results from peptide q_value calc, at the same time
        # generate input dict for proteins for q_value calculation function
        prot_q_value_dict = {"seq_level": {}, "psm_level": {}}
        for prot in results_dict[PRIDE_ID]["protein_dict"].keys():
            contaminants = False
            prot_id = prot.split(" ")[0]
            if "HVO" not in prot_id:
                contaminants = True
            if "decoy_" in prot:
                is_decoy = True
            else:
                is_decoy = False
            for seq in results_dict[PRIDE_ID]["protein_dict"][prot].keys():
                if seq == "datasets":
                    continue
                # seq_length = len(seq)
                seq_length = len(seq.split("#")[0])
                seq_q_value = seq_calc_q_value_dict[seq_length][seq]["combined PEP"]
                results_dict[PRIDE_ID]["protein_dict"][prot][seq][
                    "seq_q_value"
                ] = seq_q_value

                if seq_q_value <= SEQ_Q_VALUE_THRESHOLD:
                    if PRIDE_ID == "all":
                        if "Hex" in seq:
                            if (
                                seq_length
                                not in fdr_dict["glycopeptides_seq_level"].keys()
                            ):
                                fdr_dict["glycopeptides_seq_level"][seq_length] = {}
                            fdr_dict["glycopeptides_seq_level"][seq_length][seq] = (
                                seq_q_value,
                                is_decoy,
                            )
                        else:
                            if seq_length not in fdr_dict["peptides_seq_level"].keys():
                                fdr_dict["peptides_seq_level"][seq_length] = {}
                            fdr_dict["peptides_seq_level"][seq_length][seq] = (
                                seq_q_value,
                                is_decoy,
                            )
                    counts = len(
                        set(
                            results_dict[PRIDE_ID]["protein_dict"][prot][seq][
                                "spec_title"
                            ]
                        )
                    )
                    if is_decoy is False and contaminants is False:
                        results_dict[PRIDE_ID]["peptides"]["safe"].add(seq)
                        if counts >= num_specs:
                            results_dict[PRIDE_ID]["peptides"]["safe_num_specs"].add(seq)
                            if PRIDE_ID == "all":
                                if "Hex" in seq:
                                    if (
                                        seq_length
                                        not in fdr_dict[
                                            "glycopeptides_seq_level_2specs"
                                        ].keys()
                                    ):
                                        fdr_dict["glycopeptides_seq_level_2specs"][
                                            seq_length
                                        ] = {}
                                    fdr_dict["glycopeptides_seq_level_2specs"][
                                        seq_length
                                    ][seq] = (seq_q_value, is_decoy)
                                else:
                                    if (
                                        seq_length
                                        not in fdr_dict[
                                            "peptides_seq_level_2specs"
                                        ].keys()
                                    ):
                                        fdr_dict["peptides_seq_level_2specs"][
                                            seq_length
                                        ] = {}
                                    fdr_dict["peptides_seq_level_2specs"][seq_length][
                                        seq
                                    ] = (seq_q_value, is_decoy)
                    min_bayes_pep = min(
                        results_dict[PRIDE_ID]["protein_dict"][prot][seq]["bayes_pep"]
                    )
                    if min_bayes_pep == 0.0:
                        min_bayes_pep = np.nextafter(0, 1)
                    log_seq_bayes = math.log10(min_bayes_pep)
                    if prot not in prot_q_value_dict["seq_level"].keys():
                        prot_q_value_dict["seq_level"][prot] = {
                            "Bayes PEP": log_seq_bayes,
                            "Is decoy": is_decoy,
                        }
                    else:
                        prot_q_value_dict["seq_level"][prot][
                            "Bayes PEP"
                        ] += log_seq_bayes

                for bayes_pep in results_dict[PRIDE_ID]["protein_dict"][prot][seq][
                    "bayes_pep"
                ]:
                    if bayes_pep == 0.0:
                        bayes_pep = np.nextafter(0, 1)
                    log_psm_bayes = math.log10(bayes_pep)
                    if prot not in prot_q_value_dict["psm_level"].keys():
                        prot_q_value_dict["psm_level"][prot] = {
                            "Bayes PEP": log_seq_bayes,
                            "Is decoy": is_decoy,
                        }
                    else:
                        prot_q_value_dict["psm_level"][prot][
                            "Bayes PEP"
                        ] += log_seq_bayes

        print("calculating q-values on protein level")
        prot_calc_q_value_dict = calculate_q_value_by_group(
            prot_q_value_dict, sliding=False, picked_fdr=True
        )

        # read results from protein q_value calc
        for prot in results_dict[PRIDE_ID]["protein_dict"].keys():
            contaminants = False
            prot_id = prot.split(" ")[0]
            if "HVO" not in prot_id:
                contaminants = True
            if "decoy_" in prot:
                is_decoy = True
            else:
                is_decoy = False
            for level in ["psm_level", "seq_level"]:
                if prot in prot_calc_q_value_dict[level].keys():
                    prot_q_value = prot_calc_q_value_dict[level][prot]["combined PEP"]
                    prot_bayes_pep = prot_calc_q_value_dict[level][prot]["Bayes PEP"]
                else:
                    prot_q_value = 1
                    prot_bayes_pep = 1
                # count number of spectra for each prot (for seq FDR > 1%)
                # collect samples for simple protein inference model
                counts = 0
                samples = set()
                for seq in results_dict[PRIDE_ID]["protein_dict"][prot].keys():
                    if seq in [
                        "datasets",
                        "prot_q_value_seq",
                        "prot_q_value_psm",
                        "samples",
                    ]:
                        continue
                    if (
                        results_dict[PRIDE_ID]["protein_dict"][prot][seq]["seq_q_value"]
                        > 0.01
                    ):
                        continue
                    psm_set = set(
                        results_dict[PRIDE_ID]["protein_dict"][prot][seq]["spec_title"]
                    )
                    counts += len(psm_set)
                    for psm in psm_set:
                        ms_filename = ".".join(psm.split(".")[:-3])
                        samples.add(ms_filename2sample.get(ms_filename, ms_filename))

                if PRIDE_ID == "all":
                    if level == "seq_level":
                        fdr_dict["proteins_seq_level"][prot] = (prot_bayes_pep, is_decoy)
                        if counts >= num_specs:
                            fdr_dict["proteins_seq_level_2specs"][prot] = (
                                prot_bayes_pep,
                                is_decoy,
                            )
                    else:
                        fdr_dict["proteins_psm_level"][prot] = (prot_bayes_pep, is_decoy)
                if prot_q_value <= 0.01 and is_decoy is False and contaminants is False:
                    if level == "seq_level":
                        if len(prot.split("<|>")) > 1:
                            results_dict[PRIDE_ID]["protein_groups"]["safe_seq"].add(
                                prot
                            )
                            if counts >= num_specs:
                                results_dict[PRIDE_ID]["protein_groups"][
                                    "safe_seq_num_spec"
                                ].add(prot)
                                if prot_q_value <= PROT_Q_VALUE_THRESHOLD:
                                    results_dict[PRIDE_ID]["protein_groups"][
                                        "safe_seq_num_spec_0005"
                                    ].add(prot)
                        else:
                            results_dict[PRIDE_ID]["proteins"]["safe_seq"].add(prot)
                            if counts >= num_specs:
                                results_dict[PRIDE_ID]["proteins"][
                                    "safe_seq_num_spec"
                                ].add(prot)
                                if prot_q_value <= PROT_Q_VALUE_THRESHOLD:
                                    results_dict[PRIDE_ID]["proteins"][
                                        "safe_seq_num_spec_0005"
                                    ].add(prot)
                    elif counts >= num_specs:
                        if len(prot.split("<|>")) > 1:
                            results_dict[PRIDE_ID]["protein_groups"]["safe_psm"].add(
                                prot
                            )
                        else:
                            results_dict[PRIDE_ID]["proteins"]["safe_psm"].add(prot)
                if level == "seq_level":
                    results_dict[PRIDE_ID]["protein_dict"][prot][
                        "prot_q_value_seq"
                    ] = prot_q_value
                else:
                    results_dict[PRIDE_ID]["protein_dict"][prot][
                        "prot_q_value_psm"
                    ] = prot_q_value
                results_dict[PRIDE_ID]["protein_dict"][prot]["samples"] = samples
        print(
            "Number of confident protein identifications for {0}: {1}".format(
                PRIDE_ID,
                len(results_dict[PRIDE_ID]["proteins"]["safe_seq_num_spec_0005"]),
            )
        )

    # save results in a pkl
    pickle.dump(results_dict, open(pkl_name, "wb"))
    print("pickled results: ", pkl_name)

    pickle.dump(fdr_dict, open(fdr_pkl_name, "wb"))
    print("pickled fdr_dict: ", fdr_pkl_name)


def calculate_q_value_by_group(input_dict, sliding=True, picked_fdr=False):
    if picked_fdr:
        dict_for_calc = {}
        for group in input_dict.keys():
            dict_for_calc[group] = {}
            for prot in input_dict[group].keys():
                is_decoy = input_dict[group][prot]["Is decoy"]
                if is_decoy:
                    decoy = prot
                    decoy_score = input_dict[group][prot]["Bayes PEP"]
                    target = prot.replace("decoy_", "")
                    if target in input_dict[group].keys():
                        target_score = input_dict[group][target]["Bayes PEP"]
                    else:
                        target_score = 1
                else:
                    target = prot
                    target_score = input_dict[group][prot]["Bayes PEP"]
                    decoy = "decoy_" + prot
                    if decoy in input_dict[group].keys():
                        decoy_score = input_dict[group][decoy]["Bayes PEP"]
                    else:
                        decoy_score = 1
                if target_score < decoy_score:
                    dict_for_calc[group][target] = {
                        "Bayes PEP": target_score,
                        "Is decoy": input_dict[group][target]["Is decoy"],
                    }
                elif decoy_score < target_score:
                    dict_for_calc[group][decoy] = {
                        "Bayes PEP": decoy_score,
                        "Is decoy": input_dict[group][decoy]["Is decoy"],
                    }
                else:
                    print("Target and decoy have the same score, continuing with both")
                    dict_for_calc[group][decoy] = {
                        "Bayes PEP": decoy_score,
                        "Is decoy": input_dict[group][decoy]["Is decoy"],
                    }
                    dict_for_calc[group][target] = {
                        "Bayes PEP": target_score,
                        "Is decoy": input_dict[group][target]["Is decoy"],
                    }
    else:
        dict_for_calc = input_dict

    for group in dict_for_calc.keys():
        # print('calculating q-values for:', group)
        psms_sorted_by_bayes_pep = sorted(
            dict_for_calc[group].items(), key=lambda kv: kv[1]["Bayes PEP"]
        )
        sorted_decoy_bools = [
            kv_tuple[1]["Is decoy"] for kv_tuple in psms_sorted_by_bayes_pep
        ]

        if sliding:
            for i, (n_decoys, __, current_win_size) in enumerate(
                sliding_window(sorted_decoy_bools, window_size=249)
            ):

                n_false_positives = 2 * n_decoys
                group_q_value = n_false_positives / current_win_size
                current_seq = psms_sorted_by_bayes_pep[i][0]
                dict_for_calc[group][current_seq]["combined PEP"] = group_q_value

            assert i + 1 == len(sorted_decoy_bools), (
                "sliding_window() "
                "did not return a sliding window for all PSMs! This "
                "should never happen!"
            )

        else:
            q_value_top2bottom = []
            decoys = 0
            targets = 0
            for i, decoy_bool in enumerate(sorted_decoy_bools):
                if decoy_bool is True:
                    decoys += 1
                if decoy_bool is False:
                    targets += 1
                if targets == 0:
                    q_value_top2bottom.append(1)
                else:
                    q_value_top2bottom.append(decoys / targets)

            min_obs_so_far = 1
            q_value_bottom2top = []
            for j, qvalue in enumerate(reversed(q_value_top2bottom)):
                if qvalue < min_obs_so_far:
                    min_obs_so_far = qvalue
                q_value_bottom2top.append(min_obs_so_far)

            for k, group_q_value in enumerate(reversed(q_value_bottom2top)):
                current_seq = psms_sorted_by_bayes_pep[k][0]
                dict_for_calc[group][current_seq]["combined PEP"] = group_q_value

    return dict_for_calc


def sliding_window(elements, window_size=249, flexible=True):
    """
    Sliding window generator.
    Gives you sliding window functionality without using container
    types (list, deque etc.) to speed it up. Only works for lists of
    numbers. Yields the sum of all numbers in the sliding window
    (= the number of decoys in the sliding window in our case), the
    central number of the sliding window (required for the test only),
    and the current length of the sliding window (= total number of
    PSMs in the sliding window). Used for PEP calculation:
    PEP_of_PSM = (n_decoys_in_window * 2) / n_total_PSMs_in_window
    """
    if flexible:
        window_size = adjust_window_size(window_size, len(elements))

    if window_size % 2 == 0:
        print(
            "Warning! Window size must be uneven (to determine a "
            "central value). Adjusted window size from {0} to {1}"
            ".".format(window_size, window_size + 1)
        )
        window_size += 1

    half_window_size = int((window_size - 1) / 2)

    start_gen, stop_gen = itertools.tee(elements)  # get 2 generators, one that tells
    # us which number to subtract from the back and another one that tells us
    # which number to add at the front of the sliding window

    n_decoys = 0  # keep track of the number of decoys in current sliding window
    current_win_size = 0  # keep track of current window size
    previous_start_i, previous_stop_i = 0, 0  # remember where our sliding window
    # was one iteration earlier (start and stop positions of the sliding window)

    for center_i, center_value in enumerate(elements):
        start_i = center_i - half_window_size
        if start_i < 0:
            start_i = 0
        stop_i = center_i + half_window_size + 1

        if start_i != previous_start_i:
            n_decoys -= next(start_gen)
            current_win_size -= 1

        if stop_i != previous_stop_i:
            for i in range(stop_i - previous_stop_i):
                try:
                    n_decoys += next(stop_gen)
                    current_win_size += 1
                except StopIteration:
                    break  # cause StopIteration silently ends for-loops, will be fixed in py3.6 :)

        previous_start_i, previous_stop_i = start_i, stop_i
        yield n_decoys, center_value, current_win_size


def adjust_window_size(desired_window_size, iter_len, minimum=29):
    """
    Taken from Ursgal combined_PEP function

    Dynamically adjusts the sliding window size depending on the total
    length of values. When there are few values (below 1/5 of the
    window size), the window size is decreased.

    """
    if desired_window_size < iter_len // 5:
        adjusted_window_size = desired_window_size
    else:
        adjusted_window_size = desired_window_size // 5
        if adjusted_window_size < minimum:
            adjusted_window_size = minimum
    if adjusted_window_size != desired_window_size:

        print(
            "Adjusted window size from {0} to {1} because there "
            "are only {2} PSMs.".format(
                desired_window_size, adjusted_window_size, iter_len
            )
        )
    return adjusted_window_size


if __name__ == "__main__":
    if sys.argv[1] == "True":
        skip = True
    else:
        skip = False
    num_specs = int(sys.argv[2])
    main(os.path.dirname(sys.argv[0]), skip_old=skip, num_specs=num_specs)
