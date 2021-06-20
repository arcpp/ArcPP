#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import os
import csv
import pickle
import pprint
import glob
import re

# import plotly
# import plotly.graph_objs as go
from collections import defaultdict as ddict
from itertools import combinations

PRIDE_ids = [
    "PXD021874",
    "PXD007061",
    "PXD013046",
    "PXD011012",
    "PXD006877",
    "PXD011218",
    "PXD009116",
    "PXD011056",
    "PXD011015",
    "PXD011050",
    "PXD014974",
    "PXD010824",
]

MOD_POS_PATTERN = re.compile(r"(?P<modname>.*):(?P<pos>[0-9]*)$")


def main(input_file=None, arcpp_pep_file=None):
    """
    Post-process glycopeptide identifications in order to ensure
    that glycopeptide identifications are substantiated by at least
    two PSMs and two replicates.
    Furthermore, split results into N-glycans, non-canonical N-glycans and O-glycans.

    Usage:
        python <script_name.py> <path_to_glycopeptide_PSMs> <path_to_ArcPP_peptides_file>
    """
    params = {
        "psm_defining_colnames": [
            "Spectrum Title",
            "Sequence",
            "Modifications",
            "Charge",
            "Is decoy",
        ],
    }
    uc = ursgal.UController(params=params)

    arcpp_glycopeps = {}
    with open(arcpp_pep_file, "r") as arcpp_in:
        arcpp_csv = csv.DictReader(arcpp_in)
        for line_dict in arcpp_csv:
            glycopep, glycan_type = line_to_pep_unimod_glyc(line_dict)
            if glycopep not in arcpp_glycopeps.keys():
                arcpp_glycopeps[glycopep] = []
            arcpp_glycopeps[glycopep].append(line_dict)

    replicate_lookup = {}
    for pride in PRIDE_ids:
        replicate_lookup[pride] = {}
        file_description = os.path.join(
            "file_descriptions", "{0}_file_descriptions.csv".format(pride)
        )
        with open(file_description, "r") as descr_in:
            descr_csv = csv.DictReader(descr_in)
            for line_dict in descr_csv:
                strain = line_dict["Strain"]
                file_name = line_dict["Raw file name"].split(".")[0]
                replicate = line_dict["Replicate"]
                replicate_lookup[pride][file_name] = {
                    "strain": strain,
                    "rep": replicate,
                }

    true_n_glycopeps = ddict(dict)
    non_standard_n_glycopeps = ddict(dict)
    o_glycopeps = ddict(dict)
    all_strain = set()
    with open(input_file, "r") as glyco_in:
        glyco_csv = csv.DictReader(glyco_in)
        fieldnames = glyco_csv.fieldnames
        for line_dict in glyco_csv:
            protein = line_dict["Protein ID"]
            if protein.startswith("sp|"):
                continue
            # peptide = line_dict['Sequence']
            spec_id = line_dict["Spectrum Title"]
            file_name = line_dict["Spectrum Title"].split(".")[0]
            dataset = line_dict["Dataset"]
            strain = replicate_lookup[dataset][file_name]["strain"]
            rep = replicate_lookup[dataset][file_name]["rep"]
            glycopep, glycan_type = line_to_pep_unimod_glyc(line_dict)
            if strain not in true_n_glycopeps.keys():
                true_n_glycopeps[strain] = {}
                all_strain.add(strain)
            if "n_glycan" in glycan_type:
                if glycopep not in true_n_glycopeps[strain].keys():
                    true_n_glycopeps[strain][glycopep] = {
                        "frag_ions": set(),
                        "specs": set(),
                        "reps": set(),
                        "line_dicts": [],
                    }
                true_n_glycopeps[strain][glycopep]["specs"].add(spec_id)
                true_n_glycopeps[strain][glycopep]["reps"].add("#".join([dataset, rep]))
                true_n_glycopeps[strain][glycopep]["line_dicts"].append(line_dict)
                true_n_glycopeps[strain][glycopep]["frag_ions"].add(
                    line_dict["MS2 Glycopep Frag Ions Present"]
                )
            elif "o_glycan" in glycan_type and len(glycan_type) == 1:
                if strain not in o_glycopeps.keys():
                    o_glycopeps[strain] = {}
                    all_strain.add(strain)
                if glycopep not in o_glycopeps[strain].keys():
                    o_glycopeps[strain][glycopep] = {
                        "frag_ions": set(),
                        "specs": set(),
                        "reps": set(),
                        "line_dicts": [],
                    }
                o_glycopeps[strain][glycopep]["specs"].add(spec_id)
                o_glycopeps[strain][glycopep]["reps"].add("#".join([dataset, rep]))
                o_glycopeps[strain][glycopep]["line_dicts"].append(line_dict)
                o_glycopeps[strain][glycopep]["frag_ions"].add(
                    line_dict["MS2 Glycopep Frag Ions Present"]
                )
            elif "true_non_standard_n_glycan" in glycan_type:
                if strain not in non_standard_n_glycopeps.keys():
                    non_standard_n_glycopeps[strain] = {}
                    all_strain.add(strain)
                if glycopep not in non_standard_n_glycopeps[strain].keys():
                    non_standard_n_glycopeps[strain][glycopep] = {
                        "frag_ions": set(),
                        "specs": set(),
                        "reps": set(),
                        "line_dicts": [],
                    }
                non_standard_n_glycopeps[strain][glycopep]["specs"].add(spec_id)
                non_standard_n_glycopeps[strain][glycopep]["reps"].add(
                    "#".join([dataset, rep])
                )
                non_standard_n_glycopeps[strain][glycopep]["line_dicts"].append(
                    line_dict
                )
                non_standard_n_glycopeps[strain][glycopep]["frag_ions"].add(
                    line_dict["MS2 Glycopep Frag Ions Present"]
                )
            else:
                print(glycan_type)

    count_true_n_glycopeps = set()
    count_true_n_glycopeps_arcpp = set()
    count_o_glycopeps = set()
    count_o_glycopeps_arcpp = set()
    count_non_standard_glycopeps = set()
    count_non_standard_glycopeps_arcpp = set()
    output_line_dicts_n = []
    output_line_dicts_non_standard_n = []
    output_line_dicts_o = []
    for strain in all_strain:
        print(strain)
        for glycopep in true_n_glycopeps[strain].keys():
            if "True" not in true_n_glycopeps[strain][glycopep]["frag_ions"]:
                continue
            if len(true_n_glycopeps[strain][glycopep]["specs"]) < 2:
                continue
            if len(true_n_glycopeps[strain][glycopep]["reps"]) < 2:
                continue
            count_true_n_glycopeps.add(glycopep)
            if glycopep not in arcpp_glycopeps.keys():
                continue
            count_true_n_glycopeps_arcpp.add(glycopep)
            output_line_dicts_n.extend(true_n_glycopeps[strain][glycopep]["line_dicts"])

        for glycopep in o_glycopeps[strain].keys():
            if "True" not in o_glycopeps[strain][glycopep]["frag_ions"]:
                continue
            if len(o_glycopeps[strain][glycopep]["specs"]) < 2:
                continue
            if len(o_glycopeps[strain][glycopep]["reps"]) < 2:
                continue
            count_o_glycopeps.add(glycopep)
            if glycopep not in arcpp_glycopeps.keys():
                continue
            count_o_glycopeps_arcpp.add(glycopep)
            output_line_dicts_o.extend(o_glycopeps[strain][glycopep]["line_dicts"])

        for glycopep in non_standard_n_glycopeps[strain].keys():
            if "True" not in non_standard_n_glycopeps[strain][glycopep]["frag_ions"]:
                continue
            if len(non_standard_n_glycopeps[strain][glycopep]["specs"]) < 2:
                continue
            if len(non_standard_n_glycopeps[strain][glycopep]["reps"]) < 2:
                continue
            count_non_standard_glycopeps.add(glycopep)
            if glycopep not in arcpp_glycopeps.keys():
                continue
            count_non_standard_glycopeps_arcpp.add(glycopep)
            output_line_dicts_non_standard_n.extend(
                non_standard_n_glycopeps[strain][glycopep]["line_dicts"]
            )

    print(
        """
        True N-glyco: {0}
        True N-glyco ArcPP: {1}

        O-glyco: {2}
        O-glyco ArcPP: {3}

        Non-standard N-glyco: {4}
        Non-standard N-glyco ArcPP: {5}
    """.format(
            len(count_true_n_glycopeps),
            len(count_true_n_glycopeps_arcpp),
            len(count_o_glycopeps),
            len(count_o_glycopeps_arcpp),
            len(count_non_standard_glycopeps),
            len(count_non_standard_glycopeps_arcpp),
        )
    )

    csv_kwargs = {}
    if sys.platform == "win32":
        csv_kwargs["lineterminator"] = "\n"
    else:
        csv_kwargs["lineterminator"] = "\r\n"
    csv_out_name = "ArcPP_N_glyco_filtered_peptides_2rep.csv"
    with open(csv_out_name, "w") as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for out_dict in output_line_dicts_n:
            csv_writer.writerow(out_dict)

    csv_out_name = "ArcPP_only_O_glyco_filtered_peptides_2rep.csv"
    with open(csv_out_name, "w") as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for out_dict in output_line_dicts_o:
            csv_writer.writerow(out_dict)

    csv_out_name = "ArcPP_only_non_canonical_n_glyco_filtered_peptides_2rep.csv"
    with open(csv_out_name, "w") as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for out_dict in output_line_dicts_non_standard_n:
            csv_writer.writerow(out_dict)


def line_to_pep_unimod_glyc(line_dict):
    """
    Organize unimods, collect glycans
    """
    mod_pattern = re.compile(r""":(?P<pos>[0-9]*$)""")
    monosaccharide_compositions = (
        ursgal.chemical_composition_kb.monosaccharide_compositions
    )
    new_unimods = []
    glycan_composition = ddict(int)
    peptide = line_dict["Sequence"]
    glycan_type = set()
    for unimod in line_dict["Modifications"].split(";"):
        match = mod_pattern.search(unimod)
        if match is not None:
            pos = int(match.group("pos"))
            mod = unimod[: match.start()]
            length = 0
            glycomod = False
            if mod in monosaccharide_compositions.keys():
                glycan_composition[mod] += 1
                glycomod = True
            for monosacch in monosaccharide_compositions.keys():
                if monosacch in mod:
                    for k, v in glycan_to_dict(mod).items():
                        glycan_composition[k] += v
                        length += v
                    glycomod = True
                    break
            if glycomod is False:
                if "iTRAQ4plex" in mod or "Label:" in mod or "Oxidation" in mod:
                    continue
                new_unimods.append(unimod)
            else:
                if line_dict["Sequence"][pos - 1] != "N":
                    glycan_type.add("o_glycan")
                elif (
                    len(glycan_composition.keys()) == 1
                    and list(glycan_composition.keys())[0] == "Hex"
                ):
                    glycan_type.add("non_standard_n_glycan")
                elif length <= 3:
                    glycan_type.add("non_standard_n_glycan")
                else:
                    n_glycan = False
                    for pos, aa in enumerate(line_dict["Sequence"]):
                        if aa == "N":
                            if pos + 1 == len(line_dict["Sequence"]):
                                # glycan_type.add('non_standard_n_glycan')
                                continue
                            elif pos + 1 >= len(line_dict["Sequence"]) - 2:
                                if (
                                    "S" not in line_dict["Sequence Post AA"]
                                    and "T" not in line_dict["Sequence Post AA"]
                                ):
                                    # glycan_type.add('non_standard_n_glycan')
                                    continue
                                else:
                                    glycan_type.add("n_glycan")
                                    n_glycan = True
                            else:
                                if line_dict["Sequence"][pos + 2] not in ["S", "T"]:
                                    # glycan_type.add('non_standard_n_glycan')
                                    continue
                                else:
                                    glycan_type.add("n_glycan")
                                    n_glycan = True
                    if n_glycan is False:
                        glycan_type.add("true_non_standard_n_glycan")
    glycan = line_dict.get("Glycan", "")
    if glycan != "":
        for k, v in glycan_to_dict(glycan).items():
            glycan_composition[k] += v
    # if new_unimods != []:
    peptide_unimod_glycan = "{0}#{1}#{2}".format(
        peptide,
        ";".join(new_unimods),
        "".join(["{0}({1})".format(k, v) for k, v in glycan_composition.items()]),
    )
    return peptide_unimod_glycan, glycan_type


def glycan_to_dict(glycan):
    """
    Converts a glycan (unimod style: Hex(2)HexNAc(5)) into a dict
    with key=monosaccharide and value=count
    """
    pattern = re.compile(r"""(?P<monosacch>[A-z0-9]*)(?P<count>\([A-z0-9]*\))""")
    glycan_dict = {}
    for glyc_match in pattern.finditer(glycan):
        monosacch = glyc_match.group("monosacch")
        if monosacch == "End":
            # count = glyc_match.group('count').strip('(').strip(')')
            # if count == 'None':
            #     count = None
            continue
        elif glyc_match.group("count") == "":
            count = 1
        else:
            count = int(glyc_match.group("count").strip("(").strip(")"))
        if monosacch not in glycan_dict.keys():
            glycan_dict[monosacch] = 0
        glycan_dict[monosacch] += count
    # sp_tuple = tuple(sorted(glycan_dict.items()))
    return glycan_dict


if __name__ == "__main__":
    main(input_file=sys.argv[1], arcpp_pep_file=sys.argv[2])
