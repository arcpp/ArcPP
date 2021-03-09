#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import os
import csv
import pickle
import copy
import pprint
import plotly
import plotly.graph_objs as go
import operator
import statistics
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests
from itertools import combinations
from operator import itemgetter
import copy
import numpy as np

db_fasta = 'Haloferax_volcanii_FPF_ArcPP_20190606_unipro_cRAP_target_decoy_trypsin.fasta'

PRIDE_order = [
    'PXD007061',
    'PXD013046',
    'PXD011012',
    'PXD006877',
    'PXD011218',
    'PXD009116',
    'PXD011056',
    'PXD000202',
    'PXD011015',
    'PXD011050',
    'PXD014974',
    'PXD010824',
    'all',
]
PRIDE_whole_cells = [
    'PXD007061',
    'PXD013046',
    'PXD011012',
    'PXD006877',
    'PXD011218',
    'PXD009116',
    'PXD011056',
    # 'PXD000202',
    # 'PXD011015',
    # 'PXD011050',
    # 'PXD014974',
    # 'PXD010824',
    # 'all',
]
PRIDE_COLORS = {
    'PXD007061':'rgb(166,206,227)',
    'PXD013046':'rgb(31,120,180)',
    'PXD011012':'rgb(178,223,138)',
    'PXD006877':'rgb(51,160,44)',
    'PXD011218':'rgb(251,154,153)',
    'PXD009116':'rgb(227,26,28)',
    'PXD011056':'rgb(253,191,111)',
    'PXD000202':'rgb(255,127,0)',
    'PXD011015':'rgb(202,178,214)',
    'PXD011050':'rgb(106,61,154)',
    'PXD014974':'rgb(255,255,153)',
    'PXD010824':'rgb(177,89,40)',
    'all':'rgb(0,0,0)',
}

# arcog_letters = ['J','A','K','L','B','X','D','Y','V','T','M','N','Z','W','U','O','C','G','E','F','H','I','P','Q']
# arcog_letters = ['J','U','F','H','O','Q','I','V','E','C','L','T','R','D','P','M','G','K','S','N','X','A','B','Y','Z','W']
arcog_letters = ['J','U','F','H','O','Q','I','V','E','C','L','T','D','P','K','M','G','N','X','A','B','Y','Z','W']
arcog_letters_2 = ['J','U','F','H','O','Q','I','V','E','C','L','T','R','D','P','K','M','G','S','N','X','Not classified','A','B','Y','Z','W']

csv_kwargs = {}
if sys.platform == 'win32':
    csv_kwargs['lineterminator'] = '\n'
else:
    csv_kwargs['lineterminator'] = '\r\n'

    
def main(results_pkl=None, fdr_pkl=None, protein_pkl=None, peptide_pkl=None, plot_types=[]):
    results_dict = None
    fdr_dict = None
    prot_info_dict = None
    peptide_info_dict = None
    for plot in plot_types:
        if plot in [
            'plot_spec_number',
            'plot_comparison_org',
            'plot_comparison_detailed',
            'plot_psm_ident_rate',
            'plot_peptide_ident_rate',
            'plot_total_peptides_proteins',
            'plot_prot_coverage',
            'plot_dataset_overlap',
            'write_prot2file',
            'write_results2csv',
            'plot_prot_properties',
            'plot_peptide_distributions',
            'plot_predictions',
            'plot_fdr_results',
        ] and results_dict is None:
            with open(results_pkl, 'rb') as results_in: 
                results_dict = pickle.load(results_in)
        if plot in ['plot_fdr_results'] and fdr_dict is None:
            with open(fdr_pkl, 'rb') as fdr_in: 
                fdr_dict = pickle.load(fdr_in)
        if plot in [
            'write_results2csv',
            'plot_prot_properties',
            'plot_predictions',
            'plot_dataset_overlap',
        ] and prot_info_dict is None:
            with open(protein_pkl, 'rb') as prot_info_pkl:
                prot_info_dict = pickle.load(prot_info_pkl)
        if plot in [
            'plot_peptide_distributions',
        ] and peptide_info_dict is None:
            with open(peptide_pkl, 'rb') as pep_info_pkl:
                peptide_info_dict = pickle.load(pep_info_pkl)

    if 'plot_spec_number' in plot_types:
        print('plotting spec number')
        plot_spec_number(
            results_dict,
            color_by=None,
            output_filename='barplot_spectra.pdf'
        )

    if 'plot_comparison_org' in plot_types:
        print('plotting comaprison')
        plot_comparison_org(
            results_dict,
            output_filename='barplot_comparison.pdf',
            # safe_idents=True,
        )

    if 'plot_comparison_detailed' in plot_types:
        print('plotting detailed comaprison')
        plot_comparison_detailed(
            results_dict,
            output_filename='barplot_comparison_detail.pdf',
        )

    if 'plot_fdr_results' in plot_types:
        print('plotting FDR results')
        inference_dict = protein_inference(results_dict)
        plot_fdr_results(
            fdr_dict,
            output_filename='fdr_results_proteins.pdf',
            level = 'proteins',
            inference_dict=inference_dict,
            results_dict=results_dict,
            # safe_idents=True,
        )

    if 'plot_psm_ident_rate' in plot_types:
        print('plotting total PSMs and ident rate')
        plot_psm_ident_rate(
            results_dict,
            output_filename='total_psms_ident_rate.pdf'
        )

    if 'plot_peptide_ident_rate' in plot_types:
        print('plotting total PSMs and ident rate')
        plot_peptide_ident_rate(
            results_dict,
            output_filename='total_peptide_ident_rate.pdf'
        )

    if 'plot_total_peptides_proteins' in plot_types:
        print('plotting total peptides and proteins')
        # inference_dict = protein_inference(results_dict)
        plot_total_peptides_proteins(
            results_dict,
            # inference_dict=inference_dict,
            output_filename='total_sequences_proteins_new.pdf'
        )

    if 'plot_prot_coverage' in plot_types:
        print('plotting protein and proteome coverage')
        inference_dict = protein_inference(results_dict)
        plot_prot_coverage(
            results_dict,
            inference_dict=inference_dict,
            output_filename='prot_coverage_new.pdf'
        )

    if 'plot_dataset_overlap' in plot_types:
        print('plotting overlap between datasets')
        inference_dict = protein_inference(results_dict)
        plot_dataset_overlap(
            results_dict,
            output_filename='dataset_overlap.pdf',
            inference_dict=inference_dict,
            prot_info_dict=prot_info_dict,
        )

    if 'write_prot2file' in plot_types:
        print('writing proteins')
        write_prot2file(
            results_dict,
            # output_filename='dataset_overlap.pdf'
        )

    if 'write_results2csv' in plot_types:
        print('writing results to csv')
        inference_dict = protein_inference(results_dict)
        write_results2csv(
            results_dict,
            inference_dict=inference_dict,
            output_filename='ArcPP_results.csv',
            prot_info_dict=prot_info_dict,
        )

    if 'plot_prot_properties' in plot_types:
        print('plotting properties of identified proteins')
        inference_dict = protein_inference(results_dict)
        plot_prot_properties(
            results_dict,
            inference_dict=inference_dict,
            output_filename='ArcPP_properties.csv',
            prot_info_dict=prot_info_dict,
        )

    if 'plot_peptide_distributions' in plot_types:
        print('plotting distributions of identified peptides')
        plot_peptide_distributions(
            results_dict,
            output_filename='ArcPP_peptide_distribution.csv',
            peptide_info_dict=peptide_info_dict,
        )

    if 'plot_predictions' in plot_types:
        print('plotting graphs for predictions of identified peptides')
        inference_dict = protein_inference(results_dict)
        plot_predictions(
            results_dict,
            output_filename='Predictions_graph.pdf',
            prot_info_dict=prot_info_dict,
            inference_dict=inference_dict,
            spec_filter=2,
        )


def plot_predictions(results_dict, output_filename=None, inference_dict=None, prot_info_dict=None, spec_filter=2):
    '''
    Based on various prediction engines, proteins were sorted into different categories according
    to their predicted secretion pathway/processing.
    For these categories, different plots show the corresponding ArcPP results:
    - number/percentage of identified proteins within each category
    - number of idetified N- and C-termini within each category
    - N-terminal protein processing for proteins in each category
    - ratio of TM/Cyt proteins for each dataset
    '''

    uc = ursgal.UController(verbose=False)
    venndiagram_main = uc.unodes['venndiagram_1_1_0']['class'].import_engine_as_python_function()
    
    #get all termini from database
    prot_term_dict = {}
    with open(db_fasta, 'r') as db_input:
        total_proteome_length = 0
        for fasta_id, db_sequence in ursgal.ucore.parse_fasta(db_input):
            if 'decoy_' in fasta_id or 'spurious' in fasta_id:
                continue
            prot_id = fasta_id.split(' ')[0]
            if 'HVO' not in prot_id:
                continue
            if len(db_sequence) <= 51:
                prot_term_dict[fasta_id] = {
                    'N-term': db_sequence,
                    'C-term': db_sequence,
                }
            else:
                prot_term_dict[fasta_id] = {
                    'N-term': db_sequence[:51],
                    'C-term': db_sequence[-51:],
                }

    plot_1_dict = {
        'predicted': {},
        'identified': {},
    }

    plot_term_dict = {
        'N-term': {
            'all proteins': {},
            'identifiable': {},
            'identified': {},
        },
        'C-term': {
            'all proteins': {},
            'identifiable': {},
            'identified': {},
        },
    }

    plot_4_dict = {
        'terminal M': {},
        'acetylated terminal M': {},
        'cleaved N-term': {},
        'acetylated cleaved N-term': {},
    }

    plot_5_dict = {
        'Cyt': {},
        # 'Cytosolic': {},
        '>=2 TM' : {},
        # 'Transmembrane': {},
        # 'Secreted': {},
        'identified': {},
    }

    print('collecting info')
    prot2cov = plot_prot_coverage(results_dict, output_filename='tmp.pdf', inference_dict=inference_dict)
    for k, fasta_id in enumerate(prot_term_dict.keys()):
        if k % 100 == 0:
            print('protein:', k)
        org_fasta_id = fasta_id.split('_NumberOfIdenticalSequences')[0]
        decision_tree = prot_info_dict[org_fasta_id]['decision_tree']
        if decision_tree not in plot_1_dict['predicted'].keys():
            plot_1_dict['predicted'][decision_tree] = set()
            plot_1_dict['identified'][decision_tree] = set()
        for cat in plot_4_dict.keys():
            if decision_tree not in plot_4_dict[cat].keys():
                plot_4_dict[cat][decision_tree] = set()
        plot_1_dict['predicted'][decision_tree].add(fasta_id)
        for term in ['N-term', 'C-term']:
            if decision_tree not in plot_term_dict[term]['all proteins'].keys():
                plot_term_dict[term]['all proteins'][decision_tree] = set()
            if decision_tree not in plot_term_dict[term]['identifiable'].keys():
                plot_term_dict[term]['identifiable'][decision_tree] = set()
            if decision_tree not in plot_term_dict[term]['identified'].keys():
                plot_term_dict[term]['identified'][decision_tree] = set()
            plot_term_dict[term]['all proteins'][decision_tree].add(fasta_id)
            term_seq = prot_term_dict[fasta_id][term]
            try_count = term_seq.count('K') + term_seq.count('R')
            gluc_count = term_seq.count('D') + term_seq.count('E')
            # if try_count <= 2 or gluc_count <= 3:
            if term == 'N-term':
                if 'K' not in term_seq[:6] and 'R' not in term_seq[:6]:
                    if try_count > 0 or len(term_seq) <= 50:
                        plot_term_dict[term]['identifiable'][decision_tree].add(fasta_id)
                if 'D' not in term_seq[:6] and 'E' not in term_seq[:6]:
                    if gluc_count > 0 or len(term_seq) <= 50:
                        plot_term_dict[term]['identifiable'][decision_tree].add(fasta_id)
            if term == 'C-term':
                if 'K' not in term_seq[-6:] and 'R' not in term_seq[-6:]:
                    if try_count > 0 or len(term_seq) <= 50:
                        plot_term_dict[term]['identifiable'][decision_tree].add(fasta_id)
                if 'D' not in term_seq[-6:] and 'E' not in term_seq[-6:]:
                    if gluc_count > 0 or len(term_seq) <= 50:
                        plot_term_dict[term]['identifiable'][decision_tree].add(fasta_id)
        for prot in sorted(results_dict['all']['protein_dict'].keys()):
            if prot not in prot2cov.keys():
                continue
            if fasta_id in prot:
                plot_1_dict['identified'][decision_tree].add(prot)
                for seq in results_dict['all']['protein_dict'][prot].keys():
                    if seq in ['datasets','prot_q_value_psm', 'prot_q_value_seq', 'samples']:
                        continue
                    if seq not in results_dict['all']['peptides']['safe_num_specs']:
                        continue

                    term_M = 0
                    term_cleaved = 0
                    term_acetyl_M = 0
                    term_acetyl_cleaved = 0
                    start, stop = results_dict['all']['protein_dict'][prot][seq]['start_stop']
                    start = start.split('<|>')
                    stop = stop.split('<|>')
                    for n, p in enumerate(prot.split('<|>')):
                        org_p = p.split('_NumberOfIdenticalSequences')[0]
                        min_start = start[n].split(';')[0]
                        max_stop = stop[n].split(';')[-1]
                        if int(min_start) == 1:
                            plot_term_dict['N-term']['identified'][decision_tree].add(prot)
                            for mod in results_dict['all']['protein_dict'][prot][seq]['modifications']:
                                if 'Acetyl' in mod:
                                    term_acetyl_M += 1
                                else:
                                    term_M += 1
                                if seq[0] != 'M':
                                    print(seq)
                        elif int(min_start) == 2:# or int(min_start) == 3:
                            plot_term_dict['N-term']['identified'][decision_tree].add(prot)
                            for mod in results_dict['all']['protein_dict'][prot][seq]['modifications']:
                                if 'Acetyl' in mod:
                                    term_acetyl_cleaved += 1
                                else:
                                    term_cleaved += 1
                        if int(max_stop) == len(prot_info_dict[org_p]['sequence']):
                            plot_term_dict['C-term']['identified'][decision_tree].add(prot)
                    #PSM filter
                    if term_M >= spec_filter:
                        plot_4_dict['terminal M'][decision_tree].add(prot)
                    if term_cleaved >= spec_filter:
                        plot_4_dict['cleaved N-term'][decision_tree].add(prot)
                    if term_acetyl_M >= spec_filter:
                        plot_4_dict['acetylated terminal M'][decision_tree].add(prot)
                    if term_acetyl_cleaved >= spec_filter:
                        plot_4_dict['acetylated cleaved N-term'][decision_tree].add(prot)

        for PRIDE_ID in PRIDE_order:
            if PRIDE_ID == 'all':
                name = 'All datasets'
            else:
                name = PRIDE_ID
                # continue
            if k == 0:
                for cat in plot_5_dict.keys():
                    if name not in plot_5_dict[cat].keys():
                        plot_5_dict[cat][name] = set()

            proteins = results_dict[PRIDE_ID]['proteins']['safe_seq_num_spec_0005']
            protein_groups = results_dict[PRIDE_ID]['protein_groups']['safe_seq_num_spec_0005']

            for group in protein_groups:
                group_tmp = set()
                inference = inference_dict.get(group, None)
                if inference == set():
                    for prot in group.split('<|>'):
                        org_prot = prot.split('_NumberOfIdenticalSequences_')[0]
                        for cat in plot_5_dict.keys():
                            if prot_info_dict[org_prot]['decision_tree'] == cat:
                                group_tmp.add(cat)
                if len(group_tmp) == 1:
                    plot_5_dict[group_tmp.sorted()[0]][name].add(group)
                    plot_5_dict['identified'][name].add(group)
                else:
                    continue

            for protein in proteins:
                org_prot = protein.split('_NumberOfIdenticalSequences_')[0]
                for cat in plot_5_dict.keys():
                    if prot_info_dict[org_prot]['decision_tree'] == cat:
                        plot_5_dict[cat][name].add(protein)
                        plot_5_dict['identified'][name].add(protein)

    plot_2_dict = plot_term_dict['N-term']
    plot_3_dict = plot_term_dict['C-term']
    print('preparing to plot')

    color_dict = {
        'predicted': '#d9f0a3',
        'identified': '#2ca02c',
        'all proteins': '#d9f0a3',
        'identifiable': '#78c679',
        'terminal M': '#542788',
        'acetylated terminal M': '#998ec3',
        'cleaved N-term': '#b35806',
        'acetylated cleaved N-term': '#f1a340',
        'Cyt': '#a6bddb',
        '>=2 TM': '#02818a',
        # 'Cytosolic': '#a6bddb',
        # 'Transmembrane': '#02818a',
    }

    name_dict = {
        0: 'proteins',
        1: 'N_terms',
        2: 'C_terms',
        3: 'N_term_processing',
        4: 'Cyt-TM'
    }

    cat_sort = [
        'all proteins',
        'predicted',
        'identifiable',
        'identified',
        'terminal M',
        'acetylated terminal M',
        'cleaved N-term',
        'acetylated cleaved N-term',
        'Cyt',
        # 'Cytosolic',
        '>=2 TM',
        # 'Transmembrane',
    ]
    pred_sort = [
        'Cyt',
        '>=2 TM',
        # 'Cytosolic',
        # 'Transmembrane',
        # 'Secreted',
        # 'Cyt pI>=7',
        # 'Cyt pI<7',
        # '>1 TM pI>7',
        # '>1 TM pI<7',
        '1 TM',
        'TM N-term',
        'TM C-term',
        'Sec; SPI',
        'Tat; SPI',
        'Sec; lipobox',
        'Tat; lipobox',
        'Pil; Sec SPIII',
        'PXD007061',
        'PXD013046',
        'PXD011012',
        'PXD006877',
        'PXD011218',
        'PXD009116',
        'PXD011056',
        'PXD000202',
        'PXD011015',
        'PXD011050',
        'PXD014974',
        'PXD010824',
        'All datasets',
    ]

    for i, plot in enumerate([
        plot_1_dict,
        plot_2_dict,
        plot_3_dict, 
        plot_4_dict,
        plot_5_dict,
    ]):
        trace_list = []
        for cat in cat_sort:
            if cat not in plot.keys():
                continue
            if i == 4 and cat == 'identified':
                continue
            x_bar = []
            y_bar = []
            for decision_tree in pred_sort:
                if decision_tree not in plot[cat].keys():
                    continue
                x_bar.append(decision_tree)
                count = len(plot[cat][decision_tree])
                y_bar.append(count)

            if i == 4:
                fieldnames = [
                    'PRIDE ID',
                    '#Proteins (Cyt)',
                    '#Proteins (>=2 TM)',
                    'Ratio >=2 TM/Cyt (in %)',
                ]
                csv_out_name = output_filename.replace('.pdf', '_{0}_{1}.csv'.format(name_dict[i], cat.strip('>=')))
                with open(csv_out_name, 'w') as csv_out:
                    csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
                    csv_writer.writeheader()
                    for n, x in enumerate(x_bar):
                        out_dict = {}
                        out_dict['PRIDE ID'] = x
                        out_dict['#Proteins ({0})'.format(cat)] = y_bar[n]
                        csv_writer.writerow(out_dict)

            # if i == 3:
            #     fieldnames = [
            #         'Category',
            #         '#Proteins (terminal M)',
            #         '#Proteins (acetylated terminal M)',
            #         '#Proteins (cleaved N-term)',
            #         '#Proteins (acetylated cleaved N-term)',
            #     ]
            #     csv_out_name = output_filename.replace('.pdf', '_{0}_{1}.csv'.format(name_dict[i], cat))
            #     with open(csv_out_name, 'w') as csv_out:
            #         csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
            #         csv_writer.writeheader()
            #         for n, x in enumerate(x_bar):
            #             out_dict = {}
            #             out_dict['Category'] = x
            #             out_dict['#Proteins ({0})'.format(cat)] = y_bar[n]
            #             csv_writer.writerow(out_dict)

            trace_bar = go.Bar(
                x = x_bar,
                y = y_bar,
                name = cat,
                marker=dict(
                    color=color_dict[cat],
                )
             )
            trace_list.append(trace_bar)

            if i == 0 and cat == 'identified':
                x_scatter = []
                y_scatter = []
                fieldnames = [
                    'Category',
                    '#Proteins predicted',
                    '#Proteins identified',
                    'Ratio identified/predicted (in %)',
                ]
                csv_out_name = output_filename.replace('.pdf', '_proteins_detailed.csv')
                with open(csv_out_name, 'w') as csv_out:
                    csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
                    csv_writer.writeheader()
                    for decision_tree in pred_sort:
                        if decision_tree not in plot[cat].keys():
                            continue
                        out_dict = {}
                        out_dict['Category'] = decision_tree
                        x_scatter.append(decision_tree)
                        total_count = len(plot['predicted'][decision_tree])
                        specific_count = len(plot[cat][decision_tree])
                        y_scatter.append(100*specific_count/total_count)
                        out_dict['#Proteins predicted'] = total_count
                        out_dict['#Proteins identified'] = specific_count
                        out_dict['Ratio identified/predicted (in %)'] = 100*specific_count/total_count
                        csv_writer.writerow(out_dict)

                trace_scatter = go.Scatter(
                    x = x_scatter,
                    y = y_scatter,
                    name = cat,
                    mode='markers',
                    marker=dict(
                        color = '#4d4d4d',
                        symbol='cross',
                        size=5,
                    ),
                    yaxis='y2'
                )
                trace_list.append(trace_scatter)

            elif cat in ['all proteins', 'identifiable']:
                color_dict_scatter = {
                    'all proteins': '#4d4d4d',
                    'identifiable': '#999999',
                }
                x_scatter = []
                y_scatter = []
                fieldnames = [
                    'Category',
                    '#terms predicted',
                    '#terms identifiable',
                    '#terms identified',
                    'Ratio identified/predicted (in %)',
                    'Ratio identified/identifiable (in %)',
                ]
                csv_out_name = output_filename.replace('.pdf', '_{0}_{1}.csv'.format(name_dict[i], cat))
                with open(csv_out_name, 'w') as csv_out:
                    csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
                    csv_writer.writeheader()
                    for decision_tree in pred_sort:
                        if decision_tree not in plot[cat].keys():
                            continue
                        out_dict = {}
                        out_dict['Category'] = decision_tree
                        x_scatter.append(decision_tree)
                        total_count = len(plot[cat][decision_tree])
                        specific_count = len(plot['identified'][decision_tree])
                        y_scatter.append(100*specific_count/total_count)
                        out_dict['#terms predicted'] = total_count
                        out_dict['#terms identified'] = specific_count
                        out_dict['#terms identifiable'] = total_count
                        out_dict['Ratio identified/predicted (in %)'] = 100*specific_count/total_count
                        out_dict['Ratio identified/identifiable (in %)'] = 100*specific_count/total_count
                        csv_writer.writerow(out_dict)

                trace_scatter = go.Scatter(
                    x = x_scatter,
                    y = y_scatter,
                    name = cat,
                    mode='markers',
                    marker=dict(
                        color = color_dict_scatter[cat],
                        symbol='cross',
                        size=5,
                    ),
                    yaxis='y2'
                )
                trace_list.append(trace_scatter)

            elif cat in [
                'terminal M',
                'acetylated terminal M',
                'cleaved N-term',
                'acetylated cleaved N-term',
                'Cyt',
                '>=2 TM',
                # 'Cytosolic',
                # 'Transmembrane',
            ]:
                color_dict_scatter = {
                    'terminal M': '#252525',
                    'acetylated terminal M': '#969696',
                    'cleaved N-term': '#636363',
                    'acetylated cleaved N-term': '#bdbdbd',
                    'Cyt': '#252525',
                    '>=2 TM':'#252525',
                    # 'Cytosolic': '#252525',
                    # 'Transmembrane':'#252525',
                }
                x_scatter = []
                y_scatter = []
                for decision_tree in pred_sort:
                    if decision_tree not in plot[cat].keys():
                        continue
                    x_scatter.append(decision_tree)
                    if cat in ['Cyt', '>=2 TM']:
                    # if cat in ['Cytosolic', 'Transmembrane']:
                        total_count = len(plot_5_dict['identified'][decision_tree])
                    else:
                        total_count = len(plot_2_dict['identified'][decision_tree])
                    specific_count = len(plot[cat][decision_tree])
                    # if decision_tree not in [
                    #     'Cyt',
                    #     '>=2 TM',
                    #     '1 TM',
                    #     'TM N-term',
                    #     'TM C-term',
                    # ] and i==3:
                        # print(cat, decision_tree, plot[cat][decision_tree])
                    if total_count == 0:
                        y_scatter.append(0)
                    else:
                        if cat == 'terminal M':
                            specific_count = len(
                                plot['acetylated cleaved N-term'][decision_tree] | plot['cleaved N-term'][decision_tree]
                            )
                        elif cat == 'acetylated terminal M':
                            specific_count = len(
                                plot['acetylated cleaved N-term'][decision_tree] | plot['acetylated terminal M'][decision_tree]
                            )
                        elif cat in ['cleaved N-term', 'acetylated cleaved N-term']:
                            del x_scatter[-1]
                            continue
                        y_scatter.append(100*specific_count/total_count)
                        print('total', decision_tree, ':', total_count)
                        print('specific', decision_tree, ':', specific_count)
                        print('ratio', decision_tree, ':', 100*specific_count/total_count)

                trace_scatter = go.Scatter(
                    x = x_scatter,
                    y = y_scatter,
                    name = cat,
                    mode='markers',
                    marker=dict(
                        color = color_dict_scatter[cat],
                        symbol='cross',
                        size=5,
                    ),
                    yaxis='y2'
                )
                trace_list.append(trace_scatter)


        layout = go.Layout(
            barmode='group',
            yaxis=dict(
                title='Proteins',
                # rangemode='tozero',
                showline=True,
                linecolor='rgb(0, 0, 0)',
                titlefont=dict(
                    size=16,
                    color='rgb(0, 0, 0)'
                ),
                tickfont=dict(
                    size=14,
                    color='rgb(0, 0, 0)',
                ),
                ticks='outside',
                ticklen=2,
                # dtick=500,
                dtick=100,
                tickwidth=0.25,
                tickcolor='rgb(0, 0, 0)',
                # tickangle=-90,
                side='left',
                autorange=True,
                # range=[0,1.5]
            ),
            yaxis2=dict(
                # title='Identification rate in %',
                title='Percent of total',
                # rangemode='tozero',
                showgrid=False,
                showline=True,
                linecolor='rgb(0, 0, 0)',
                titlefont=dict(
                    size=16,
                    color='rgb(0, 0, 0)'
                ),
                dtick=10,
                tickfont=dict(
                    size=14,
                    color='rgb(0, 0, 0)'
                ),
                ticks='outside',
                ticklen=2,
                tickwidth=0.25,
                tickcolor='rgb(0, 0, 0)',
                overlaying='y',
                side='right',
                autorange=False,
                range=[0,102]
            ),
            xaxis=dict(
                showgrid=False,
                showline=False,
                zeroline=True,
                ticks='outside',
                ticklen=2,
                tickwidth=0.25,
                title='',
                titlefont=dict(
                    size=16,
                    color='rgb(0, 0, 0)'
                ),
                tickangle=-90,
                tickfont=dict(
                    size=14,
                    color='rgb(0, 0, 0)'
                )
            ),
            showlegend=False,
            legend=dict(
                orientation= "v",
                x= 1.3,
                xanchor= 'left',
                y= 1,
                # font: {
                #   size: 9
                # },
                #     borderwidth: 1
            ),
            autosize=False,
            # width=500,
            # height=400,
            width=300,
            height=350,
            margin=go.layout.Margin(
                l=65,
                r=50,
                # b=200,
                b=150,
                t=10,
                # pad=4
            ),
        )

        final_output_filename = output_filename.replace('.pdf', '_{0}.pdf'.format(name_dict[i]))
        fig = go.Figure(data=trace_list, layout=layout)
        plotly.io.write_image(fig, final_output_filename)
        plotly.io.write_html(fig, final_output_filename.replace('.pdf', '.html'))
        print('plotted: ', final_output_filename)


    for to_plot in ['Cyt', '>=2 TM']:
    # for to_plot in ['Cytosolic', 'Transmembrane']:
        venn_data = []
        for cat in [
            'terminal M',
            'acetylated terminal M',
            'cleaved N-term',
            'acetylated cleaved N-term',
        ]:
            venn_dict = {}
            venn_dict['label'] = cat
            venn_dict['data'] = plot_4_dict[cat][to_plot]
            venn_data.append(venn_dict)
        if to_plot == 'Cyt':
        # if to_plot == 'Cytosolic':
            name = to_plot
        else:
            name = 'TM'
        final_output_filename = output_filename.replace('.pdf', '_venn_{0}_{1}PSMs.svg'.format(name, spec_filter))
        venndiagram_main(
            data=venn_data,
            output_file = final_output_filename
        )

def write_results2csv(results_dict, output_filename=None, inference_dict=None, prot_info_dict=None):
    '''
    Results are written into three separate csv files:
    One contains protein identifications, one peptide identifications and one all PSMs.
    Only safe identifications are saved in the files.
    In addition, results for specific proteins of interest can be written into a file.
    '''
    #protein centric
    fieldnames = [
        'HVO ID',
        'Uniprot ID',
        'Description',
        'Protein q-value',
        'Peptides',
        'PSMs',
        'Sequence Coverage',
        'Dataset',
        'Total Datasets',
        'Predicted Processing',
        'Molecular Weight',
        'pI',
        'Hydrophobicity',
        'arCOG',
        'arCOGlet',
    ]

    prot_db_dict = {}
    with open(db_fasta, 'r') as db_input:
        for fasta_id, db_sequence in ursgal.ucore.parse_fasta(db_input):
            if 'decoy_' in fasta_id or 'spurious' in fasta_id:
                continue
            prot_id = fasta_id.split(' ')[0]
            if 'HVO' not in prot_id:
                continue
            prot_db_dict[fasta_id] = False

    prot2cov = plot_prot_coverage(results_dict, output_filename='tmp.pdf', inference_dict=inference_dict, all_datasets=True)

    count_nonidentified = 0
    output_filename_proteins = output_filename.replace('.csv', '_proteins.csv')
    with open(output_filename_proteins, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for PRIDE_ID in PRIDE_order:
            out_dict = {}    
            if PRIDE_ID == 'all':
                out_dict['Dataset'] = 'Total'
            else:
                out_dict['Dataset'] = PRIDE_ID
            for prot in sorted(results_dict[PRIDE_ID]['protein_dict'].keys()):
                if prot not in results_dict[PRIDE_ID]['proteins']['safe_seq_num_spec_0005'] and \
                    prot not in results_dict[PRIDE_ID]['protein_groups']['safe_seq_num_spec_0005']:
                    continue
                datasets = 0
                for pid in PRIDE_order:
                    # if pid == 'all':
                    #     continue
                    if prot in results_dict[pid]['proteins']['safe_seq_num_spec_0005'] or \
                        prot in results_dict[pid]['protein_groups']['safe_seq_num_spec_0005']:
                        datasets += 1
                prot_id = prot.split(' ')
                if prot not in prot2cov[PRIDE_ID].keys():
                    continue
                # out_dict['Protein ID'] = prot
                out_dict['HVO ID'] = prot_id[0]
                out_dict['Uniprot ID'] = prot_id[1].strip('[').strip(']')
                out_dict['Description'] = ' '.join(prot_id[2:])
                out_dict['Protein q-value'] = results_dict[PRIDE_ID]['protein_dict'][prot]['prot_q_value_seq']
                # datasets = sorted(results_dict['all']['protein_dict'][prot]['datasets'])
                # out_dict['Datasets'] = ';'.join(datasets)
                out_dict['Sequence Coverage'] = prot2cov[PRIDE_ID][prot][0]
                out_dict['Total Datasets'] = datasets
                out_dict['Peptides'] = 0
                out_dict['PSMs'] = 0
                for seq in results_dict[PRIDE_ID]['protein_dict'][prot].keys():
                    if seq in ['datasets','prot_q_value_psm', 'prot_q_value_seq', 'samples']:
                        continue
                    elif seq in results_dict[PRIDE_ID]['peptides']['safe_num_specs']:
                        out_dict['Peptides'] += 1
                        out_dict['PSMs'] += len(results_dict[PRIDE_ID]['protein_dict'][prot][seq]['spec_title'])
                if out_dict['Peptides'] == 0:
                    continue
                org_prot = prot.split('_NumberOfIdenticalSequences_')[0]
                if '<|>' in org_prot:
                    prediction_list = []
                    mw_list = []
                    pi_list = []
                    hyd_list = []
                    arcog_list = []
                    arcoglet_list = []
                    for p in org_prot.split('<|>'):
                        prediction_list.append(prot_info_dict[p]['decision_tree'])
                        mw_list.append(str(prot_info_dict[p]['molecular_weight']))
                        pi_list.append(str(prot_info_dict[p]['pi']))
                        hyd_list.append(str(prot_info_dict[p]['hydrophobicity']))
                        arcog_list.append(str(prot_info_dict[p]['arcog']))
                        arcoglet_list.append(str(prot_info_dict[p]['arcoglet']))
                    out_dict['Predicted Processing'] = '<|>'.join(prediction_list)
                    out_dict['Molecular Weight'] = '<|>'.join(mw_list)
                    out_dict['pI'] = '<|>'.join(pi_list)
                    out_dict['Hydrophobicity'] = '<|>'.join(hyd_list)
                    out_dict['arCOG'] = '<|>'.join(arcog_list)
                    out_dict['arCOGlet'] = '<|>'.join(arcoglet_list)
                else:
                    out_dict['Predicted Processing'] = prot_info_dict[org_prot]['decision_tree']
                    out_dict['Molecular Weight'] = prot_info_dict[org_prot]['molecular_weight']
                    out_dict['pI'] = prot_info_dict[org_prot]['pi']
                    out_dict['Hydrophobicity'] = prot_info_dict[org_prot]['hydrophobicity']
                    out_dict['arCOG'] = prot_info_dict[org_prot]['arcog']
                    out_dict['arCOGlet'] = prot_info_dict[org_prot]['arcoglet']
                if out_dict['arCOG'] == '':
                    out_dict['arCOG'] = 'nd'
                    out_dict['arCOGlet'] = 'nd'
                csv_writer.writerow(out_dict)
                if '<|>' not in prot:
                    prot_db_dict[prot] = True
        for prot in prot_db_dict.keys():
            if prot_db_dict[prot] is True:
                continue
            if 'spurious' in prot:
                continue
            count_nonidentified+=1
            out_dict = {} 
            out_dict['Dataset'] = 'None'
            out_dict['Total Datasets'] = 0
            prot_id = prot.split(' ')
            out_dict['HVO ID'] = prot_id[0]
            out_dict['Uniprot ID'] = prot_id[1].strip('[').strip(']')
            out_dict['Description'] = ' '.join(prot_id[2:])
            out_dict['Protein q-value'] = 1
            out_dict['Sequence Coverage'] = 0
            out_dict['Peptides'] = 0
            out_dict['PSMs'] = 0
            org_prot = prot.split('_NumberOfIdenticalSequences_')[0]
            out_dict['Predicted Processing'] = prot_info_dict[org_prot]['decision_tree']
            out_dict['Molecular Weight'] = prot_info_dict[org_prot]['molecular_weight']
            out_dict['pI'] = prot_info_dict[org_prot]['pi']
            out_dict['Hydrophobicity'] = prot_info_dict[org_prot]['hydrophobicity']
            out_dict['arCOG'] = prot_info_dict[org_prot]['arcog']
            out_dict['arCOGlet'] = prot_info_dict[org_prot]['arcoglet']
            if out_dict['arCOG'] == '':
                out_dict['arCOG'] = 'nd'
                out_dict['arCOGlet'] = 'nd'
            csv_writer.writerow(out_dict)

    #peptide centric
    fieldnames = [
        'Sequence',
        'HVO ID',
        'Uniprot ID',
        'Description',
        'Sequence Start',
        'Sequence Stop',
        'Sequence q-value',
        'PSMs',
        'Dataset',
        'Files'
    ]

    output_filename_peptides = output_filename.replace('.csv', '_peptides.csv')
    with open(output_filename_peptides, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for PRIDE_ID in PRIDE_order:
            out_dict = {}    
            if PRIDE_ID == 'all':
                # out_dict['Dataset'] = 'Total'
                continue
            else:
                out_dict['Dataset'] = PRIDE_ID
            for prot in sorted(results_dict[PRIDE_ID]['protein_dict'].keys()):
                prot_id = prot.split(' ')
                if prot not in prot2cov[PRIDE_ID].keys():
                    continue
                out_dict['HVO ID'] = prot_id[0]
                out_dict['Uniprot ID'] = prot_id[1].strip('[').strip(']')
                out_dict['Description'] = ' '.join(prot_id[2:])
                for seq in results_dict[PRIDE_ID]['protein_dict'][prot].keys():
                    if seq in ['datasets','prot_q_value_psm', 'prot_q_value_seq', 'samples']:
                        continue
                    if seq not in results_dict[PRIDE_ID]['peptides']['safe_num_specs']:
                        continue
                    out_dict['Sequence'] = seq
                    out_dict['Sequence q-value'] = results_dict[PRIDE_ID]['protein_dict'][prot][seq]['seq_q_value']
                    out_dict['PSMs'] = len(results_dict[PRIDE_ID]['protein_dict'][prot][seq]['spec_title'])
                    files = set()
                    for spec_title in results_dict[PRIDE_ID]['protein_dict'][prot][seq]['spec_title']:
                        files.add(
                            '.'.join(spec_title.split('.')[:-3])
                        )
                    out_dict['Files'] = ';'.join(sorted(files))
                    start, stop = results_dict[PRIDE_ID]['protein_dict'][prot][seq]['start_stop']
                    out_dict['Sequence Start'] = start
                    out_dict['Sequence Stop'] = stop
                    csv_writer.writerow(out_dict)

    #PSM centric
    fieldnames = [
        'Sequence',
        'Modifications',
        'Spectrum Title',
        'Dataset',
        'HVO ID',
        'Uniprot ID',
        'Description',
        'Sequence Start',
        'Sequence Stop',
        'PSM q-value',
        'Charge',
    ]

    output_filename_PSMs = output_filename.replace('.csv', '_PSMs.csv')
    with open(output_filename_PSMs, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for PRIDE_ID in PRIDE_order:
            out_dict = {}    
            if PRIDE_ID == 'all':
                continue
            else:
                out_dict['Dataset'] = PRIDE_ID
            for prot in sorted(results_dict[PRIDE_ID]['protein_dict'].keys()):
                prot_id = prot.split(' ')
                if prot not in prot2cov[PRIDE_ID].keys():
                    continue
                out_dict['HVO ID'] = prot_id[0]
                out_dict['Uniprot ID'] = prot_id[1].strip('[').strip(']')
                out_dict['Description'] = ' '.join(prot_id[2:])
                for seq in results_dict[PRIDE_ID]['protein_dict'][prot].keys():
                    if seq in ['datasets','prot_q_value_psm', 'prot_q_value_seq', 'samples']:
                        continue
                    out_dict['Sequence'] = seq
                    start, stop = results_dict[PRIDE_ID]['protein_dict'][prot][seq]['start_stop']
                    out_dict['Sequence Start'] = start
                    out_dict['Sequence Stop'] = stop
                    for n, spec_title in enumerate(results_dict[PRIDE_ID]['protein_dict'][prot][seq]['spec_title']):
                        out_dict['Spectrum Title'] = spec_title
                        out_dict['Modifications'] = results_dict[PRIDE_ID]['protein_dict'][prot][seq]['modifications'][n]
                        out_dict['PSM q-value'] = results_dict[PRIDE_ID]['protein_dict'][prot][seq]['psm_q_value'][n]
                        out_dict['Charge'] = results_dict[PRIDE_ID]['protein_dict'][prot][seq]['charge'][n]
                        csv_writer.writerow(out_dict)

    #single protein
    fieldnames = [
        'Spectrum Title',
        'Sequence',
        'Modifications',
        'Protein ID',
        'Sequence Start',
        'Sequence Stop',
        'PSM q-value',
        'Charge',
    ]

    for poi in [
        # 'HVO_1091',
        # 'HVO_2923',
        # 'HVO_1562',
        # 'HVO_1957',
        # 'HVO_0850',
        # 'HVO_1524 '
    ]:
        output_filename_PSMs = output_filename.replace('.csv', '_PSMs_{0}.csv'.format(poi))
        with open(output_filename_PSMs, 'w') as csv_out:
            csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
            csv_writer.writeheader()
            for prot in sorted(results_dict['all']['protein_dict'].keys()):
                if poi not in prot:
                    continue
                for seq in results_dict['all']['protein_dict'][prot].keys():
                    out_dict = {}
                    out_dict['Protein ID'] = prot
                    if seq in ['datasets','prot_q_value_psm', 'prot_q_value_seq', 'samples']:
                        continue
                    out_dict['Sequence'] = seq
                    start, stop = results_dict['all']['protein_dict'][prot][seq]['start_stop']
                    out_dict['Sequence Start'] = start
                    out_dict['Sequence Stop'] = stop
                    for n, spec_title in enumerate(results_dict['all']['protein_dict'][prot][seq]['spec_title']):
                        out_dict['Spectrum Title'] = spec_title
                        out_dict['Modifications'] = results_dict['all']['protein_dict'][prot][seq]['modifications'][n]
                        out_dict['PSM q-value'] = results_dict['all']['protein_dict'][prot][seq]['psm_q_value'][n]
                        out_dict['Charge'] = results_dict['all']['protein_dict'][prot][seq]['charge'][n]
                        csv_writer.writerow(out_dict)

def plot_dataset_overlap(results_dict, output_filename=None, inference_dict=None, prot_info_dict=None):
    '''
    Analyzes the overlap between datasets and plots various graphs
    '''
    color_dict={
        'overlap' : ['#a6a6a6','#e08214', '#fdb863', '#d9ef8b', '#a6d96a', '#66bd63', '#1a9850', '#1b7837', '#000000'],
        'overlap_missing' : ['#fc8d59', '#d7301f', '#990000'],
        1: '#e08214',
        7: '#1b7837',
        0: '#a6a6a6',
        'Total proteome': '#000000',
    }
    marker_dict={
        'overlap' : ['circle','diamond', 'triangle-up', 'triangle-down', 'cross', 'octagon', 'pentagon', 'star', 'square'],
    }
    datadict = {
        'x': [],
        'y': [],
        'text': [],
        'name': 'overlap',
    }
    datadict_single = {
        'x': [],
        'y': [],
        'text': [],
        'name': 'overlap',
    }
    prot2cov = plot_prot_coverage(results_dict, output_filename='tmp.pdf', inference_dict=inference_dict)
    intersection_dict = {}
    dataset_lookup = {}
    all_whole_cell_proteins = set()
    prot_db_set = set()
    total_arcog_dict = {}
    with open(db_fasta, 'r') as db_input:
        total_proteome_length = 0
        for fasta_id, db_sequence in ursgal.ucore.parse_fasta(db_input):
            if 'decoy_' in fasta_id or 'spurious' in fasta_id:
                continue
            prot_id = fasta_id.split(' ')[0]
            if 'HVO' not in prot_id:
                continue
            prot_db_set.add(fasta_id)
            funct_set = set()
            org_prot = fasta_id.split('_NumberOfIdenticalSequences_')[0]
            for l in prot_info_dict[org_prot]['arcoglet'].split(';'):
                funct_set.add(l)
            for l in arcog_letters + ['R', 'S', 'total']:
                if l not in total_arcog_dict.keys():
                    total_arcog_dict[l] = 0
                if l in funct_set:
                    total_arcog_dict[l] += 1
                if l == 'total':
                    total_arcog_dict[l] += 1
            if len(funct_set) == 1 and list(funct_set)[0] == '':
                if 'Not classified' not in total_arcog_dict.keys():
                    total_arcog_dict['Not classified'] = 0
                total_arcog_dict['Not classified'] += 1

    print('Total:', total_arcog_dict['total'])

    datadict_list = []
    for n in range(1, len(PRIDE_whole_cells)+1):
        datadict['x'].append(n)
        n_set = set()
        dataset_lookup[n] = {}
        for l, combo in enumerate(combinations(PRIDE_whole_cells, n)):
            c_set = copy.deepcopy(results_dict[combo[0]]['proteins']['safe_seq_num_spec_0005'])
            for group in copy.deepcopy(results_dict[combo[0]]['protein_groups']['safe_seq_num_spec_0005']):
                if group not in prot2cov.keys():
                    continue
                else:
                    c_set.add(group)
            for c in combo:
                comp_set = copy.deepcopy(results_dict[c]['proteins']['safe_seq_num_spec_0005'])
                for group in copy.deepcopy(results_dict[c]['protein_groups']['safe_seq_num_spec_0005']):
                    if group not in prot2cov.keys():
                        continue
                    else:
                        comp_set.add(group)
                c_set &= comp_set
            dataset_lookup[n][combo] = c_set
            all_whole_cell_proteins |= c_set

            not_c_set = set()
            for PRIDE_ID in PRIDE_whole_cells:
                if PRIDE_ID in combo:
                    continue
                not_c_set |= copy.deepcopy(results_dict[PRIDE_ID]['proteins']['safe_seq_num_spec_0005'])
                for group in copy.deepcopy(results_dict[PRIDE_ID]['protein_groups']['safe_seq_num_spec_0005']):
                    if group not in prot2cov.keys():
                        continue
                    else:
                        not_c_set.add(group)

            n_set |= (c_set - not_c_set)
            if n == 1:
                datadict_single['x'].append(combo[0])
                datadict_single['y'].append(len(c_set - not_c_set))

        datadict['y'].append(len(n_set))
        intersection_dict[n] = n_set
        print(n, len(n_set))
        if n ==1 or n == 7:
            pie_colors = {
                1 : ['#f5f5f5', '#e08214'],
                7 : ['#f5f5f5', '#1b7837'],
            }
            unknown_funct_prots = set()
            arcog_dict = {}
            for l in arcog_letters + ['R', 'S', 'Not classified', 'total']:
                arcog_dict[l] = 0
            for protein_group in n_set:
                funct_set = set()
                for prot in protein_group.split('<|>'):
                    org_prot = prot.split('_NumberOfIdenticalSequences_')[0]
                    for l in prot_info_dict[org_prot]['arcoglet'].split(';'):
                        funct_set.add(l)
                for l in arcog_letters + ['R', 'S', 'total']:
                    if l in funct_set:
                        arcog_dict[l] += 1
                    if l == 'total':
                        arcog_dict[l] += 1
                if len(funct_set) == 1 and list(funct_set)[0] == '':
                    arcog_dict['Not classified'] += 1

            values = []
            labels = []
            for l in ['Not classified'] + arcog_letters + ['R', 'S']:
                if l in ['A', 'B', 'Y', 'Z', 'W']:
                    continue
                values.append(arcog_dict[l])
                labels.append(l)
            fig = go.Figure(data=[go.Pie(
                labels=labels,
                values=values,
                sort=False,
                showlegend=False,
            )])
            fig.update_traces(
                hoverinfo='label+percent+value',
                textinfo='none',
                textfont_size=20,
                marker=dict(
                    colors=['#ffffff', '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000'],
                    line=dict(
                        color='#000000',
                        width=0.5)
                ),
            )
            final_output_filename = 'dataset_overlap_pie_arcog_{0}.pdf'.format(n)
            plotly.io.write_image(fig, final_output_filename)
            plotly.io.write_html(fig, final_output_filename.replace('.pdf', '.html'))
            print('plotted: ', final_output_filename)

            with open('{0}_intersection_unknown_function_proteins.csv'.format(n), 'w') as txt_out:
                for prot in unknown_funct_prots:
                    print(prot, file=txt_out)

            return_list = fisher_test_with_correction(
                arcog_dict=arcog_dict,
                total_arcog_dict=total_arcog_dict,
                n=n,
                return_total=True,
            )
            if n == 1:
                datadict_list.append(return_list[1])
                datadict_list.append(return_list[0])
            else:
                datadict_list.append(return_list[0])

    fieldnames = [
        '#PSMs',
        'Protein ID',
        'Protein Description',
        'Datasets'
    ]
    for i in intersection_dict.keys():
        with open('{0}_intersection_proteins.csv'.format(i), 'w') as csv_out:
            csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
            csv_writer.writeheader()
            for prot in intersection_dict[i]:
                out_dict = {}
                specs = set()
                for seq in results_dict['all']['protein_dict'][prot].keys():
                    if seq in ['datasets','prot_q_value_psm', 'prot_q_value_seq', 'samples']:
                        continue
                    else:
                        specs |= set(results_dict['all']['protein_dict'][prot][seq]['spec_title'])
                dataset_combos = []
                for combo in dataset_lookup[i].keys():
                    if prot in dataset_lookup[i][combo]:
                        dataset_combos.append('&'.join(combo))
                spec_count = len(specs)
                out_dict['#PSMs'] = spec_count
                prot_id = []
                prot_descr = []
                if '<|>' in prot:
                    split_prot = prot.split('<|>')
                    for p in split_prot:
                        prot_id.append(' '.join(p.split(' ')[:2]))
                        prot_descr.append(' '.join(p.split(' ')[2:]))
                else:
                    prot_id.append(' '.join(prot.split(' ')[:2]))
                    prot_descr.append(' '.join(prot.split(' ')[2:]))
                out_dict['Protein ID'] = '<|>'.join(prot_id)
                out_dict['Protein Description'] = '<|>'.join(prot_descr)
                out_dict['Datasets'] = ';'.join(dataset_combos)
                csv_writer.writerow(out_dict)


    common_identified_prots = set()
    for identified_prot in intersection_dict[len(PRIDE_whole_cells)]:
        for p in identified_prot.split('<|>'):
            prot_id = p.split(' ')[0]
            common_identified_prots.add(prot_id)

    essential_arcog = set()
    essential_cog = set()
    with open('essential_sulfolobus_homologs.txt', 'r') as sulfolobus_in:
        for line in sulfolobus_in:
            arcog, cog = line.split('; ')
            arcog = arcog.split(': ')[1]
            for a in arcog.split(', '):
                if 'HVO_' not in a:
                    continue
                essential_arcog.add(a)
            cog = cog.split(': ')[1]
            for c in cog.split(', '):
                if 'HVO_' not in c:
                    continue
                essential_cog.add(c.strip())

    print('''
        total essenial arCOG: {0} - identified: {1} ({2}%)
        total essenial   COG: {3} - identified: {4} ({5}%)
        total essenial  both: {6} - identified: {7} ({8}%)
    '''.format(
        len(essential_arcog),
        len(essential_arcog & common_identified_prots),
        len(essential_arcog & common_identified_prots) / len(essential_arcog) * 100,
        len(essential_cog),
        len(essential_cog & common_identified_prots),
        len(essential_cog & common_identified_prots) / len(essential_cog) * 100,
        len(essential_cog|essential_arcog),
        len((essential_cog|essential_arcog) & common_identified_prots),
        len((essential_cog|essential_arcog) & common_identified_prots) / len(essential_cog|essential_arcog) * 100,
    ))

    not_found_small = set()
    for not_found_essential in essential_arcog - common_identified_prots:
        for org_prot in prot_info_dict.keys():
            if not_found_essential in org_prot:
                if prot_info_dict[org_prot]['molecular_weight'] < 15000:
                    not_found_small.add(not_found_essential)
    print('not found because small:', len(not_found_small))

    essential_haloferax = set()
    with open('essential_haloferax_failed_deletions.txt', 'r') as haloferax_in:
        for line in haloferax_in:
            essential_haloferax.add(line.strip())

    print('''
        total essenial Hv: {0} - identified: {1} ({2}%)
    '''.format(
        len(essential_haloferax),
        len(essential_haloferax & common_identified_prots),
        len(essential_haloferax & common_identified_prots) / len(essential_haloferax) * 100,
    ))

    not_found_small = set()
    for not_found_essential in essential_haloferax - common_identified_prots:
        for org_prot in prot_info_dict.keys():
            if not_found_essential in org_prot:
                if prot_info_dict[org_prot]['molecular_weight'] < 15000:
                    not_found_small.add(not_found_essential)
    print('not found because small:', len(not_found_small))


    with open('{0}_intersection_proteins_whole_proteome_datasets.csv'.format(0), 'w') as txt_out:
        for prot in prot_db_set - all_whole_cell_proteins:
            print(prot, file=txt_out)

    # datadict['x'].append('whole cell datasets <br> 0.5% FDR; >=2 PSMs')
    all_whole_cell_proteins_with_groups = set()
    for prot in all_whole_cell_proteins:
        for p in prot.split('<|>'):
            all_whole_cell_proteins_with_groups.add(p)
    whole_cell_missing = prot_db_set - all_whole_cell_proteins_with_groups
    # datadict['y'].append(len(whole_cell_missing))
    # datadict['text'].append(len(whole_cell_missing))
    # datadict['x'].append('all datasets <br> 0.5% FDR; >=2 PSMs')
    datadict['x'].insert(0, 0)
    all_proteins = copy.deepcopy(results_dict['all']['proteins']['safe_seq_num_spec_0005'])
    protein_groups = copy.deepcopy(results_dict['all']['protein_groups']['safe_seq_num_spec_0005'])
    for group in protein_groups:
        if group not in prot2cov.keys():
            continue
        else:
            for prot in group.split('<|>'):
                all_proteins.add(prot)
            all_proteins.add(group)
    all_missing = prot_db_set - all_proteins
    intersection_dict[0] = all_missing
    intersection_dict[8] = all_proteins
    # datadict['y'].append(len(all_missing))
    datadict['y'].insert(0, len(all_missing))
    # datadict['text'].append(len(all_missing))
    # pprint.pprint(whole_cell_missing - all_missing)
    # datadict['x'].append('all datasets <br> 1% FDR; >=1 PSM')
    unsafe_proteins = copy.deepcopy(results_dict['all']['proteins']['safe_seq'])
    protein_groups = copy.deepcopy(results_dict['all']['protein_groups']['safe_seq'])
    for group in protein_groups:
        if group not in prot2cov.keys():
            continue
        else:
            for prot in group.split('<|>'):
                unsafe_proteins.add(prot)
    # datadict['y'].append(len(prot_db_set - unsafe_proteins))
    # datadict['text'].append(len(prot_db_set - unsafe_proteins))


    arcog_dict = {}
    for l in arcog_letters + ['R', 'S', 'Not classified', 'total']:
        arcog_dict[l] = 0
    for prot in all_missing:
        funct_set = set()
        org_prot = prot.split('_NumberOfIdenticalSequences_')[0]
        for l in prot_info_dict[org_prot]['arcoglet'].split(';'):
            funct_set.add(l)
        for l in arcog_letters + ['R', 'S', 'total']:
            if l in funct_set:
                arcog_dict[l] += 1
            if l == 'total':
                arcog_dict[l] += 1
        if len(funct_set) == 1 and list(funct_set)[0] == '':
            arcog_dict['Not classified'] += 1
    return_list = fisher_test_with_correction(
        arcog_dict=arcog_dict,
        total_arcog_dict=total_arcog_dict,
        n=0,
    )
    datadict_list.insert(1,return_list[0])

    values = []
    labels = []
    for l in ['Not classified'] + arcog_letters + ['R', 'S']:
        if l in ['A', 'B', 'Y', 'Z', 'W']:
            continue
        values.append(arcog_dict[l])
        labels.append(l)
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        sort=False,
        showlegend=False,
    )])
    fig.update_traces(
        hoverinfo='label+percent+value',
        textinfo='none',
        textfont_size=20,
        marker=dict(
            colors=['#ffffff', '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000'],
            line=dict(
                color='#000000',
                width=0.5)
        ),
    )
    final_output_filename = 'dataset_overlap_pie_arcog_{0}.pdf'.format(0)
    plotly.io.write_image(fig, final_output_filename)
    plotly.io.write_html(fig, final_output_filename.replace('.pdf', '.html'))
    print('plotted: ', final_output_filename)

    layout = go.Layout(
        barmode='group',
        # bargap=0.5,
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=False,
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            dtick=1,
            title='Number of datasets<br>sharing protein identifications',
            # title = '',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickangle=0,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        yaxis=dict(
            showgrid=True,
            showline=True,
            ticks='outside',
            ticklen=3,
            tickwidth=0.25,
            dtick=250,
            title='Proteins',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            ),
            # tickangle=-90,
            # autorange = False,
            range = [0, 1200],
        ),
        showlegend=False,
        autosize=False,
        width=350,
        height=350,
        margin=go.layout.Margin(
            l=55,
            r=50,
            b=150,
            t=10,
        ),
    )

    plot_barplot(
        datadict_list=[datadict],
        output_filename=output_filename,
        layout=layout,
        color_dict=color_dict
    )
    fieldnames = [
        '#Datasets',
        '#Proteins'
    ]
    csv_out_name = output_filename.replace('.pdf', '.csv')
    with open(csv_out_name, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for n, x in enumerate(datadict['x']):
            out_dict = {
                '#Datasets': x,
                '#Proteins': datadict['y'][n],
            }
            csv_writer.writerow(out_dict)

    plot_barplot(
        datadict_list=[datadict_single],
        output_filename=output_filename.replace('.pdf', '_single.pdf'),
        layout=layout,
        color_dict=color_dict
    )

    layout = go.Layout(
        barmode='group',
        # orientation=-90,
        # bargap=0.5,
        xaxis=dict(
            showgrid=True,
            showline=True,
            zeroline=False,
            ticks='outside',
            ticklen=3,
            tickwidth=0.25,
            dtick=5,
            title='Percent of identified proteins',
            # title = '',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickangle=0,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        yaxis=dict(
            showgrid=False,
            showline=True,
            ticks='outside',
            ticklen=3,
            tickwidth=0.25,
            # dtick=250,
            title='arCOG category',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            ),
            # tickangle=-90,
            autorange = True,
            # range = [0, 1200],
        ),
        showlegend=True,
        autosize=False,
        width=600,
        height=1000,
        margin=go.layout.Margin(
            l=100,
            r=50,
            b=150,
            t=10,
        ),
    )

    plot_barplot(
        datadict_list=datadict_list,
        output_filename=output_filename.replace('.pdf', '_stats.pdf'),
        layout=layout,
        color_dict=color_dict,
        orientation='h',
    )


    arcog_dict = {}
    # all_arcog_letters = arcog_letters + ['R', 'S', 'Not classified']
    all_arcog_letters = arcog_letters_2
    for letter in all_arcog_letters:
        arcog_dict[letter] = {
            'total' : set()
        }
    for n, prot in enumerate(prot_info_dict):
        if 'spurious' in prot:
            continue
        arcoglet = prot_info_dict[prot]['arcoglet']
        arcog = arcoglet.split(';')
        for letter in arcog:
            if letter == '':
                letter = 'Not classified'
            arcog_dict[letter]['total'].add(prot)

    trace_list = []
    fieldnames = [
        'arCOGlet',
    ]
    out_dict_list = []
    for intersection in sorted(intersection_dict.keys()):
        fieldnames.append('{0} datasets, ratio identified/theoretical'.format(intersection))
        arcog_x_list = []
        arcog_y_list = []
        arcog_text_list = []
        for l in all_arcog_letters:
            if l in ['A', 'B', 'Y', 'Z', 'W']:
                continue
            total_arcog = len(arcog_dict[l]['total'])
            arcog_x_list.append(l)
            arcog_count = 0
            for prot_group in intersection_dict[intersection]:
                if '<|>' in prot_group:
                    inference = inference_dict.get(prot_group, None)
                    if inference != set():
                        continue
                for prot in prot_group.split('<|>'):
                    org_prot = prot.split('_NumberOfIdenticalSequences_')[0]
                    for x in prot_info_dict[org_prot]['arcoglet'].split(';'):
                        if l == x or (x == '' and l == 'Not classified'):
                            arcog_count += 1
                            break
            try:
                arcog_y_list.append(arcog_count/total_arcog) 
            except:
                print(l, arcog_count, total_arcog)
            arcog_text_list.append(
                '{0} of {1}'.format(arcog_count, total_arcog)
            )
        trace_scatter = go.Scatter(
            x = arcog_x_list,
            y = arcog_y_list,
            text = arcog_text_list,
            name = intersection,
            mode='markers',
            marker=dict(
                color = color_dict['overlap'][intersection],
                symbol=marker_dict['overlap'][intersection],
                size=5,
            ),
        )
        trace_list.append(trace_scatter)
        for n, cat_bin in enumerate(arcog_x_list):
            try:
                out_dict = out_dict_list[n]
                out_dict['{0} datasets, ratio identified/theoretical'.format(intersection)] = arcog_y_list[n]
                out_dict_list[n] = out_dict
            except:
                out_dict = {}
                out_dict['arCOGlet'] = cat_bin
                out_dict['{0} datasets, ratio identified/theoretical'.format(intersection)] = arcog_y_list[n]
                out_dict_list.append(out_dict)

    csv_out_name = output_filename.replace('.pdf', '_scatter_{0}.csv'.format('arCOGlet'))
    with open(csv_out_name, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for l in out_dict_list:
            csv_writer.writerow(l)

    layout = go.Layout(
        legend=dict(
            orientation= "v",
            font=dict(
                family="sans-serif",
                size=10,
                color="black"
            ),
        ),
        yaxis=dict(
            title='Identified/theoretical proteins',
            showline=True,
            zeroline=True,
            linecolor='rgb(0, 0, 0)',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)',
            ),
            ticks='outside',
            ticklen=2,
            # dtick=0.25,
            tickwidth=0.25,
            tickcolor='rgb(0, 0, 0)',
            # ticksuffix='k',
            # tickangle=-90,
            side='left',
            autorange=False,
            range=[0,1],
        ),
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=False,
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            title='arCOG single letter code',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickangle=-90,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        showlegend=False,
        autosize=False,
        width=525,
        height=400,
        margin=go.layout.Margin(
            l=65,
            r=55,
            b=150,
            t=40,
        ),
    )

    fig = go.Figure(data=trace_list, layout=layout)
    plotly.io.write_image(
        fig,
        output_filename.replace('.pdf', '_{0}_scatter.pdf'.format('arcoglet'))
    )
    plotly.io.write_html(
        fig,
        output_filename.replace('.pdf', '_{0}_scatter.html'.format('arcoglet'))
    )


def fisher_test_with_correction(arcog_dict={}, total_arcog_dict={}, n=0, return_total=False):
    '''
    Performs a Fisher test to test if two given arCOG letter distributions are significantly different.
    A Bonferroni correction is performed to correct for multiple testing.
    '''
    fisher_test_dict = {}
    p_value_list = []
    datadict = {
        'x': [],
        'y': [],
        'text': [],
        'name': '',
    }
    datadict_total = {
        'x': [],
        'y': [],
        'text': [],
        'name': '',
    }
    for l in arcog_letters + ['R', 'S', 'Not classified']:
        fisher_test_dict[l] = {
            'arCOGlet': l,
            'No Proteins Intersection {0}'.format(n): arcog_dict[l],
            'No Proteins Total DB': total_arcog_dict[l],
            'p-value': 1,
            'Corrected p-value': 1, 
        }
        oddsratio, pvalue = stats.fisher_exact(
            [
                [arcog_dict[l],
                arcog_dict['total'] - arcog_dict[l]],
                [total_arcog_dict[l],
                total_arcog_dict['total'] - total_arcog_dict[l]]
            ]
        )
        fisher_test_dict[l]['p-value'] = pvalue
        p_value_list.append(pvalue)
    corrected_p_values = multipletests(p_value_list, method='fdr_bh')[1]
    for i, l in enumerate(arcog_letters + ['R', 'S', 'Not classified']):
        fisher_test_dict[l]['Corrected p-value'] = corrected_p_values[i]

    fieldnames = [
        'arCOGlet',
        'No Proteins Intersection {0}'.format(n),
        'No Proteins Total DB',
        'p-value',
        'Corrected p-value', 
    ]
    with open('{0}_intersection_arcog_stats.csv'.format(n), 'w') as stats_out:
        csv_writer = csv.DictWriter(stats_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        letters = arcog_letters + ['R', 'S', 'Not classified']
        for l in reversed(letters):
            if l in ['W', 'Z', 'Y', 'B', 'A']:
                continue
            out_dict = fisher_test_dict[l]
            csv_writer.writerow(out_dict)
            datadict['y'].append(l)
            datadict['x'].append(
                arcog_dict[l]/arcog_dict['total']*100
            )
            datadict['name'] = n
            if return_total:
                datadict_total['y'].append(l)
                datadict_total['x'].append(
                    total_arcog_dict[l]/total_arcog_dict['total']*100
                )
                datadict_total['name'] = 'Total proteome'

    return [datadict, datadict_total]

def plot_prot_coverage(results_dict, output_filename=None, inference_dict=None, all_datasets=False):
    '''
    Plot a boxplot for the protein sequence coverage of all proteins identified in the different data sets.
    In addition, a scatter plot for the total proteome sequence coverage is overlaying the boxplot.
    '''
    trace_list = []
    prot_db_dict = {}
    fasta_dict = {}
    N_termini = set()
    C_termini = set()

    with open(db_fasta, 'r') as db_input:
        total_proteome_length = 0
        for fasta_id, db_sequence in ursgal.ucore.parse_fasta(db_input):
            if 'decoy_' in fasta_id or 'spurious' in fasta_id:
                continue
            prot_id = fasta_id.split(' ')[0]
            if 'HVO' not in prot_id:
                continue
            prot_db_dict[fasta_id] = len(db_sequence)
            total_proteome_length += len(db_sequence)
            fasta_dict[fasta_id] = db_sequence
        print('Total proteome:', total_proteome_length)

    prot_db_set = set(prot_db_dict.keys())
    prot2cov = {}
    x_scatter = []
    y_scatter = []
    fieldnames = [
        'PRIDE ID',
        'Protein sequence coverage (median)',
        'Protein sequence coverage (lower quartile)',
        'Protein sequence coverage (upper quartile)',
        'Total proteome sequence coverage',
    ]
    csv_out_name = output_filename.replace('.pdf', '.csv')
    with open(csv_out_name, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for PRIDE_ID in PRIDE_order:
            out_dict = {}
            out_dict['PRIDE ID'] = PRIDE_ID
            if all_datasets:
                prot2cov[PRIDE_ID] = {}
            print(PRIDE_ID)
            if PRIDE_ID == 'all':
                name = 'All datasets'
                with open('0_intersection_proteins_all_datasets.csv', 'w') as txt_out:
                    nowhere_found = (prot_db_set - results_dict[PRIDE_ID]['proteins']['safe_seq_num_spec_0005'])
                    for prot in nowhere_found:
                        print(prot, file=txt_out)
            else:
                name = PRIDE_ID 
            x_scatter.append(name)

            y_box = []
            proteins = results_dict[PRIDE_ID]['proteins']['safe_seq_num_spec_0005']
            protein_groups = results_dict[PRIDE_ID]['protein_groups']['safe_seq_num_spec_0005']
            protein2group = {}
            proteome_aas = 0
            unique_map = 0
            clear_map = 0
            ambiguous_map = 0
            for group in protein_groups:
                if group not in inference_dict.keys() and PRIDE_ID == 'all':
                    continue
                if inference_dict[group][PRIDE_ID] == set():
                    if PRIDE_ID == 'all':
                        unique_map += 1
                    pass
                elif len(inference_dict[group][PRIDE_ID]) == 1:
                    protein = sorted(inference_dict[group][PRIDE_ID])[0]
                    if protein not in protein2group.keys():
                        protein2group[protein] = set()
                    protein2group[protein].add(group)
                    if PRIDE_ID == 'all':
                        clear_map += 1
                    continue
                else:
                    if PRIDE_ID == 'all':
                        ambiguous_map += 1
                    continue

                group_prot_lengths = []
                for prot in group.split('<|>'):
                    group_prot_lengths.append(prot_db_dict[prot])
                seq_position = set()
                for seq in results_dict[PRIDE_ID]['protein_dict'][group].keys():
                    if seq not in results_dict[PRIDE_ID]['peptides']['safe_num_specs']:
                        continue
                    start_group, stop_group = results_dict[PRIDE_ID]['protein_dict'][group][seq]['start_stop']
                    start = start_group.split('<|>')[0]
                    stop = stop_group.split('<|>')[0]
                    for aa in range(int(start), int(stop)+1):
                        seq_position.add(aa)
                protein_cov = 100*len(seq_position)/statistics.median(group_prot_lengths)
                y_box.append(protein_cov)
                proteome_aas += len(seq_position)
                if all_datasets:
                    prot2cov[PRIDE_ID][prot] = [protein_cov, seq_position]
                if PRIDE_ID == 'all':
                    prot2cov[group] = [protein_cov, seq_position]

            for prot in proteins:
                seq_position = set()
                start_stop_set = set()
                if prot in protein2group.keys():
                    for group in protein2group[prot]:
                        prot_index = group.split('<|>').index(prot)
                        for seq in results_dict[PRIDE_ID]['protein_dict'][group].keys():
                            if seq in ['datasets', 'prot_q_value_psm', 'prot_q_value_seq', 'samples']:
                                continue
                            start_group, stop_group = results_dict[PRIDE_ID]['protein_dict'][group][seq]['start_stop']
                            start_stop_set.add(
                                (
                                    start_group.split('<|>')[prot_index],
                                    stop_group.split('<|>')[prot_index]
                                )
                            )
                seq_set = set(results_dict[PRIDE_ID]['protein_dict'][prot].keys())
                for seq in seq_set:
                    if seq in ['datasets', 'prot_q_value_psm', 'prot_q_value_seq', 'samples']:
                        continue
                    if seq not in results_dict[PRIDE_ID]['peptides']['safe_num_specs']:
                        continue
                    start, stop = results_dict[PRIDE_ID]['protein_dict'][prot][seq]['start_stop']
                    start_stop_set.add((start, stop))
                for start, stop in start_stop_set:
                    try:
                        start = int(start)
                        stop = int(stop)
                    except:
                        start = int(start.split(';')[0])
                        stop = int(stop.split(';')[0])
                    for aa in range(start, stop+1):
                        seq_position.add(aa)
                    if start <= 2 and PRIDE_ID == 'all':
                        N_termini.add(prot)
                    if stop == prot_db_dict[prot] and PRIDE_ID == 'all':
                        C_termini.add(prot)
                protein_cov = 100*len(seq_position)/prot_db_dict[prot]
                y_box.append(protein_cov)
                proteome_aas += len(seq_position)
                if all_datasets:
                    prot2cov[PRIDE_ID][prot] = [protein_cov, seq_position]
                if PRIDE_ID == 'all':
                    prot2cov[prot] = [protein_cov, seq_position]

            y_scatter.append(100*proteome_aas/total_proteome_length)
            out_dict['Total proteome sequence coverage'] = proteome_aas/total_proteome_length

            trace = go.Box(
                y=y_box,
                name=name,
                marker_color = '#2ca02c'
            )
            print(
                'median {0}:'.format(PRIDE_ID),
                statistics.median(y_box)
            )
            out_dict['Protein sequence coverage (median)'] = np.quantile(y_box, 0.5)
            out_dict['Protein sequence coverage (lower quartile)'] = np.quantile(y_box, 0.25)
            out_dict['Protein sequence coverage (upper quartile)'] = np.quantile(y_box, 0.75)
            csv_writer.writerow(out_dict)
            trace_list.append(trace)

    trace_scatter = go.Scatter(
        x = x_scatter,
        y = y_scatter,
        name = 'Proteome coverage',
        mode='markers',
        marker=dict(
            color = 'rgb(0,0,0)',
            symbol='cross',
            size=5,
        ),
        # yaxis='y2'
    )

    trace_list.append(trace_scatter)

    layout = go.Layout(
        legend=dict(orientation= "h"),
        yaxis=dict(
            title='Sequence coverage in %',
            showline=True,
            linecolor='rgb(0, 0, 0)',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)',
            ),
            ticks='outside',
            ticklen=2,
            # dtick=0.25,
            tickwidth=0.25,
            tickcolor='rgb(0, 0, 0)',
            # ticksuffix='k',
            # tickangle=-90,
            side='left',
            autorange=False,
            range=[0,101]
        ),
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=True,
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            title='',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickangle=-90,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        showlegend=False,
        autosize=False,
        width=600,
        height=400,
        margin=go.layout.Margin(
            l=55,
            r=50,
            b=150,
            t=10,
        ),
    )

    print('N-term:', len(N_termini))
    print('C-term:', len(C_termini))
    fig = go.Figure(data=trace_list, layout=layout)
    plotly.io.write_image(fig, output_filename)

    return prot2cov

def plot_total_peptides_proteins(results_dict, inference_dict=None, output_filename='test.pdf'):
    '''
    Plot a barplot with the total number of identified peptide sequences and proteins for all data sets
    '''
    trace_list = []
    color_dict = {
        'safe': '#0b5d65',
        # 'safe_num_specs': 'rgb(4,90,141)',
        'safe_num_specs': '#20b7c6',
        'safe_seq': '#175317',
        'safe_seq_num_spec': 'rgb(44,162,95)',
        # 'safe_seq_num_spec_0005': 'rgb(0,109,44)', #'#2ca02c',
        'safe_seq_num_spec_0005': '#2ca02c',
    }
    name_dict = {
        'safe': 'Peptide sequences 1% FDR',
        'safe_num_specs': 'Peptide sequences 1% FDR, >1 PSM',
        'safe_seq': 'Proteins 1% FDR',
        'safe_seq_num_spec': 'Proteins 1% FDR, >1 PSM',
        'safe_seq_num_spec_0005': 'Proteins 0.5% FDR, >1PSM',
    }
    fieldnames = [
        'PRIDE ID',
        '#Peptide sequences (peptide q-value <= 1%)',
        '#Peptide sequences (peptide q-value <= 1% and #PSMs >= 2)',
        '#Proteins (protein q-value <= 0.5%)',
        '#Proteins (protein q-value <= 0.5% and #PSMs >= 2)',
    ]
    csv_out_name = output_filename.replace('.pdf', '.csv')
    with open(csv_out_name, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for level in['safe', 'safe_num_specs']:
        # for level in['safe_num_specs']:
            out_dict = {}
            out_dict['PRIDE ID'] = level
            csv_writer.writerow(out_dict)
            x_seq = []
            y_seq = []
            text_seq = []
            for PRIDE_ID in PRIDE_order:
                if PRIDE_ID == 'all':
                    y_name = 'All datasets'
                else:
                    y_name = PRIDE_ID
                out_dict['PRIDE ID'] = PRIDE_ID
                x_seq.append(y_name)
                y_seq.append(len(results_dict[PRIDE_ID]['peptides'][level])/1000)
                # text_seq.append(len(results_dict[PRIDE_ID]['peptides'][level]))
                if level == 'safe':
                    out_dict['#Peptide sequences (peptide q-value <= 1%)'] = len(
                        results_dict[PRIDE_ID]['peptides'][level]
                    )
                else:
                    out_dict['#Peptide sequences (peptide q-value <= 1% and #PSMs >= 2)'] = len(
                        results_dict[PRIDE_ID]['peptides'][level]
                    )
                csv_writer.writerow(out_dict)
            
            trace_bar = go.Bar(
                x = x_seq,
                y = y_seq,
                text = text_seq,
                name = name_dict[level],
                textposition = 'outside',
                textangle = -90,
                outsidetextfont=dict(
                    size=9,
                    color='rgb(0,0,0)'
                ),
                marker=dict(
                    color=color_dict[level],
                    # line=dict(
                    #     color='rgb(8,48,107)',
                    #     width=1.5,
                    # )
                )
             )
            trace_list.append(trace_bar)

        trace_bar = go.Bar(
            x = x_seq,
            y = [0 for x in x_seq],
            text = [],
            name = '',
            marker=dict(
                color=color_dict[level],
            )
         )
        trace_list.append(trace_bar)
        # trace_bar = go.Bar(
        #     x = x_seq,
        #     y = [0 for x in x_seq],
        #     text = [],
        #     name = '',
        #     marker=dict(
        #         color=color_dict[level],
        #     )
        #  )
        # trace_list.append(trace_bar)
        # # trace_bar = go.Bar(
        # #     x = x_seq,
        # #     y = [0 for x in x_seq],
        # #     text = [],
        # #     name = '',
        # #     marker=dict(
        # #         color=color_dict[level],
        # #     )
        # #  )
        # # trace_list.append(trace_bar)


        trace_bar = go.Bar(
            x = x_seq,
            y = [0 for x in x_seq],
            text = [],
            name = '',
            marker=dict(
                color=color_dict[level],
            ),
            yaxis='y2'
         )
        trace_list.append(trace_bar)
        # trace_bar = go.Bar(
        #     x = x_seq,
        #     y = [0 for x in x_seq],
        #     text = [],
        #     name = '',
        #     marker=dict(
        #         color=color_dict[level],
        #     ),
        #     yaxis='y2'
        #  )
        # trace_list.append(trace_bar)

        for level in ['safe_seq', 'safe_seq_num_spec_0005']:
        # for level in ['safe_seq_num_spec_0005']:
            print('protein level:', level)
            inference_dict = protein_inference(results_dict, safe_word=level)
            x_prot = []
            y_prot = []
            text_prot = []
            maps = 0
            no_map = 0
            multimap = 0
            for PRIDE_ID in PRIDE_order:
                if PRIDE_ID == 'all':
                    y_name = 'All datasets'
                else:
                    y_name = PRIDE_ID
                out_dict['PRIDE ID'] = PRIDE_ID
                proteins = results_dict[PRIDE_ID]['proteins'][level]
                protein_groups = results_dict[PRIDE_ID]['protein_groups'][level]

                for group in protein_groups:
                    if group not in inference_dict.keys() and PRIDE_ID == 'all':
                        continue
                    if inference_dict[group][PRIDE_ID] == set():
                        proteins.add(group)
                        if PRIDE_ID == 'all':
                            no_map += 1
                    elif len(inference_dict[group][PRIDE_ID]) >= 2:
                        if PRIDE_ID == 'all':
                            multimap += 1
                    else:
                        if PRIDE_ID == 'all':
                            maps += 1
                        continue   

                x_prot.append(y_name)
                y_prot.append(len(proteins)/1000)
                # text_prot.append(len(proteins))
                if level == 'safe_seq':
                    out_dict['#Proteins (protein q-value <= 0.5%)'] = len(
                        proteins
                    )
                else:
                    out_dict['#Proteins (protein q-value <= 0.5% and #PSMs >= 2)'] = len(
                        proteins
                    )
                csv_writer.writerow(out_dict)

            print('groups that map:', maps)
            print('not mapping groups:', no_map)
            print('mapping to multiple prots:', multimap)

            trace_bar = go.Bar(
                x = x_prot,
                y = y_prot,
                text = text_prot,
                name = name_dict[level],
                textposition = 'outside',
                textangle = -90,
                outsidetextfont=dict(
                    size=9,
                    color='rgb(0,0,0)'
                ),
                marker=dict(
                    color=color_dict[level],
                    # line=dict(
                    #     color='rgb(8,48,107)',
                    #     width=1.5,
                    # )
                ),
                yaxis='y2'
             )
            trace_list.append(trace_bar)

    layout = go.Layout(
        barmode='group',
        yaxis=dict(
            title='Peptide sequences',
            # rangemode='tozero',
            showline=True,
            linecolor='#20b7c6',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='#20b7c6',
            ),
            ticks='outside',
            ticklen=2,
            # dtick=0.25,
            tickwidth=0.25,
            tickcolor='#20b7c6',
            ticksuffix='k',
            # tickangle=-90,
            side='left',
            # autorange=False,
            range=[0,50]
        ),
        yaxis2=dict(
            title='Proteins',
            # rangemode='tozero',
            showgrid=False,
            showline=True,
            linecolor='#2ca02c',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='#2ca02c'
            ),
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            tickcolor='#2ca02c',
            overlaying='y',
            side='right',
            autorange=False,
            range=[0,3.5],
            ticksuffix='k',
        ),
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=True,
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            title='',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickangle=-90,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        showlegend=False,
        legend=dict(
            orientation= "v",
        ),
        autosize=False,
        width=600,
        height=400,
        margin=go.layout.Margin(
            l=55,
            r=50,
            b=150,
            t=10,
        ),
    )

    fig = go.Figure(data=trace_list, layout=layout)
    plotly.io.write_image(fig, output_filename)

def plot_prot_properties(
    results_dict,
    inference_dict=None,
    prot_info_dict=None,
    # peptide_info_dict=None,
    output_filename='test.pdf'
):
    '''
    Plot box plots and scatter plots for different protein
    properties (length, hydrophobicity, pI).
    Barplots show the distribution of identified proteins for each dataset,
    while scatterplots indicate, for each dataset, the percentage of theoretical
    proteins that have been identified within equally sized bins. 
    '''

    box_plot_dict = {
        'molecular_weight': [],
        'pi': [],
        'hydrophobicity': [],
    }
    y_axis_dict = {
        'name': {
            'molecular_weight': 'Molecular weight in Da',
            'pi': 'pI',
            'hydrophobicity': 'Hydrophobicity',
            'arcoglet': 'arCOG one letter code'
        },
        'range':{
            'molecular_weight': [0, 240000],
            'pi': [2.5, 12.5],
            'hydrophobicity': [-1.8, 2.2],
        }
    }
    for cat in box_plot_dict.keys():
        fieldnames = [
            'PRIDE ID',
            '(median)',
            '(lower quartile)',
            '(upper quartile)',
        ]
        csv_out_name = output_filename.replace('.csv', '_box_{0}.csv'.format(cat))
        with open(csv_out_name, 'w') as csv_out:
            csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
            csv_writer.writeheader()
            out_dict = {}
            for PRIDE_ID in PRIDE_order:
                out_dict['PRIDE ID'] = PRIDE_ID
                if PRIDE_ID == 'all':
                    name = 'All datasets'
                else:
                    name = PRIDE_ID
                proteins = results_dict[PRIDE_ID]['proteins']['safe_seq_num_spec_0005']
                protein_groups = results_dict[PRIDE_ID]['protein_groups']['safe_seq_num_spec_0005']

                y_box = []
                text = []
                for group in protein_groups:
                    inference = inference_dict.get(group, None)
                    group_tmp = []
                    if inference == set():
                        for prot in group.split('<|>'):
                            org_prot = prot.split('_NumberOfIdenticalSequences_')[0]
                            group_tmp.append(prot_info_dict[org_prot][cat])
                    else:
                        continue
                    mean = statistics.mean(group_tmp)
                    y_box.append(mean)
                    text,append(group)

                for protein in proteins:
                    org_prot = protein.split('_NumberOfIdenticalSequences_')[0]
                    y_box.append(prot_info_dict[org_prot][cat])
                    text.append(protein)

                trace = go.Box(
                    y=y_box,
                    name=name,
                    text=text,
                    marker_color = '#2ca02c',
                    jitter=0.5,
                    pointpos=-1.8,
                    boxpoints='all',
                    marker=dict(
                        color = '#2ca02c',
                        symbol='cross',
                        size=1,
                    ),
                    line_width=0.5
                )
                box_plot_dict[cat].append(trace)
                out_dict['(median)'] = np.quantile(y_box, 0.5)
                out_dict['(lower quartile)'] = np.quantile(y_box, 0.25)
                out_dict['(upper quartile)'] = np.quantile(y_box, 0.75)
                csv_writer.writerow(out_dict)

        name = 'All proteins'
        y_box = []
        text = []
        protein_dict_list = []
        for protein in prot_info_dict.keys():
            if 'spurious' in protein:
                continue
            y_box.append(prot_info_dict[protein][cat])
            text.append(protein)
            protein_dict_list.append({
                'protein': protein,
                cat: prot_info_dict[protein][cat],
                'arcoglet': prot_info_dict[protein]['arcoglet'],
            })

        trace = go.Box(
            y=y_box,
            name=name,
            marker_color = '#2ca02c',
            text=text,
            jitter=0.5,
            pointpos=-1.8,
            boxpoints='all',
            marker=dict(
                color = '#2ca02c',
                symbol='cross',
                size=1,
            ),
            line_width=0.5
        )
        box_plot_dict[cat].append(trace)


        layout = go.Layout(
            legend=dict(orientation= "h"),
            yaxis=dict(
                title=y_axis_dict['name'][cat],
                showline=True,
                zeroline=False,
                linecolor='rgb(0, 0, 0)',
                titlefont=dict(
                    size=16,
                    color='rgb(0, 0, 0)'
                ),
                tickfont=dict(
                    size=14,
                    color='rgb(0, 0, 0)',
                ),
                ticks='outside',
                ticklen=2,
                # dtick=0.25,
                tickwidth=0.25,
                tickcolor='rgb(0, 0, 0)',
                # ticksuffix='k',
                # tickangle=-90,
                side='left',
                autorange=False,
                range=y_axis_dict['range'][cat]
            ),
            xaxis=dict(
                showgrid=False,
                showline=True,
                zeroline=False,
                ticks='outside',
                ticklen=2,
                tickwidth=0.25,
                title='',
                titlefont=dict(
                    size=16,
                    color='rgb(0, 0, 0)'
                ),
                tickangle=-90,
                tickfont=dict(
                    size=14,
                    color='rgb(0, 0, 0)'
                )
            ),
            showlegend=False,
            autosize=False,
            width=500,
            height=350,
            margin=go.layout.Margin(
                l=65,
                r=55,
                # b=200
                b=150,
                t=10,
                # pad=4
            ),
        )

        fig = go.Figure(data=box_plot_dict[cat], layout=layout)
        plotly.io.write_image(
            fig,
            output_filename.replace('.csv', '_{0}_no_spurious.pdf'.format(cat))
        )
        plotly.io.write_html(
            fig,
            output_filename.replace('.csv', '_{0}_no_spurious.html'.format(cat))
        )

        y_axis_lists = {}
        x_axis_list = []
        for n, protein_dict in enumerate(sorted(protein_dict_list, key=itemgetter(cat))):
            if n % 206 == 0 or n == len(protein_dict_list)-1:
                if n != 0:
                    max_v = protein_dict[cat]
                    if cat != 'molecular_weight':
                        x_axis_list.append('{0:.2f} to {1:.2f}'.format(
                            round(min_v, 2),
                            round(max_v, 2),
                        ))
                    else:
                        x_axis_list.append('{0} to {1}'.format(
                            min_v,
                            max_v,
                        ))
                    for PRIDE_ID in pract.keys():
                        y_axis_lists[PRIDE_ID].append(
                            len(pract[PRIDE_ID])/len(theo)
                        )
                min_v = protein_dict[cat]
                theo = []
                pract = {}
                print('>>>>> n =', n)
            protein = protein_dict['protein']
            theo.append(protein)
            for PRIDE_ID in PRIDE_order:
                if PRIDE_ID not in pract.keys():
                    pract[PRIDE_ID] = set()
                if PRIDE_ID not in y_axis_lists.keys():
                    y_axis_lists[PRIDE_ID] = []
                proteins = results_dict[PRIDE_ID]['proteins']['safe_seq_num_spec_0005']
                protein_groups = results_dict[PRIDE_ID]['protein_groups']['safe_seq_num_spec_0005']
                for group in protein_groups:
                    inference = inference_dict.get(group, None)
                    if inference == set():
                        if protein in group:
                            pract[PRIDE_ID].add(protein)
                for p in proteins:
                    if protein in p:
                        pract[PRIDE_ID].add(protein)

        if cat == 'molecular_weight':
            arcog_dict = {}
            all_arcog_letters = arcog_letters + ['R', 'S', 'Not classified']
            for letter in all_arcog_letters:
                arcog_dict[letter] = {
                    'total' : set()
                }
            print('total identified:', len(protein_dict_list))
            for n, protein_dict in enumerate(protein_dict_list):
                protein = protein_dict['protein']
                arcog = protein_dict['arcoglet'].split(';')
                for letter in arcog:
                    if letter == '':
                        letter = 'Not classified'
                    for PRIDE_ID in PRIDE_order:
                        if PRIDE_ID not in arcog_dict[letter].keys():
                            arcog_dict[letter][PRIDE_ID] = set()
                        proteins = results_dict[PRIDE_ID]['proteins']['safe_seq_num_spec_0005']
                        protein_groups = results_dict[PRIDE_ID]['protein_groups']['safe_seq_num_spec_0005']
                        for group in protein_groups:
                            inference = inference_dict.get(group, None)
                            if inference == set():
                                if protein in group:
                                    arcog_dict[letter][PRIDE_ID].add(protein)
                                    # arcog_dict[letter]['total'].add(protein)
                        for p in proteins:
                            if protein in p:
                                arcog_dict[letter][PRIDE_ID].add(protein)
                        arcog_dict[letter]['total'].add(protein)

            total_arcog = 0
            trace_list = []
            fieldnames = [
                'arCOGlet',
            ]
            out_dict_list = []
            for PRIDE_ID in PRIDE_order:
                fieldnames.append('{0} ratio identified/theoretical'.format(PRIDE_ID))
                if PRIDE_ID == 'all':
                    name = 'All datasets'
                    for l in all_arcog_letters:
                        total_arcog += len(arcog_dict[l]['total'])
                else:
                    name = PRIDE_ID
                arcog_y_list = []
                arcog_text_list = []
                for l in all_arcog_letters:
                    pride_arcog = len(arcog_dict[l].get(PRIDE_ID, []))
                    if pride_arcog == 0:
                        arcog_y_list.append(0)
                    else:
                        arcog_y_list.append(pride_arcog/len(arcog_dict[l]['total'])) 
                    arcog_text_list.append(
                        '{0} of {1}'.format(pride_arcog, len(arcog_dict[l]['total']))
                    )
                trace_scatter = go.Scatter(
                    x = all_arcog_letters,
                    y = arcog_y_list,
                    text = arcog_text_list,
                    name = name,
                    mode='markers',
                    marker=dict(
                        color = PRIDE_COLORS[PRIDE_ID],
                        symbol='square',
                        size=3,
                    ),
                )
                trace_list.append(trace_scatter)
                for n, cat_bin in enumerate(all_arcog_letters):
                    try:
                        out_dict = out_dict_list[n]
                        out_dict['{0} ratio identified/theoretical'.format(PRIDE_ID)] = arcog_y_list[n]
                        out_dict_list[n] = out_dict
                    except:
                        out_dict = {}
                        out_dict['arCOGlet'] = cat_bin
                        out_dict['{0} ratio identified/theoretical'.format(PRIDE_ID)] = arcog_y_list[n]
                        out_dict_list.append(out_dict)

            csv_out_name = output_filename.replace('.csv', '_scatter_{0}.csv'.format('arCOGlet'))
            with open(csv_out_name, 'w') as csv_out:
                csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
                csv_writer.writeheader()
                for l in out_dict_list:
                    csv_writer.writerow(l)

            print('total arcog:', total_arcog)
            layout = go.Layout(
                legend=dict(
                    orientation= "v",
                    font=dict(
                        family="sans-serif",
                        size=10,
                        color="black"
                    ),
                ),
                yaxis=dict(
                    title='Identified/theoretical proteins',
                    showline=True,
                    zeroline=True,
                    linecolor='rgb(0, 0, 0)',
                    titlefont=dict(
                        size=16,
                        color='rgb(0, 0, 0)'
                    ),
                    tickfont=dict(
                        size=14,
                        color='rgb(0, 0, 0)',
                    ),
                    ticks='outside',
                    ticklen=2,
                    # dtick=0.25,
                    tickwidth=0.25,
                    tickcolor='rgb(0, 0, 0)',
                    # ticksuffix='k',
                    # tickangle=-90,
                    side='left',
                    autorange=False,
                    range=[0,1],
                ),
                xaxis=dict(
                    showgrid=False,
                    showline=False,
                    zeroline=False,
                    ticks='outside',
                    ticklen=2,
                    tickwidth=0.25,
                    title=y_axis_dict['name']['arcoglet'],
                    titlefont=dict(
                        size=16,
                        color='rgb(0, 0, 0)'
                    ),
                    tickangle=-90,
                    tickfont=dict(
                        size=14,
                        color='rgb(0, 0, 0)'
                    )
                ),
                showlegend=True,
                autosize=False,
                width=750,
                height=350,
                margin=go.layout.Margin(
                    l=65,
                    r=55,
                    b=150,
                    t=40,
                ),
            )

            fig = go.Figure(data=trace_list, layout=layout)
            plotly.io.write_image(
                fig,
                output_filename.replace('.csv', '_{0}_scatter_no_spurious.pdf'.format('arcoglet'))
            )
            plotly.io.write_html(
                fig,
                output_filename.replace('.csv', '_{0}_scatter_no_spurious.html'.format('arcoglet'))
            )

        trace_list = []
        fieldnames = [
            cat,
        ]
        out_dict_list = []
        for PRIDE_ID in PRIDE_order:
            fieldnames.append('{0} ratio identified/theoretical'.format(PRIDE_ID))
            if PRIDE_ID == 'all':
                name = 'All datasets'
            else:
                name = PRIDE_ID
            trace_scatter = go.Scatter(
                x = x_axis_list,
                y = y_axis_lists[PRIDE_ID],
                name = name,
                mode='markers',
                marker=dict(
                    color = PRIDE_COLORS[PRIDE_ID],
                    symbol='square',
                    size=3,
                ),
            )
            trace_list.append(trace_scatter)
            for n, cat_bin in enumerate(x_axis_list):
                try:
                    out_dict = out_dict_list[n]
                    out_dict['{0} ratio identified/theoretical'.format(PRIDE_ID)] = y_axis_lists[PRIDE_ID][n]
                    out_dict_list[n] = out_dict
                except:
                    out_dict = {}
                    out_dict[cat] = cat_bin
                    out_dict['{0} ratio identified/theoretical'.format(PRIDE_ID)] = y_axis_lists[PRIDE_ID][n]
                    out_dict_list.append(out_dict)
        
        csv_out_name = output_filename.replace('.csv', '_scatter_{0}.csv'.format(cat))
        with open(csv_out_name, 'w') as csv_out:
            csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
            csv_writer.writeheader()
            for l in out_dict_list:
                csv_writer.writerow(l)

        layout = go.Layout(
            legend=dict(
                orientation= "h",
                font=dict(
                    family="sans-serif",
                    size=10,
                    color="black"
                ),
            ),
            yaxis=dict(
                title='Identified/theoretical proteins',
                showline=True,
                zeroline=True,
                linecolor='rgb(0, 0, 0)',
                titlefont=dict(
                    size=16,
                    color='rgb(0, 0, 0)'
                ),
                tickfont=dict(
                    size=14,
                    color='rgb(0, 0, 0)',
                ),
                ticks='outside',
                ticklen=2,
                # dtick=0.25,
                tickwidth=0.25,
                tickcolor='rgb(0, 0, 0)',
                # ticksuffix='k',
                # tickangle=-90,
                side='left',
                autorange=False,
                range=[0,1],
            ),
            xaxis=dict(
                showgrid=False,
                showline=False,
                zeroline=False,
                ticks='outside',
                ticklen=2,
                tickwidth=0.25,
                title=y_axis_dict['name'][cat],
                titlefont=dict(
                    size=16,
                    color='rgb(0, 0, 0)'
                ),
                tickangle=-90,
                tickfont=dict(
                    size=14,
                    color='rgb(0, 0, 0)'
                )
            ),
            showlegend=True,
            autosize=False,
            width=500,
            height=350,
            margin=go.layout.Margin(
                l=65,
                r=55,
                b=150,
                t=40,
            ),
        )

        fig = go.Figure(data=trace_list, layout=layout)
        plotly.io.write_image(
            fig,
            output_filename.replace('.csv', '_{0}_scatter_no_spurious.pdf'.format(cat))
        )
        plotly.io.write_html(
            fig,
            output_filename.replace('.csv', '_{0}_scatter_no_spurious.html'.format(cat))
        )

def plot_peptide_distributions(
    results_dict,
    peptide_info_dict=None,
    # peptide_info_dict=None,
    output_filename='test.pdf'
):
    '''
    Plot the number of identified peptide sequences in bins across different
    properties (length, hydrophobicity, pI).
    Barplots are plotted with one trace for peptide q-value <= 1% and one trace
    for q-value <= 1% and #PSMs >= 2.
    '''
    y_axis_dict = {
        'name': {
            'length': 'Peptide length',
            'pi': 'pI',
            'hydrophobicity': 'Hydrophobicity',
        },
        'range':{
            'length': [5, 55],
            'pi': [2.5, 12.5],
            'hydrophobicity': [-1.8, 2.2],
        }
    }
    for cat in y_axis_dict['name'].keys():
        trace_list = []
        bins = []
        bin_elements = {}
        if cat == 'length':
            for x in range(6,51):
                bins.append(x)
            for n, b in enumerate(bins):
                if n == len(bins)-1:
                    continue
                bin_elements[b] = []
                for element in range(b, bins[n+1]):
                    bin_elements[b].append(element)
        elif cat == 'pi':
            for x in np.linspace(2.8, 13.04, 51):
                bins.append(round(x, 2))
            for n, b in enumerate(bins):
                if n == len(bins)-1:
                    continue
                bin_elements['{0}-{1}'.format(b, bins[n+1])] = []
                for element in np.arange(b, bins[n+1], 0.01):
                    bin_elements['{0}-{1}'.format(b, bins[n+1])].append(round(element, 2))
        elif cat == 'hydrophobicity':
            for x in np.linspace(-4.23, 3.34, 51):
                bins.append(round(x, 2))
            for n, b in enumerate(bins):
                if n == len(bins)-1:
                    continue
                bin_elements['{0}-{1}'.format(b, bins[n+1])] = []
                for element in np.arange(b, bins[n+1], 0.01):
                    bin_elements['{0}-{1}'.format(b, bins[n+1])].append(round(element, 2))

        fieldnames = [
            'Sequence length',
            '#Sequences (for peptide q-value <= 1%)',
            '#Sequences (for peptide q-value <= 1% and #PSMs >= 2)',
        ]
        csv_out_name = output_filename.replace('.csv', '_{0}.csv'.format(cat))
        with open(csv_out_name, 'w') as csv_out:
            csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
            csv_writer.writeheader()
            out_dict_list = []
            plot_bins = {}
            peptides = results_dict['all']['peptides']['safe']
            for peptide in peptides:
                if cat == 'length':
                    v = len(peptide)
                else:
                    try:
                        v = peptide_info_dict[peptide][cat]
                    except:
                        continue
                for b in bin_elements.keys():
                    if b not in plot_bins.keys():
                        plot_bins[b] = 0
                    if v in bin_elements[b]:
                        plot_bins[b] += 1

            x_axis = []
            y_axis = []
            for n, b in enumerate(bins):
                out_dict = {}
                if n == len(bins)-1:
                    continue
                if cat == 'length':
                    bin_name = b
                else:
                    bin_name = '{0}-{1}'.format(b, bins[n+1])
                out_dict['Sequence length'] = bin_name
                out_dict['#Sequences (for peptide q-value <= 1%)'] = plot_bins[bin_name]
                out_dict_list.append(out_dict)
                x_axis.append(bin_name)
                y_axis.append(plot_bins[bin_name])

            trace_bar = go.Bar(
                x = x_axis,
                y = y_axis,
                name = 'FDR<=1%',
                marker=dict(
                    color='#1f77b4',
                ),
             )
            trace_list.append(trace_bar)

            plot_bins = {}
            peptides = results_dict['all']['peptides']['safe_num_specs']
            for peptide in peptides:
                if cat == 'length':
                    v = len(peptide)
                else:
                    try:
                        v = round(peptide_info_dict[peptide][cat], 2)
                    except:
                        continue
                for b in bin_elements.keys():
                    if b not in plot_bins.keys():
                        plot_bins[b] = 0
                    if v in bin_elements[b]:
                        plot_bins[b] += 1

            x_axis = []
            y_axis = []
            for n, b in enumerate(bins):
                if n == len(bins)-1:
                    continue
                out_dict = out_dict_list[n]
                if cat == 'length':
                    bin_name = b
                else:
                    bin_name = '{0}-{1}'.format(b, bins[n+1])
                x_axis.append(bin_name)
                y_axis.append(plot_bins[bin_name])
                out_dict['#Sequences (for peptide q-value <= 1% and #PSMs >= 2)'] = plot_bins[bin_name]
                csv_writer.writerow(out_dict)

            trace_bar = go.Bar(
                x = x_axis,
                y = y_axis,
                name = 'FDR<=1%; >=2PSMs',
                marker=dict(
                    color='#20b7c6',
                ),
             )
            trace_list.append(trace_bar)


        layout = go.Layout(
            barmode='group',
            legend=dict(orientation= "v"),
            yaxis=dict(
                title='Peptide sequences',
                # rangemode='tozero',
                showline=True,
                linecolor='rgb(0, 0, 0)',
                titlefont=dict(
                    size=16,
                    color='rgb(0, 0, 0)'
                ),
                tickfont=dict(
                    size=14,
                    color='rgb(0, 0, 0)',
                ),
                ticks='outside',
                ticklen=2,
                # dtick=0.25,
                tickwidth=0.25,
                tickcolor='rgb(0, 0, 0)',
                # ticksuffix='k',
                # tickangle=-90,
                side='left',
                autorange=True,
                # range=[0,50]
            ),
            xaxis=dict(
                showgrid=False,
                showline=False,
                zeroline=True,
                ticks='outside',
                ticklen=2,
                tickwidth=0.25,
                dtick=10,
                title=y_axis_dict['name'][cat],
                titlefont=dict(
                    size=16,
                    color='rgb(0, 0, 0)'
                ),
                tickangle=-90,
                tickfont=dict(
                    size=14,
                    color='rgb(0, 0, 0)'
                )
            ),
            showlegend=True,
            autosize=False,
            width=700,
            height=400,
            margin=go.layout.Margin(
                l=65,
                r=50,
                b=150,
                t=10,
                # pad=4
            ),
        )

        fig = go.Figure(data=trace_list, layout=layout)
        plotly.io.write_image(
            fig,
            output_filename.replace('.csv', '_{0}.pdf'.format(cat))
        )
        plotly.io.write_html(
            fig,
            output_filename.replace('.csv', '_{0}.html'.format(cat))
        )

def plot_psm_ident_rate(results_dict, output_filename='test.pdf'):
    '''
    Plot a barplot with the number of PSMs identified for each dataset
    and an overlaying scatter plot with the PSM assignment rate
    '''
    trace_list = []
    x_bar = []
    y_bar = []
    x_scatter = []
    y_scatter = []
    fieldnames = [
        'PRIDE ID',
        '#PSMs (in Million)',
        '#Spectra (in Million)',
        'PSM assignment rate',
    ]
    csv_out_name = output_filename.replace('.pdf', '.csv')
    with open(csv_out_name, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for PRIDE_ID in PRIDE_order:
            out_dict = {}
            out_dict['PRIDE ID'] = PRIDE_ID
            if PRIDE_ID == 'all':
                y_name = 'All datasets'
            else:
                y_name = PRIDE_ID
            print(PRIDE_ID)
            x_bar.append(y_name)
            y_bar.append(len(results_dict[PRIDE_ID]['spectra']['all'])/1000000)
            x_scatter.append(y_name)
            y_scatter.append(
                len(results_dict[PRIDE_ID]['spectra']['all'])/results_dict[PRIDE_ID]['num_spectra']
            )
            out_dict['#PSMs (in Million)'] = len(results_dict[PRIDE_ID]['spectra']['all'])/1000000
            out_dict['#Spectra (in Million)'] = results_dict[PRIDE_ID]['num_spectra']/1000000
            out_dict['PSM assignment rate'] = len(results_dict[PRIDE_ID]['spectra']['all'])/results_dict[PRIDE_ID]['num_spectra']
            csv_writer.writerow(out_dict)
    trace_bar = go.Bar(
        x = x_bar,
        y = y_bar,
        name = 'PSMs in Mio',
        # textposition = 'outside',
        # textangle = -90,
        # outsidetextfont=dict(
        #     size=8,
        #     color='rgb(0,0,0)'
        # ),
        marker=dict(
            color='#ff7f0e',
            # line=dict(
            #     color='rgb(8,48,107)',
            #     width=1.5,
            # )
        )
     )
    trace_list.append(trace_bar)

    trace_scatter = go.Scatter(
        x = x_scatter,
        y = y_scatter,
        name = 'MS2 identification rate',
        mode='markers',
        marker=dict(
            color = 'rgb(0,0,0)',
            symbol='cross',
            size=5,
        ),
        yaxis='y2'
    )
    trace_list.append(trace_scatter) 

    layout = go.Layout(
        legend=dict(orientation= "v"),
        yaxis=dict(
            title='Number of PSMs',
            # rangemode='tozero',
            showline=True,
            linecolor='#ff7f0e',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='#ff7f0e',
            ),
            ticks='outside',
            ticklen=2,
            dtick=0.5,
            tickwidth=0.25,
            tickcolor='#ff7f0e',
            # tickangle=-90,
            ticksuffix='M',
            side='left',
            autorange=False,
            range=[0,4]
        ),
        yaxis2=dict(
            title='PSM assignment rate',
            # rangemode='tozero',
            showgrid=False,
            showline=True,
            linecolor='rgb(0,0,0)',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(0,0,0)'
            ),
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            tickcolor='rgb(0,0,0)',
            overlaying='y',
            side='right'
        ),
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=True,
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            title='PRIDE ID',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickangle=-90,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        showlegend=False,
        autosize=False,
        width=500,
        height=400,
        margin=go.layout.Margin(
            l=65,
            r=50,
            # b=200
            b=150,
            t=10,
            # pad=4
        ),
    )

    fig = go.Figure(data=trace_list, layout=layout)
    plotly.io.write_image(fig, output_filename)

def plot_peptide_ident_rate(results_dict, output_filename='test.pdf'):
    '''
    Plots a barplot of identified peptide sequences for each dataset.
    This is overlayed with a scatter plot of the peptides/spectra ratio.
    '''
    trace_list = []
    x_bar = []
    y_bar = []
    x_scatter = []
    y_scatter = []
    fieldnames = [
        'PRIDE ID',
        '#Peptide sequences',
        '#Spectra (in Million)',
        'Peptide assignment rate',
    ]
    csv_out_name = output_filename.replace('.pdf', '.csv')
    with open(csv_out_name, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for PRIDE_ID in PRIDE_order:
            out_dict = {}
            out_dict['PRIDE ID'] = PRIDE_ID
            if PRIDE_ID == 'all':
                y_name = 'All datasets'
            else:
                y_name = PRIDE_ID
            print(PRIDE_ID)
            x_bar.append(y_name)
            y_bar.append(len(results_dict[PRIDE_ID]['peptides']['safe_num_specs']))
            x_scatter.append(y_name)
            y_scatter.append(
                len(results_dict[PRIDE_ID]['peptides']['safe_num_specs'])/results_dict[PRIDE_ID]['num_spectra']
            )
            out_dict['#Peptide sequences'] = len(results_dict[PRIDE_ID]['peptides']['safe_num_specs'])
            out_dict['#Spectra (in Million)'] = results_dict[PRIDE_ID]['num_spectra']/1000000
            out_dict['Peptide assignment rate'] = len(
                results_dict[PRIDE_ID]['peptides']['safe_num_specs']
            )/results_dict[PRIDE_ID]['num_spectra']
            csv_writer.writerow(out_dict)
    trace_bar = go.Bar(
        x = x_bar,
        y = y_bar,
        name = 'Peptide sequences',
        # textposition = 'outside',
        # textangle = -90,
        # outsidetextfont=dict(
        #     size=8,
        #     color='rgb(0,0,0)'
        # ),
        marker=dict(
            color='#20b7c6',
            # line=dict(
            #     color='rgb(8,48,107)',
            #     width=1.5,
            # )
        )
     )
    trace_list.append(trace_bar)

    trace_scatter = go.Scatter(
        x = x_scatter,
        y = y_scatter,
        name = 'Peptide assignment rate',
        mode='markers',
        marker=dict(
            color = 'rgb(0,0,0)',
            symbol='cross',
            size=5,
        ),
        yaxis='y2'
    )
    trace_list.append(trace_scatter) 

    layout = go.Layout(
        legend=dict(orientation= "v"),
        yaxis=dict(
            title='Peptide sequences',
            # rangemode='tozero',
            showline=True,
            linecolor='#20b7c6',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='#20b7c6',
            ),
            ticks='outside',
            ticklen=2,
            # dtick=0.25,
            tickwidth=0.25,
            # ticksuffix='k',
            tickcolor='#20b7c6',
            # tickangle=-90,
            side='left',
            autorange=True,
            # range=[0,1.0]
        ),
        yaxis2=dict(
            title='Peptide assignment rate',
            # rangemode='tozero',
            showgrid=False,
            showline=True,
            linecolor='rgb(0,0,0)',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(0,0,0)'
            ),
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            tickcolor='rgb(0,0,0)',
            overlaying='y',
            side='right',
            range=[0,0.08]
        ),
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=True,
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            title='PRIDE ID',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickangle=-90,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        showlegend=False,
        autosize=False,
        width=500,
        height=400,
        margin=go.layout.Margin(
            l=65,
            r=55,
            # b=200
            b=150,
            t=10,
            # pad=4
        ),
    )

    fig = go.Figure(data=trace_list, layout=layout)
    plotly.io.write_image(fig, output_filename)

def write_prot2file(results_dict, output_filename='test.txt'):
    with open('nonfunctional_and_spurious.txt', 'w') as nf_n_spur:
        for prot in results_dict['all']['protein_dict'].keys():
            if len(prot.split('<|>')) > 1:
                for p in prot.split('<|>'):
                    if p not in results_dict['all']['protein_groups']['safe_seq_num_spec_0005']:
                            continue
                    else:
                        prot = p
            else:
                if prot not in results_dict['all']['proteins']['safe_seq_num_spec_0005']:
                    continue
            if 'nonfunctional' in prot or 'spurious' in prot:
                print(prot, file=nf_n_spur)
                for seq in results_dict['all']['protein_dict'][prot].keys():
                    if seq in ['datasets', 'prot_q_value_psm', 'prot_q_value_seq', 'samples']:
                        continue
                    print('   {0}'.format(seq), file=nf_n_spur)
                    for spec in sorted(results_dict['all']['protein_dict'][prot][seq]['spec_title']):
                        print('       {0}'.format(spec), file=nf_n_spur)
                print('',file=nf_n_spur)
                print('',file=nf_n_spur)

    with open('trpA_pyrE2_proteins.txt', 'w') as arta:
        for prot in results_dict['all']['protein_dict'].keys():
            if len(prot.split('<|>')) > 1:
                if prot not in results_dict['all']['protein_groups']['safe_seq_num_spec_0005']:
                    continue
            else:
                if prot not in results_dict['all']['proteins']['safe_seq_num_spec_0005']:
                    continue
            for artA_substrate in [
                'HVO_0333',
                'HVO_0789',

                # 'HVO_0225',
                # 'HVO_0322',
                # 'HVO_0595',
                # 'HVO_1634',

                # 'HVO_A0174',
                # 'HVO_2603',
                # 'HVO_0437',

                # 'HVO_1225',
                # 'HVO_1222',
                # 'HVO_1224',
                # 'HVO_1223'

                # 'HVO_1203',
                # 'HVO_1200',
                # 'HVO_1210',

                # 'HVO_0620',
                # 'HVO_1160',
                # 'HVO_2451',
                # 'HVO_A0632'

                # 'HVO_0405',
                # 'HVO_1095',
                # 'HVO_1110',
                # 'HVO_2006',
                # 'HVO_2072',
                # 'HVO_2160',
                # 'HVO_2533',
                # 'HVO_A0263',
                # 'HVO_B0206',
            ]:
                if artA_substrate in prot:
                    print(prot, file=arta)
                    for seq in sorted(results_dict['all']['protein_dict'][prot].keys()):
                        if seq in ['prot_q_value_psm', 'prot_q_value_seq']:
                            continue
                        print('   {0}'.format(seq), file=arta)
                        if seq == 'datasets':
                            for dataset in sorted(results_dict['all']['protein_dict'][prot]['datasets']):
                                print('       {0}'.format(dataset), file=arta)
                            continue
                        # try:
                        #     for spec in results_dict['all']['protein_dict'][prot][seq]['spec_title']:
                        #         print('       {0}'.format(spec), file=arta)
                        # except:
                        #     print(seq)
                        #     exit()
                    print('',file=arta)
                    print('',file=arta)

def plot_fdr_results(fdr_dict, results_dict=None, output_filename='test.pdf', level='peptides', inference_dict=None):
    '''
    Plot scatterplots for the FDR on peptide (level='peptides') or protein (level='proteins') level.
    For peptides, FDRs are calculated for groups of peptide sequences with the same length.
    Proteins are sorted by are sorted by their q-value
    '''
    if level == 'peptides':
        datadict_list = []
        title_x = 'Peptide length'
        fieldnames = [
            'Peptide length',
            '#Targets',
            '#Decoys',
            'FDR (in %)'
        ]
        csv_out_name = output_filename.replace('.pdf', '.csv')
        with open(csv_out_name, 'w') as csv_out:
            csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
            csv_writer.writeheader()
            for cat in ['peptides_psm_level', 'peptides_seq_level', 'peptides_seq_level_2specs']:
                out_dict = {}
                out_dict['Peptide length'] = cat
                csv_writer.writerow(out_dict)
                print(cat)
                if cat == 'peptides_seq_level':
                    name = 'Peptide 1% FDR'
                if cat == 'peptides_psm_level':
                    name = 'PSM 1% FDR'
                if cat == 'peptides_seq_level_2psecs':
                    name = 'Peptide 1% FDR, 2 PSMs'
                datadict = {
                    'x': [],
                    'y': [],
                    'text': [],
                    'name': name,
                    'mode': 'markers',
                }
                for seq_length in sorted(fdr_dict[cat].keys()):
                    out_dict['Peptide length'] = seq_length
                    targets = 0
                    decoys = 0
                    for seq, seq_tuple in fdr_dict[cat][seq_length].items():
                        pep, is_decoy = seq_tuple
                        if is_decoy == 'TRUE' or is_decoy == 'true' or is_decoy is True:
                            decoy = True
                        elif is_decoy == 'FALSE' or is_decoy == 'false' or is_decoy is False:
                            decoy = False
                        if decoy:
                            decoys += 1
                        else:
                            targets += 1
                    if targets == 0 and decoys == 0:
                        continue
                    fdr = 100*decoys/(targets+decoys)
                    datadict['x'].append(seq_length)
                    datadict['y'].append(fdr)
                    print(seq_length, targets, decoys)
                    out_dict['#Targets'] = targets
                    out_dict['#Decoys'] = decoys
                    out_dict['FDR (in %)'] = fdr
                    csv_writer.writerow(out_dict)
                datadict_list.append(datadict)
        datadict_list.append({
            'x': [x for x in range(0,max(datadict['x'])+1)],
            'y': [1 for y in range(0,max(datadict['x'])+1)],
            'text': [],
            'name': "FDR 1%",
            'mode': 'lines',
        })

        plot_scatterplot(
            datadict_list,
            output_filename=output_filename,
            titles={
                'x': title_x,
                'y': 'FDR in %'
            },
            showlegend=True,
        )

    if level == 'proteins':
        prot2cov = plot_prot_coverage(results_dict, output_filename='tmp.pdf', inference_dict=inference_dict)
        datadict_list = []
        fieldnames = [
            '#Accepted proteins',
            '#Targets',
            '#Decoys',
            'FDR (in %)'
        ]
        csv_out_name = output_filename.replace('.pdf', '.csv')
        with open(csv_out_name, 'w') as csv_out:
            csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
            csv_writer.writeheader()
            for plot_level in ['proteins_psm_level','proteins_seq_level', 'proteins_seq_level_2specs']:
                out_dict = {}
                out_dict['#Accepted proteins'] = plot_level
                csv_writer.writerow(out_dict)
                print('plotting', plot_level)
                if plot_level == 'proteins_seq_level':
                    name = 'Sum(log<sub>10</sub>(Bayes PEP))\nfor sequences'
                elif plot_level == 'proteins_psm_level':
                    name = 'Sum(log<sub>10</sub>(Bayes PEP))\nfor PSMs'
                else:
                    name = 'Sum(log<sub>10</sub>(Bayes PEP))\nfor sequences, 2 PSMs'
                # new_output_filename = output_filename.replace('.pdf', '{0}.pdf'.format(plot_level))
                title_x = 'Accepted proteins'
                datadict = {
                    'x': [],
                    'y': [],
                    'text': [],
                    'name': name,
                    'mode': 'markers',
                }
                n = 0
                targets = 0
                decoys = 0
                safe = 0
                discarded = 0
                fdr_value_top2bottom = []
                td_top_2_bottom = []
                for prot, prot_tuple in sorted(fdr_dict[plot_level].items(), key=operator.itemgetter(1), reverse=False):
                    n += 1
                    score, is_decoy = prot_tuple
                    # print(score)
                    if is_decoy == 'TRUE' or is_decoy == 'true' or is_decoy is True:
                        decoy = True
                    elif is_decoy == 'FALSE' or is_decoy == 'false' or is_decoy is False:
                        decoy = False
                    if decoy:
                        decoys += 1
                    elif prot not in prot2cov.keys():
                        continue
                    else:
                        targets += 1
                    fdr = 100*decoys/(targets+decoys)
                    td_top_2_bottom.append(decoy)
                    fdr_value_top2bottom.append(fdr)

                min_obs_so_far = 100
                fdr_value_bottom2top = []
                for fdr in reversed(fdr_value_top2bottom):
                    if fdr < min_obs_so_far:
                        min_obs_so_far = fdr
                    fdr_value_bottom2top.append(min_obs_so_far)

                targets = 0
                decoys = 0
                for k, final_fdr in enumerate(reversed(fdr_value_bottom2top)):
                    current_bool = td_top_2_bottom[k]
                    if current_bool:
                        decoys += 1
                    else:
                        targets += 1
                    datadict['x'].append(k+1)
                    datadict['y'].append(final_fdr)
                    if final_fdr >= 0.5 and current_bool is False:
                        discarded += 1
                    elif final_fdr <= 0.5 and current_bool is False:
                        safe += 1
                    out_dict['#Accepted proteins'] = k+1
                    out_dict['#Targets'] = targets
                    out_dict['#Decoys'] = decoys
                    out_dict['FDR (in %)'] = final_fdr
                    csv_writer.writerow(out_dict)
                print('safe', safe)
                print('discarded', discarded)
                datadict_list.append(datadict)
                # datadict_list.append(datadict_decoys)

        range_max = max(datadict['x'])
        datadict_list.append({
            'x': [x for x in range(0,range_max+1)],
            'y': [1 for y in range(0,range_max+1)],
            'text': [],
            'name': "FDR 1%",
            'mode': 'lines',
        })

        datadict_list.append({
            'x': [x for x in range(0,range_max+1)],
            'y': [0.5 for y in range(0,range_max+1)],
            'text': [],
            'name': "FDR 0.5%",
            'mode': 'lines',
        })

        plot_scatterplot(
            datadict_list,
            output_filename=output_filename,
            titles={
                'x': title_x,
                'y': 'FDR in %'
            },
            showlegend=False,
        )

def plot_comparison_org(results_dict, output_filename='barplot_test.pdf', safe_idents=False):
    '''
    Plot a barplot comparing the results of the original analysis with results from the ArcPP analysis.
    For each data set, the bar represents the change in identifications (on PSM and peptide level) in %
    with 100% corresponding to the number of identifications in the original analysis.

    The key word argument "safe_idents" (True or False) can be used to only accept peptides
    with a q-value <= 0.01.
    '''
    datadict_list=[]
    color_dict = {
            'PSMs': '#ff7f0e',
            'Peptides': '#1f77b4',
        }
    fieldnames = [
        'PRIDE ID',
        'Original #PSMs',
        'ArcPP #PSMs',
        'Change #PSMs',
        'Change #PSMs (in %)',
        'Original #Peptides',
        'ArcPP #Peptides',
        'Change #Peptides',
        'Change #Peptides (in %)',
    ]
    csv_out_name = output_filename.replace('.pdf', '.csv')
    with open(csv_out_name, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for PRIDE_ID in PRIDE_order:
            out_dict = {}
            out_dict['PRIDE ID'] = PRIDE_ID
            for level in ['spectra','peptides',]:
                safe_word = 'all'
                if level == 'spectra':
                    name = 'PSMs'
                else:
                    name = 'Peptides'
                    if safe_idents:
                        safe_word = 'safe'
                datadict = {
                    'x': [],
                    'y': [],
                    'text': [],
                    'name': name,
                }
                if len(results_dict[PRIDE_ID]['original_results']['spectra']['all']) == 0:
                    continue
                if PRIDE_ID == 'all':
                    datadict['x'].append('All datasets')
                else:
                    datadict['x'].append(PRIDE_ID)
                new_idents = results_dict[PRIDE_ID][level][safe_word]
                out_dict['ArcPP #{0}'.format(name)] = len(new_idents)
                org_idents = results_dict[PRIDE_ID]['original_results'][level][safe_word]
                out_dict['Original #{0}'.format(name)] = len(org_idents)
                datadict['y'].append(
                    (100/len(org_idents) * len(new_idents)) - 100
                )
                out_dict['Change #{0}'.format(name)] = len(new_idents) - len(org_idents)
                out_dict['Change #{0} (in %)'.format(name)] = (100/len(org_idents) * len(new_idents)) - 100
                # datadict['text'].append(
                #     len(new_idents) - len(org_idents)
                # )
            csv_writer.writerow(out_dict)
            datadict_list.append(datadict)

    layout = go.Layout(
        barmode='group',
        # bargap=0.5,
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=True,
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            title='PRIDE ID',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickangle=-90,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        yaxis=dict(
            showgrid=True,
            showline=True,
            ticks='outside',
            ticklen=3,
            tickwidth=0.25,
            dtick=20,
            title='Change in identifications in %',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            ),
            # tickangle=-90,
            autorange = False,
            range = [-50, 180],
        ),
        showlegend=False,
        autosize=False,
        width=700,
        height=400,
        margin=go.layout.Margin(
            l=55,
            r=50,
            b=150,
            t=10,
        ),
    )
    plot_barplot(
        datadict_list=datadict_list,
        output_filename=output_filename,
        layout=layout,
        color_dict=color_dict
    )

def plot_comparison_detailed(results_dict, output_filename='barplot_test.pdf'):
    '''
    Plot a barplot comparing the results of the original analysis with results from the ArcPP analysis.
    This is a more detailed comparison than "plot_comparison" taking into account not only the final
    optimized results of the ArcPP reanalysis but also the intermediate step of merging results from
    3 search engines.
    For each data set, the bar represents the number of identifications (on PSM and peptide level) in %
    with 100% corresponding to the number of identifications in the original analysis.
    '''
    datadict_list=[]
    color_dict = {
        'peptides original_results': 'rgb(208,209,230)',
        'spectra original_results': 'rgb(253,208,162)',
        'peptides 3engines_results': 'rgb(116,169,207)',
        'spectra 3engines_results': 'rgb(253,141,60)',
        'peptides combined_results': 'rgb(43,140,190)',
        'spectra combined_results': 'rgb(230,85,13)',
        'peptides optimized_results': 'rgb(4,90,141)',
        'spectra optimized_results': 'rgb(166,54,3)',
    }
    name_dict = {
        'peptides original_results': 'original',
        'spectra original_results': 'original',
        'peptides 3engines_results': '3 engines',
        'spectra 3engines_results': '3 engines',
        'peptides combined_results': 'combined',
        'spectra combined_results': 'combined',
        'peptides optimized_results': 'optimized',
        'spectra optimized_results': 'optimized',
    }
    fieldnames = [
        'PRIDE ID',
        'Original #PSMs',
        'Original #PSMs (in %)',
        'ArcPP #PSMs 3engines',
        'ArcPP #PSMs 3engines (in %)',
        'ArcPP #PSMs optimized',
        'ArcPP #PSMs optimized (in %)',
        'Original #Peptides',
        'Original #Peptides (in %)',
        'ArcPP #Peptides 3engines',
        'ArcPP #Peptides 3engines (in %)',
        'ArcPP #Peptides optimized',
        'ArcPP #Peptides optimized (in %)',
    ]
    csv_out_name = output_filename.replace('.pdf', '.csv')
    with open(csv_out_name, 'w') as csv_out:
        csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
        csv_writer.writeheader()
        for level in ['spectra', 'peptides']:
            out_dict = {}
            if level == 'spectra':
                level_name = 'PSMs'
            else:
                level_name = 'Peptides'
            for cat in ['original_results', '3engines_results', 'optimized_results']:#'combined_results'
                name = level + ' ' + cat
                datadict = {
                    'x': [],
                    'y': [],
                    'text': [],
                    'name': name,
                }
                for PRIDE_ID in ['PXD006877', 'PXD011218']:
                    out_dict['PRIDE ID'] = PRIDE_ID
                    datadict['x'].append(PRIDE_ID)
                    org_idents = results_dict[PRIDE_ID]['original_results'][level]['all']
                    out_dict['Original #{0}'.format(level_name)] = len(org_idents)
                    if cat == 'optimized_results':
                        new_idents = results_dict[PRIDE_ID][level]['all']
                    else:
                        new_idents = results_dict[PRIDE_ID][cat][level]['all']
                    datadict['y'].append(
                        (100/len(org_idents) * len(new_idents))
                    )
                    # datadict['text'].append(
                    #     len(new_idents) - len(org_idents)
                    # )
                    if cat == 'original_results':
                        out_dict['Original #{0} (in %)'.format(level_name)] = 100/len(org_idents) * len(new_idents)
                    else:
                        cat_name = cat.split('_')[0]
                        out_dict['ArcPP #{0} {1}'.format(level_name, cat_name)] = len(new_idents)
                        out_dict['ArcPP #{0} {1} (in %)'.format(level_name, cat_name)] = 100/len(org_idents) * len(new_idents)
                    csv_writer.writerow(out_dict)
                datadict_list.append(datadict)

    layout = go.Layout(
        barmode='group',
        # bargap=0.5,
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=True,
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            title='PRIDE ID',
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            # tickangle=-90,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        yaxis=dict(
            showgrid=True,
            showline=True,
            ticks='outside',
            ticklen=3,
            tickwidth=0.25,
            dtick=20,
            title="Identifications in %",
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            ),
            # tickangle=-90,
            autorange = False,
            range = [0, 170],
        ),
        # legend=dict(
        #     x=1.0,
        #     y=1.0,
        #     bgcolor='rgba(255, 255, 255, 0)',
        #     bordercolor='rgba(255, 255, 255, 0)'
        # ),
        showlegend=False,
        autosize=False,
        width=400,
        height=300,
        margin=go.layout.Margin(
            l=55,
            r=50,
            b=150,
            t=10,
            # pad=4
        ),
    )
    plot_barplot(
        datadict_list=datadict_list,
        output_filename=output_filename,
        layout=layout,
        color_dict=color_dict
    )

def plot_spec_number(results_dict, color_by=None, output_filename='barplot_test.pdf'):
    '''
    Plots the total number of spectra for each dataset.
    The color_by keyword can be used to plot e.g. different instruments or labs
    in different colors.
    '''
    datadict_list = []
    if color_by is not None:
        for instrument in sorted(results_dict['all'][color_by]):
            datadict = {
                'x': [],
                'y': [],
                'text': [],
                'name': instrument,
            }
            for PRIDE_ID in sorted(results_dict.keys()):
                if PRIDE_ID == 'all':
                    continue
                datadict['x'].append(PRIDE_ID)
                if instrument == results_dict[PRIDE_ID]['instrument']:
                    datadict['y'].append(results_dict[PRIDE_ID]['num_spectra']/1000000)
                else:
                    datadict['y'].append(0)
            datadict_list.append(datadict)
    else:
        fieldnames = [
            'PRIDE ID',
            'Number of spectra',
        ]
        datadict = {
            'x': [],
            'y': [],
            'text': [],
            'name': 'spectra',
        }
        color_dict = {
            'spectra': 'rgb(125,125,125)'
        }
        csv_out_name = output_filename.replace('.pdf', '.csv')
        with open(csv_out_name, 'w') as csv_out:
            csv_writer = csv.DictWriter(csv_out, fieldnames=fieldnames, **csv_kwargs)
            csv_writer.writeheader()
            for PRIDE_ID in PRIDE_order:
                if PRIDE_ID == 'all':
                    continue
                print(PRIDE_ID, results_dict[PRIDE_ID]['num_spectra']/1000000)
                datadict['x'].append(PRIDE_ID)
                datadict['y'].append(results_dict[PRIDE_ID]['num_spectra']/1000000)
                out_dict = {
                    'PRIDE ID': PRIDE_ID,
                    'Number of spectra': results_dict[PRIDE_ID]['num_spectra']
                }
                csv_writer.writerow(out_dict)
            datadict_list.append(datadict)
    
    layout = go.Layout(
        barmode='group',
        bargap=0.5,
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=True,
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            title='PRIDE ID',
            titlefont=dict(
                size=7,
                color='rgb(0, 0, 0)'
            ),
            tickangle=-90,
            tickfont=dict(
                size=6,
                color='rgb(0, 0, 0)'
            )
        ),
        yaxis=dict(
            showgrid=True,
            showline=True,
            ticks='outside',
            ticklen=3,
            tickwidth=0.25,
            dtick=1,
            title='Number of spectra',
            titlefont=dict(
                size=7,
                color='rgb(0, 0, 0)'
            ),
            ticksuffix='M',
            tickfont=dict(
                size=6,
                color='rgb(0, 0, 0)'
            ),
            tickangle=-90,
            autorange = True,
            # range = [0, 200],
        ),
        # legend=dict(
        #     x=1.0,
        #     y=1.0,
        #     bgcolor='rgba(255, 255, 255, 0)',
        #     bordercolor='rgba(255, 255, 255, 0)'
        # ),
        showlegend=False,
        autosize=False,
        width=450,
        height=275,
        margin=go.layout.Margin(
            l=55,
            r=50,
            b=150,
            t=10,
            # pad=4
        ),
    )
    plot_barplot(
        datadict_list=datadict_list,
        output_filename=output_filename,
        layout=layout,
        color_dict=color_dict
    )


def plot_barplot(
    datadict_list=[],
    output_filename='test.png',
    layout=None,
    color_dict=None,
    orientation='v',
):
    '''
    Plots a barplot using Plotly.
    Given is a list of dictionaries with
    {
        'x' = [x_value_1, ...]
        'y' = [y_value_1, ...]
        'text' = [text_1, ...]
        'name' = 'name of the trace'
    }
    Each dictionary in the list is printed as a trace in the barplot.
    '''

    trace_list = []
    for datadict in datadict_list:
        trace =  go.Bar(
            x=datadict['x'],
            y=datadict['y'],
            text=datadict['text'],
            name=datadict['name'],
            textposition = 'outside',
            textangle = -90,
            outsidetextfont=dict(
                size=8,
                color='rgb(0,0,0)'
            ),
            marker=dict(
                color=color_dict[datadict['name']],
                # line=dict(
                #     color='rgb(8,48,107)',
                #     width=1.5,
                # )
            ),
            orientation=orientation,
            # opacity=0.6
        )
        trace_list.append(trace)

    # layout = go.Layout(
    #     barmode=barmode,
    #     # bargap=0.5,
    #     xaxis=dict(
    #         showgrid=False,
    #         showline=False,
    #         zeroline=True,
    #         ticks='outside',
    #         ticklen=2,
    #         tickwidth=0.25,
    #         title='PRIDE ID',
    #         titlefont=dict(
    #             size=16,
    #             color='rgb(0, 0, 0)'
    #         ),
    #         tickangle=-90,
    #         tickfont=dict(
    #             size=14,
    #             color='rgb(0, 0, 0)'
    #         )
    #     ),
    #     yaxis=dict(
    #         showgrid=True,
    #         showline=True,
    #         ticks='outside',
    #         ticklen=3,
    #         tickwidth=0.25,
    #         dtick=20,
    #         title=titles,
    #         titlefont=dict(
    #             size=16,
    #             color='rgb(0, 0, 0)'
    #         ),
    #         tickfont=dict(
    #             size=14,
    #             color='rgb(0, 0, 0)'
    #         ),
    #         # tickangle=-90,
    #         autorange = True,
    #         # range = [0, 200],
    #     ),
    #     legend=dict(
    #         x=1.0,
    #         y=1.0,
    #         bgcolor='rgba(255, 255, 255, 0)',
    #         bordercolor='rgba(255, 255, 255, 0)'
    #     ),
    #     showlegend=showlegend,
    #     autosize=False,
    #     width=700,
    #     height=400,
    #     margin=go.layout.Margin(
    #         l=55,
    #         r=50,
    #         b=150,
    #         t=10,
    #         # pad=4
    #     ),
    #     # bargap=0.15,
    #     # bargroupgap=0.1
    # )

    print('saving', output_filename)
    fig = go.Figure(data=trace_list, layout=layout)
    plotly.io.write_image(fig, output_filename)
    plotly.io.write_html(
        fig,
        output_filename.replace('.pdf', '.html')
    )
    print('Done')

def plot_scatterplot(datadict_list, output_filename='test.pdf', titles={}, showlegend=False):
    '''
        Plots a scatterplot using Plotly.
        Given is a list of dictionaries with
        {
            'x' = [x_value_1, ...]
            'y' = [y_value_1, ...]
            'text' = [text_1, ...]
            'name' = 'name of the trace'
        }
        Each dictionary in the list is printed as a trace in the barplot.
    '''
    marker_colors = {
        0: '#ff7f0e',
        1: '#1f77b4',
        2: '#20b7c6',
        3: 'rgb(111, 111, 111)',
        4: 'rgb(111, 111, 111)',
    }
    marker_shapes = {
        0: 'square',
        1: 'circle',
        2: 'triangle-up',
        3: 'square',
        4: 'square',
    }
    line_shapes = {
        0: 'dash',
        1: 'dash',
        2: 'dash',
        3: 'dot',
        4: 'dot',
    }
    trace_list = []
    for n, datadict in enumerate(datadict_list):
        trace = go.Scatter(
            x = datadict['x'],
            y = datadict['y'],
            text=datadict['text'],
            name=datadict['name'],
            mode=datadict['mode'],
            marker=dict(
                color = marker_colors[n],
                symbol=marker_shapes[n],
                size=4,
            ),
            line = dict(
                color = 'rgb(111, 111, 111)',
                width = 0.5,
                dash = line_shapes[n],
            ),
            opacity=1
        )
        trace_list.append(trace)

    layout = go.Layout(
        xaxis=dict(
            showgrid=False,
            showline=False,
            zeroline=True,
            autorange=False,
            # range=[0,3250],
            range=[0,55],
            ticks='outside',
            ticklen=2,
            tickwidth=0.25,
            # dtick=500,
            dtick=10,
            title=titles['x'],
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            # tickangle=-90,
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        yaxis=dict(
            showgrid=False,
            showline=True,
            autorange=False,
            range=[-0.1,5],
            ticks='outside',
            ticklen=3,
            tickwidth=0.25,
            title=titles['y'],
            titlefont=dict(
                size=16,
                color='rgb(0, 0, 0)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(0, 0, 0)'
            )
        ),
        legend=dict(
            x=1.0,
            y=1.0,
            bgcolor='rgba(255, 255, 255, 0)',
            bordercolor='rgba(255, 255, 255, 0)',
            orientation="h"
        ),
        showlegend=showlegend,
        autosize=False,
        width=350,
        height=350,
        margin=go.layout.Margin(
            l=50,
            r=50,
            b=150,
            t=10,
        ),
    )

    fig = go.Figure(data=trace_list, layout=layout)
    plotly.io.write_image(fig, output_filename)
    

def protein_inference(results_dict, safe_word='safe_seq_num_spec_0005'):
    '''
    A simple model for protein inference is used:
    If from a protein group, only one protein is securely
    identified in the same sample, this protein is assumed
    to be the identified one, instead of the group.
    '''
    inference_dict = {}
    for PRIDE_ID in results_dict.keys():
        if PRIDE_ID == 'all':
            continue
        proteins = results_dict[PRIDE_ID]['proteins'][safe_word]
        protein_groups = results_dict[PRIDE_ID]['protein_groups'][safe_word]
        for group in protein_groups:
            safe_prots_from_group = set()
            for prot in group.split('<|>'):
                if prot in proteins:
                    for sample in results_dict[PRIDE_ID]['protein_dict'][group]['samples']:
                        if sample in results_dict[PRIDE_ID]['protein_dict'][prot]['samples']:
                            safe_prots_from_group.add(prot)
            if group not in inference_dict.keys():
                inference_dict[group] = {'all':set()}
            if PRIDE_ID not in inference_dict[group].keys():
                inference_dict[group][PRIDE_ID] = set()
            inference_dict[group][PRIDE_ID] |= safe_prots_from_group
            inference_dict[group]['all'] |= safe_prots_from_group

    return inference_dict

if __name__ == '__main__':
    main(
        results_pkl=sys.argv[1],
        fdr_pkl=sys.argv[2],
        protein_pkl=sys.argv[3],
        peptide_pkl=sys.argv[4],
        plot_types=sys.argv[5:],
    )
