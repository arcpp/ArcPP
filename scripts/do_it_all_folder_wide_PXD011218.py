#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os

offsets = {
    'MIG_Mem_exp_rep1A': {
        'MIG_1-1.mzML': -2,
        'MIG_1-2.mzML': 0,
        'MIG_1-3.mzML': 0,
        'MIG_1-4.mzML': 2,
    },
    'MIG_Mem_exp_rep2A': {
        'MIG_2-1.mzML': -2,
        'MIG_2-2.mzML': -2,
        'MIG_2-3.mzML': 2,
        'MIG_2-4.mzML': 2,
    },
    'MIG_Mem_exp_rep3A': {
        'MIG_3-1.mzML': 2,
        'MIG_3-2.mzML': -2,
        'MIG_3-3.mzML': 2,
        'MIG_3-4.mzML': 2,
    },
    'MIG_Mem_exp_rep4A': {
        'MIG_EXP1_botttom20.mzML': 0,
        'MIG_EXP1_top20.mzML': 0,
    },
    'MIG_Mem_exp_rep5A': {
        'MIG_EXP2_botttom20.mzML': 0,
        'MIG_EXP2_top20.mzML': 2,
    },
    'MIG_Mem_exp_rep6A': {
        'MIG_EXP3_botttom20.mzML': 2,
        'MIG_EXP3_top20.mzML': 2,
    },
    'MIG_Mem_exp_rep1B': {
        'MIG_L-1M.mzML': 0,
    },
    'MIG_Mem_exp_rep2B': {
        'MIG_L-2M.mzML': 0,
    },
    'MIG_Mem_exp_rep3B': {
        'Mig_L-3M.mzML': 2,
    },
    'MIG_Mem_stat_rep1B': {
        'MIG_S-1M.mzML': 0,
    },
    'MIG_Mem_stat_rep2B': {
        'MIG_S-2M.mzML': 0,
    },
    'MIG_Mem_stat_rep3B': {
        'Mig_S-3M.mzML': 0,
    },
    'MIG_Cyt_exp_rep1B': {
        'MIG_L-1C.mzML': 0,
    },
    'MIG_Cyt_exp_rep2B': {
        'MIG_L-2C.mzML': -4,
    },
    'MIG_Cyt_exp_rep3B': {
        'Mig_L-3C.mzML': 0,
    },
    'MIG_Cyt_stat_rep1B': {
        'MIG_S-1C.mzML': 0,
    },
    'MIG_Cyt_stat_rep2B': {
        'MIG_S-2C.mzML': 0,
    },
    'MIG_Cyt_stat_rep3B': {
        'Mig_S-3C.mzML': 0,
    },
    
    'MIG_Mem_stat_rep4A': {
        'MIG_ST1_botttom20.mzML': 0,
        'MIG_ST1_top20.mzML': 2,
    },
    'MIG_Mem_stat_rep5A': {
        'MIG_ST2_botttom20.mzML': 0,
        'MIG_ST2_top20.mzML': 0,
    },
    'MIG_Mem_stat_rep6A': {
        'MIG_ST3_botttom20.mzML': 2,
        'MIG_ST3_top20.mzML': 2,
    },
    'WT_Mem_exp_rep1A': {
        'WT_1-1.mzML': 0,
        'WT_1-2.mzML': 0,
        'WT_1-3.mzML': 0,
        'WT_1-4.mzML': 2,
    },
    'WT_Mem_exp_rep2A': {
        'WT_2-1.mzML': 0,
        'WT_2-2.mzML': 0,
        'WT_2-3.mzML': 2,
        'WT_2-4.mzML': 0,
    },
    'WT_Mem_exp_rep3A': {
        'WT_3-1.mzML': 2,
        'WT_3-2.mzML': 0,
        'WT_3-3.mzML': 2,
        'WT_3-4.mzML': 0,
    },
}

def main(folder=None, enzyme=None, target_decoy_database=None):
    '''
    '''
    # # define folder with mzML_files as sys.argv[1]
    mzML_files = []
    offset_files = []
    for sample in offsets.keys():
        for mzml in glob.glob(os.path.join(folder, sample, '*.mzML')):
            mzML_files.append(mzml)
        for offset_file in offsets[sample].keys():
            offset_files.append(offset_file)
    for mzml in mzML_files:
        if os.path.basename(mzml) not in offset_files:
            print('mzML file in folder but NOT in offset dict: {0}'.format(mzml))
            exit()

    mass_spectrometer = 'QExactive+'
    search_engines  = [
        'xtandem_vengeance',
        'msfragger_20190222',
        'msgfplus_v2019_04_18',
    ]

    validation_engine = 'percolator_3_4_0'

    params = {
        'database'      : target_decoy_database,
        'enzyme'        : enzyme,
        'precursor_mass_tolerance_minus' : 8,
        'precursor_mass_tolerance_plus': 8,
        'frag_mass_tolerance' : 0.4,
        'frag_mass_tolerance_unit': 'da',
        'rounded_mass_decimals':2,
        '-xmx'     : '32g',
        'peptide_mapper_class_version' : 'UPeptideMapper_v4',
        'use_pyqms_for_mz_calculation' : True,
        'percolator_post_processing': 'mix-max',
        'psm_defining_colnames': [
            'Spectrum Title',
            'Sequence',
            'Modifications',
            'Charge',
            'Is decoy',
        ],
    }

    uc = ursgal.UController(
        profile = mass_spectrometer,
        params = params
    )

    all_result_files = []
    for n, sample in enumerate(offsets.keys()):
        validated_result_files = []
        combined_pep_result_files = []
        for search_engine in search_engines:
            results = []
            for spec_file in offsets[sample].keys():
                basename = spec_file
                dirname = os.path.join(folder)
                offset = offsets[sample][basename]
                spec_file_path = os.path.join(dirname, basename)
                if offset == 'skip':
                    continue
                uc.params['machine_offset_in_ppm'] = offset
                mgf_file = uc.convert(
                    input_file=spec_file_path,
                    engine='mzml2mgf_2_0_0',
                )

                uc.params['modifications'] = [
                    'C,fix,any,Carbamidomethyl',
                    'M,opt,any,Oxidation',
                    '*,opt,Prot-N-term,Acetyl',
                ]
                search_result = uc.search_mgf(
                    input_file = mgf_file,
                    engine = search_engine,
                )

                converted_result = uc.convert(
                    input_file=search_result,
                    guess_engine = True,
                )
                
                mapped_results = uc.execute_misc_engine(
                    input_file=converted_result,
                    engine='upeptide_mapper',
                )

                unified_search_results = uc.execute_misc_engine(
                    input_file = mapped_results,
                    engine='unify_csv'
                )

                results.append(unified_search_results)

                # validated_single_csv = uc.validate(
                #     input_file  = unified_search_results,
                #     engine      = validation_engine,
                # )
                # 
                # uc.params['csv_filter_rules'] = [
                #     # ['Is decoy', 'equals', 'false'],
                #     ['combined PEP','lte', 0.01],
                #     ['Conflicting uparam', 'contains_not', 'enzyme'],
                # ]
                # filtered_combined_results = uc.execute_misc_engine(
                #     input_file = validated_single_csv,
                #     engine='filter_csv',
                # )
        
            uc.params['prefix'] = sample
            results_one_engine = uc.execute_misc_engine(
                input_file = results,
                engine='merge_csvs',
                # merge_duplicates=True,
            )
            uc.params['prefix'] = ''

            validated_csv = uc.validate(
                input_file  = results_one_engine,
                engine      = validation_engine,
            )
            # filtered_combined_results = uc.execute_misc_engine(
            #         input_file = validated_csv,
            #         engine='filter_csv',
            #     )

            validated_result_files.append(validated_csv)

        combined_results = uc.combine_search_results(
            input_files=validated_result_files,
            engine='combine_pep_1_0_0',
        )

        uc.params['csv_filter_rules'] = [
            ['combined PEP','lte', 0.01],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ]
        filtered_combined_results = uc.execute_misc_engine(
            input_file = combined_results,
            engine='filter_csv',
        )
        all_result_files.append(filtered_combined_results)
    
    results_all_files = uc.execute_misc_engine(
        input_file = all_result_files,
        engine='merge_csvs',
        merge_duplicates=True,
    )

    uc.params.update({
        'validation_score_field': 'combined PEP',
        'bigger_scores_better': False,
        'num_compared_psms': 10,
        'accept_conflicting_psms': False,
        'threshold_is_log10': True,
        'score_diff_threshold': 1,
        'psm_defining_colnames': [
            'Spectrum Title',
            'Sequence',
        ],
    })

    sanitized_combined_results = uc.execute_misc_engine(
        input_file = results_all_files,
        engine='sanitize_csv',
    )

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(main.__doc__)
        exit()
    main(
        folder                = sys.argv[1],
        enzyme                = sys.argv[2],
        target_decoy_database = sys.argv[3],
    )
