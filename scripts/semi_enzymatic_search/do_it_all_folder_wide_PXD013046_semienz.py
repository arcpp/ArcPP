#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os

offsets = {
    'ABI_Mem_exp_rep1': {
        '10_131129073758.mzML': 2,
        '10.mzML': 0,
    },
    'ABI_Mem_exp_rep2': {
        '10A_131129150241.mzML': 0,
        '10A.mzML': 2,
    },
    'ABI_Mem_exp_rep3': {
        '10B_131129222708.mzML': 4,
        '10B.mzML': 0,
    },
    'ABI_Mem_exp_rep4': {
        '10C_131130055123.mzML': 2,
        '10C.mzML': 4,
    },
    'ABI_Cyt_exp_rep1': {
        '9_131125033003.mzML': 0,
        '9.mzML': -4,
    },
    'ABI_Cyt_exp_rep2': {
        '9A_131125110110.mzML': 0,
        '9A.mzML': -2,
    },
    'ABI_Cyt_exp_rep3': {
        '9B_131125223541.mzML': -2,
        '9B.mzML': -2,
    },
    'ABI_Cyt_exp_rep4': {
        '9C_131126060708.mzML': 0,
        '9C.mzML': -2,
    },
    'ABI_Mem_stat_rep1': {
        '12_131130131543.mzML': 0,
        '12.mzML': -2,
    },
    'ABI_Mem_stat_rep2': {
        '12A_131130204006.mzML': 0,
        '12A.mzML': 0,
    },
    'ABI_Mem_stat_rep3': {
        '12B_131201040426.mzML': 0,
        '12B.mzML': 2,
    },
    'ABI_Mem_stat_rep4': {
        '12C_131201112849.mzML': 2,
        '12C.mzML': 0,
    },
    'ABI_Cyt_stat_rep1': {
        '11_131126133808.mzML': -2,
        '11.mzML': -4,
    },
    'ABI_Cyt_stat_rep2': {
        '11A_131126210917.mzML': -2,
        '11A.mzML': 0,
    },
    'ABI_Cyt_stat_rep3': {
        '11B_131127044005.mzML': -0,
        '11B.mzML': -2,
    },
    'ABI_Cyt_stat_rep4': {
        '11C_131127121052.mzML': -2,
        '11C.mzML': -4,
    },
    'HVLON3_Mem_exp_rep1': {
        '6_top20.mzML': 2,
        '6_botttom20.mzML': 2,
    },
    'HVLON3_Mem_exp_rep2': {
        '6A_botttom20.mzML': 2,
        '6A_top20.mzML': 0,
    },
    'HVLON3_Mem_exp_rep3': {
        '6B_131101122740.mzML': 0,
        '6B.mzML': 0,
    },
    'HVLON3_Mem_exp_rep4': {
        '6C_131128164906.mzML': 2,
        '6C.mzML': 0,
    },
    'HVLON3_Cyt_exp_rep1': {
        '5_131112222612.mzML': 0,
        '5.mzML': 2,
    },
    'HVLON3_Cyt_exp_rep2': {
        '5A_131113055702.mzML': 4,
        '5A.mzML': 2,
    },
    'HVLON3_Cyt_exp_rep3': {
        '5B_131113132754.mzML': 0,
        '5B.mzML': 0,
    },
    'HVLON3_Cyt_exp_rep4': {
        '5C_131113205840.mzML': 0,
        '5C.mzML': 0,
    },
    'HVLON3_Mem_stat_rep1': {
        '8.mzML': 2,
        '8_131101201658.mzML': 2,
    },
    'HVLON3_Mem_stat_rep2': {
        '8A_131102040617.mzML': 0,
        '8A.mzML': 2,
    },
    'HVLON3_Mem_stat_rep3': {
        '8B_131102115537.mzML': 0,
        '8B.mzML': 2,
    },
    'HVLON3_Mem_stat_rep4': {
        '8C_131129001334.mzML': 0,
        '8C.mzML': -2,
    },
    'HVLON3_Cyt_stat_rep1': {
        '7.mzML': 0,
        '7_131122192428.mzML': -2,
    },
    'HVLON3_Cyt_stat_rep2': {
        '7A_131123025552.mzML': -2,
        '7A.mzML': 0,
    },
    'HVLON3_Cyt_stat_rep32': {
        '7B_131125150408.mzML': -2,
        '7B.mzML': 0,
    },
    'HVLON3_Cyt_stat_rep4': {
        '7C_131124195907.mzML': -2,
        '7C.mzML': -2,
    },
    'WT_Mem_exp_rep1': {
        '2_top20.mzML': 2,
        '2_botttom20.mzML': 0,
    },
    'WT_Mem_exp_rep2': {
        '2A_botttom20.mzML': 0,
        '2A_top20_140615200858.mzML': 0,
    },
    'WT_Mem_exp_rep3': {
        '2B_131030002035.mzML': 2,
        '2B.mzML': 2,
    },
    'WT_Mem_exp_rep4': {
        '2C_131128020026.mzML': 2,
        '2C.mzML': 0,
    },
    'WT_Cyt_exp_rep1': {
        '1.mzML': 0,
        '1_131111162318.mzML': 4,
    },
    'WT_Cyt_exp_rep2': {
        '1A_131111235403.mzML': 2,
        '1A.mzML': 2,
    },
    'WT_Cyt_exp_rep3': {
        '1B_131112072445.mzML': -2,
        '1B.mzML': -2,
    },
    'WT_Cyt_exp_rep4': {
        '1C_131112145528.mzML': 0,
        '1C.mzML': 4,
    },
    'WT_Mem_stat_rep1': {
        '4_top20.mzML': 4,
        '4_botttom20.mzML': 0,
    },
    'WT_Mem_stat_rep2': {
        '4A_botttom20.mzML': 0,
        '4A_top20.mzML': 2,
    },
    'WT_Mem_stat_rep3': {
        '4B_131030234822.mzML': -8,
        '4B.mzML': 0,
    },
    'WT_Mem_stat_rep4': {
        '4C_131128092443.mzML': 2,
        '4C.mzML': 4,
    },
    'WT_Cyt_stat_rep1': {
        '3.mzML': -4,
        '3_131114042927.mzML': -2,
    },
    'WT_Cyt_stat_rep2': {
        '3A_131114120015.mzML': -4,
        '3A.mzML': 0,
    },
    'WT_Cyt_stat_rep3': {
        '3B_131114193105.mzML': -4,
        '3B.mzML': -4,
    },
    'WT_Cyt_stat_rep4': {
        '3C_131115030157.mzML': 0,
        '3C.mzML': -2,
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
            print('mzML file in folder but NOT in offset dict: {}'.format(mzml))
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
        'csv_filter_rules' : [
            ['Is decoy', 'equals', 'false'],
            ['PEP','lte', 0.01],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ],
        'precursor_mass_tolerance_minus' : 8,
        'precursor_mass_tolerance_plus': 8,
        'frag_mass_tolerance' : 0.4,
        'frag_mass_tolerance_unit': 'da',
        'rounded_mass_decimals':2,
        '-xmx'     : '32g',
        'peptide_mapper_class_version' : 'UPeptideMapper_v4',
        'use_pyqms_for_mz_calculation' : True,
        'semi_enzyme' : True,
        'precursor_min_charge': 1,
        'precursor_max_charge': 5,
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
    semi_result_files = []
    full_result_files = []
    for n, sample in enumerate(offsets.keys()):
        all_validated_result_files = []
        semi_validated_result_files = []
        full_validated_result_files = []
        combined_pep_result_files = []
        for search_engine in search_engines:
            results = []
            for spec_file in offsets[sample].keys():
                offset = offsets[sample][spec_file]
                if offset == 'skip':
                    continue
                uc.params['machine_offset_in_ppm'] = offset
                dirname = folder
                mzml_file = os.path.join(dirname, spec_file)
                mgf_file = uc.convert(
                    input_file=mzml_file,
                    engine='mzml2mgf_2_0_0',
                )
                uc.params['modifications'] = [
                    # 'C,fix,any,Carbamidomethyl',
                    'M,opt,any,Oxidation',
                    '*,opt,Prot-N-term,Acetyl',
                ]
                
                search_result = uc.search_mgf(
                    input_file = mgf_file,
                    engine = search_engine,
                )
                uc.params['prefix'] = ''

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
                    # force = True,
                )

            results_one_engine = uc.execute_misc_engine(
                input_file = results,
                engine='merge_csvs',
                # merge_duplicates=True,
                # force=True,
            )

            all_validated_csv = uc.validate(
                input_file  = results_one_engine,
                engine      = validation_engine,
            )
            all_validated_result_files.append(all_validated_csv)

            uc.params.update({
                'csv_filter_rules' : [
                    ['Enzyme Specificity', 'contains_not', 'full'],
                ],
                'prefix': 'Semi',
            })
            semi_filtered_csv = uc.execute_misc_engine(
                input_file = results_one_engine,
                engine='filter_csv',
            )
            semi_validated_csv = uc.validate(
                input_file  = semi_filtered_csv,
                engine      = validation_engine,
            )
            semi_validated_result_files.append(semi_validated_csv)

            uc.params.update({
                'csv_filter_rules' : [
                    ['Enzyme Specificity', 'contains', 'full'],
                ],
                'prefix': 'Full',
            })
            full_filtered_csv = uc.execute_misc_engine(
                input_file = results_one_engine,
                engine='filter_csv',
            )
            full_validated_csv = uc.validate(
                input_file  = full_filtered_csv,
                engine      = validation_engine,
            )
            full_validated_result_files.append(full_validated_csv)
            uc.params.update({
                'csv_filter_rules' : [
                    ['Is decoy', 'equals', 'false'],
                    ['PEP','lte', 0.01],
                    ['Conflicting uparam', 'contains_not', 'enzyme'],
                ],
                'prefix': '',
            })

        all_combined_results = uc.combine_search_results(
            input_files=all_validated_result_files,
            engine='combine_pep_1_0_0',
        )
        semi_combined_results = uc.combine_search_results(
            input_files=semi_validated_result_files,
            engine='combine_pep_1_0_0',
        )
        full_combined_results = uc.combine_search_results(
            input_files=full_validated_result_files,
            engine='combine_pep_1_0_0',
        )

        uc.params['csv_filter_rules'] = [
            # ['Is decoy', 'equals', 'false'],
            ['combined PEP','lte', 0.01],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ]
        all_filtered_combined_results = uc.execute_misc_engine(
            input_file = all_combined_results,
            engine='filter_csv',
        )
        all_result_files.append(all_filtered_combined_results)

        semi_filtered_combined_results = uc.execute_misc_engine(
            input_file = semi_combined_results,
            engine='filter_csv',
        )
        semi_result_files.append(semi_filtered_combined_results)

        full_filtered_combined_results = uc.execute_misc_engine(
            input_file = full_combined_results,
            engine='filter_csv',
        )
        full_result_files.append(full_filtered_combined_results)
    
    for l in [all_result_files, semi_result_files, full_result_files]:
        uc.params.update({
            'psm_defining_colnames': [
                'Spectrum Title',
                'Sequence',
                'Modifications',
                'Charge',
                'Is decoy',
            ],
        })
        all_files = uc.execute_misc_engine(
            input_file = l,
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
            input_file = all_files,
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
