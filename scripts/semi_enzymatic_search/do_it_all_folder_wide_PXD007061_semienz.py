#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os

offsets = {
    '1-1.mzML': 0,
    '1-1N.mzML': 0,
    '1-3.mzML': 0,
    '1-3N.mzML': 0,
    '2-1.mzML': 0,
    '2-1N.mzML': 0,
    '2-3.mzML': 0,
    '2-3N.mzML': 0,
    '3-1.mzML': 0,
    '3-1N.mzML': 0,
    '3-3.mzML': 0,
    '3-3N.mzML': 0,
    '4-1.mzML': 0,
    '4-1N.mzML': 0,
    '4-3.mzML': 0,
    '4-3N.mzML': 0,
    '1-1_C-top20.mzML': -2,
    '1-1_M-top20.mzML': -2,
    '1-2_C-top20.mzML': 2,
    '1-2_M-top20.mzML': 0,
    '1-3_C-top20.mzML': 0,
    '1-3_M-top20.mzML': 0,
    '1-4_C-top20.mzML': 2,
    '1-4_M-top20.mzML': 0,
    '1-5_C-top20.mzML': 2,
    '1-5_M-top20.mzML': -2,
    '1-6_C-top20.mzML': -2,
    '1-6_M-top20.mzML': -4,
    '1-7_C-top20.mzML': 0,
    '1-7_M-top20.mzML': 2,
    '1-8_C-top20.mzML': 0,
    '1-8_M-top20.mzML': 0,
    '3-1_C-top20.mzML': -2,
    '3-1_M-top20.mzML': -2,
    '3-2_C-top20.mzML': -2,
    '3-2_M-top20.mzML': 0,
    '3-3_C-top20.mzML': 0,
    '3-3_M-top20.mzML': 0,
    '3-4_C-top20.mzML': 0,
    '3-4_M-top20.mzML': 0,
    '3-5_C-top20.mzML': 0,
    '3-5_M-top20.mzML': 0,
    '3-6_C-top20.mzML': -2,
    '3-6_M-top20.mzML': 0,
    '3-7_C-top20.mzML': -2,
    '3-7_M-top20.mzML': -2,
    '3-8_C-top20.mzML': -2,
    '3-8_M-top20.mzML': 2,
    '4-1_C-top20.mzML': 0,
    '4-1_M-top20.mzML': 0,
    '4-2_C-top20.mzML': 0,
    '4-2_M-top20.mzML': 0,
    '4-3_C-top20.mzML': 0,
    '4-3_M-top20.mzML': -2,
    '4-4_C-top20.mzML': 2,
    '4-4_M-top20.mzML': 2,
    '4-5_C-top20.mzML': 2,
    '4-5_M-top20.mzML': 2,
    '4-6_C-top20.mzML': 0,
    '4-6_M-top20.mzML': 2,
    '4-7_C-top20.mzML': 2,
    '4-7_M-top20.mzML': 0,
    '4-8_C-top20.mzML': 0,
    '4-8_M-top20.mzML': 0,
    'D_1-0_C-top20.mzML': -2,
    'D_1-0_M-top20.mzML': 0,
    'D_1-1_C-top20.mzML': 0,
    'D_1-1_M-top20.mzML': 0,
    'D_1-2_C-top20.mzML': 0,
    'D_1-2_M-top20.mzML': 0,
    'D_1-3_C-top20.mzML': 0,
    'D_1-3_M-top20.mzML': 0,
    'D_1-4_C-top20.mzML': 0,
    'D_1-4_M-top20.mzML': 0,
    'D_1-5_C-top20.mzML': 0,
    'D_1-5_M-top20.mzML': 0,
    'D_1-6_C-top20.mzML': 0,
    'D_1-6_M-top20.mzML': 0,
    'D_1-7_C-top20.mzML': 0,
    'D_1-7_M-top20.mzML': 0,
    'D_1-8_C-top20.mzML': 0,
    'D_1-8_M-top20.mzML': 0,
    'D_3-0_C-top20.mzML': -2,
    'D_3-0_M-top20.mzML': -2,
    'D_3-1_C-top20.mzML': 0,
    'D_3-1_M-top20.mzML': 0,
    'D_3-2_C-top20.mzML': 0,
    'D_3-2_M-top20.mzML': 0,
    'D_3-3_C-top20.mzML': 0,
    'D_3-3_M-top20.mzML': 0,
    'D_3-4_C-top20.mzML': 0,
    'D_3-4_M-top20.mzML': 0,
    'D_3-5_C-top20.mzML': 0,
    'D_3-5_M-top20.mzML': 0,
    'D_3-6_C-top20.mzML': 0,
    'D_3-6_M-top20.mzML': 0,
    'D_3-7_C-top20.mzML': 0,
    'D_3-7_M-top20.mzML': 0,
    'D_3-8_C-top20.mzML': 0,
    'D_3-8_M-top20.mzML': 0,
    'D_4-0_C-top20.mzML': 0,
    'D_4-0_M-top20.mzML': 2,
    'D_4-1_C-top20.mzML': 0,
    'D_4-1_M-top20.mzML': 2,
    'D_4-2_C-top20.mzML': 2,
    'D_4-2_M-top20.mzML': 0,
    'D_4-3_C-top20.mzML': 0,
    'D_4-3_M-top20.mzML': 0,
    'D_4-4_C-top20.mzML': 0,
    'D_4-4_M-top20.mzML': 0,
    'D_4-5_C-top20.mzML': 0,
    'D_4-5_M-top20.mzML': 0,
    'D_4-6_C-top20.mzML': 0,
    'D_4-6_M-top20.mzML': 0,
    'D_4-7_C-top20.mzML': -2,
    'D_4-7_M-top20.mzML': 0,
    'D_4-8_C-top20.mzML': -2,
    'D_4-8_M-top20.mzML': 2,
}

def main(folder=None, enzyme=None, target_decoy_database=None):
    '''
    '''
    # # define folder with mzML_files as sys.argv[1]
    mzML_files = []
    for mzml in glob.glob(os.path.join(folder, '*.mzML')):
        mzML_files.append(mzml)
    offset_files = []
    for sample in offsets.keys():
        offset_files.append(sample)
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
    for n, sample in enumerate(sorted(offsets.keys())):
        all_validated_result_files = []
        semi_validated_result_files = []
        full_validated_result_files = []
        combined_pep_result_files = []
        for search_engine in search_engines:
            spec_file = sample
            offset = offsets[spec_file]
            if offset == 'skip':
                continue
            uc.params['machine_offset_in_ppm'] = offset
            dirname = os.path.join(folder)
            mzml_file = os.path.join(dirname, spec_file)
            mgf_file = mzml_file.replace('.mzML', '.mgf')
            if os.path.exists(mgf_file) is False:
                continue
            mgf_file = uc.convert(
                input_file=mzml_file,
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

            all_validated_csv = uc.validate(
                input_file  = unified_search_results,
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
                input_file = unified_search_results,
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
                input_file = unified_search_results,
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

