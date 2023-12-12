#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os

offsets = {
    'ZJ_C1-8': {
        'Z_Jevtic_040316_C+Dfractions_01_C1-CTL_01_F1.mzML': 0,
        'Z_Jevtic_040316_C+Dfractions_01_C2-CTL_01_F2.mzML': -2,
        'Z_Jevtic_040316_C+Dfractions_01_C3-CTL_01_F3.mzML': -4,
        'Z_Jevtic_040316_C+Dfractions_01_C4-CTL_01_F4.mzML': -4,
        'Z_Jevtic_040316_C+Dfractions_01_C5-CTL_01_F5.mzML': -6,
        'Z_Jevtic_040316_C+Dfractions_01_C6-CTL_01_F6.mzML': 0,
        'Z_Jevtic_040316_C+Dfractions_01_C7-CTL_01_F7.mzML': 2,
        'Z_Jevtic_040316_C+Dfractions_01_C8-CTL_01_F8.mzML': -4,
    },
    'ZJ_D1-8': {
        'Z_Jevtic_040316_C+Dfractions_01_D1-DTL_01_F1.mzML': 0,
        'Z_Jevtic_040316_C+Dfractions_01_D2-DTL_01_F2.mzML': 0,
        'Z_Jevtic_040316_C+Dfractions_01_D3-DTL_01_F3.mzML': 4,
        'Z_Jevtic_040316_C+Dfractions_01_D4b-DTL_01_F4.mzML': 0,
        'Z_Jevtic_040316_C+Dfractions_01_D5-DTL_01_F5.mzML': 2,
        'Z_Jevtic_040316_C+Dfractions_01_D6-DTL_01_F6.mzML': 0,
        'Z_Jevtic_040316_C+Dfractions_01_D7-DTL_01_F7.mzML': -4,
        'Z_Jevtic_040316_C+Dfractions_01_D8-DTL_01_F8.mzML': 4,
    },
    'ZJ_HALVO_01': {
        'Z_Jevtic_071215_HALVO_F01_01-5ul_injection.mzML': 2,
        'Z_Jevtic_071215_HALVO_F02_01-5ul_injection.mzML': 2,
        'Z_Jevtic_071215_HALVO_F03_01-5ul_injection.mzML': 0,
        'Z_Jevtic_071215_HALVO_F04_01-5ul_injection.mzML': 0,
        'Z_Jevtic_071215_HALVO_F05_01-5ul_injection.mzML': -2,
        'Z_Jevtic_071215_HALVO_F06_01-5ul_injection.mzML': 0,
        'Z_Jevtic_071215_HALVO_F07_01-5ul_injection.mzML': -2,
        'Z_Jevtic_071215_HALVO_F08_01-5ul_injection.mzML': -2,
    },
    'ZJ_HALVO_02': {
        'Z_Jevtic_071215_HALVO_F01_02-5ul_injection.mzML': 2,
        'Z_Jevtic_071215_HALVO_F02_02-5ul_injection.mzML': 4,
        'Z_Jevtic_071215_HALVO_F03_02-5ul_injection.mzML': 0,
        'Z_Jevtic_071215_HALVO_F04_02-5ul_injection.mzML': 2,
        'Z_Jevtic_071215_HALVO_F05_02-5ul_injection.mzML': 0,
        'Z_Jevtic_071215_HALVO_F06_02-5ul_injection.mzML': -2,
        'Z_Jevtic_071215_HALVO_F07_03-5ul_injection.mzML': 0,
        'Z_Jevtic_071215_HALVO_F08_03-5ul_injection.mzML': 0,
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
        'precursor_mass_tolerance_minus' : 20,
        'precursor_mass_tolerance_plus': 20,
        'frag_mass_tolerance' : 20,
        'frag_mass_tolerance_unit': 'ppm',
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
        'pymzml_spec_id_attribute': {'index': None},
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
