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
    # define folder with mzML_files as sys.argv[1]
    mzML_files = []
    for mzml in glob.glob(os.path.join(folder, '*.mzML')):
        mzML_files.append(os.path.basename(mzml))
    offset_files = []
    for sample in offsets.keys():
        for spec_file in offsets[sample].keys():
            offset_files.append(spec_file)
    for mzml in mzML_files:
        if mzml not in offset_files:
            print('mzML file in folder but NOT in offset dict: {}'.format(mzml))
            import pprint
            pprint.pprint(offset_files)
            exit()

    mass_spectrometer = 'QExactive+'
    search_engines  = [
        'xtandem_vengeance',
        # 'msfragger_2_3',
        'msgfplus_v2019_07_03',
    ]

    validation_engine = 'percolator_3_4_0'

    params = {
        'database'      : target_decoy_database,
        'enzyme'        : enzyme,
        'precursor_mass_tolerance_minus' : 10,
        'precursor_mass_tolerance_plus': 10,
        'frag_mass_tolerance' : 10,
        'frag_mass_tolerance_unit': 'ppm',
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
        'max_missed_cleavages': 2,
        'pymzml_spec_id_attribute': {'index': None},
    }

    Hvo_Glyco = [
        '',
        'N,opt,any,Hex',
        'N,opt,any,1427',
        'N,opt,any,Hex(1)HexA(2),C18H26O17',
        'N,opt,any,Hex(1)HexA(3),C24H34O23',
        'N,opt,any,Hex(1)HexA(2)MeHexA(1)Hex(1),C31H46O28',
        'N,opt,any,Hex(1)HexA(2)MeHexA(1),C25H36O23',
        'N,opt,any,SO3Hex(1),C6H10O8S1',
        'N,opt,any,SO3Hex(1)Hex(1),C12H20O13S1',
        'N,opt,any,SO3Hex(1)Hex(2),C18H30O18S1',
        'N,opt,any,SO3Hex(1)Hex(2)dHex(1),C24H40O22S1',
    ]

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
                for n, mod in enumerate(Hvo_Glyco):
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

                    if n == 0:
                        uc.params['modifications'] = [
                            'C,fix,any,Carbamidomethyl',
                            'M,opt,any,Oxidation',
                            '*,opt,Prot-N-term,Acetyl',
                        ]
                    else:
                        uc.params['modifications'] = [
                            'C,fix,any,Carbamidomethyl',
                            'M,opt,any,Oxidation',
                            '*,opt,Prot-N-term,Acetyl',
                            'S,opt,any,Hex(2)',
                            'T,opt,any,Hex(2)',
                        ]
                        uc.params['modifications'].append(mod)
                        uc.params['prefix'] = mod.split(',')[3]

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
                    )

                    results.append(unified_search_results)

    #                 # validated_single_csv = uc.validate(
    #                 #     input_file  = unified_search_results,
    #                 #     engine      = validation_engine,
    #                 # )
    #                 # 
    #                 # uc.params['csv_filter_rules'] = [
    #                 #     # ['Is decoy', 'equals', 'false'],
    #                 #     ['combined PEP','lte', 0.01],
    #                 #     ['Conflicting uparam', 'contains_not', 'enzyme'],
    #                 # ]
    #                 # filtered_combined_results = uc.execute_misc_engine(
    #                 #     input_file = validated_single_csv,
    #                 #     engine='filter_csv',
    #                 # )
        
    #         uc.params['prefix'] = sample
    #         results_one_engine = uc.execute_misc_engine(
    #             input_file = results,
    #             engine='merge_csvs',
    #             # merge_duplicates=True,
    #         )
    #         uc.params['prefix'] = ''

    #         validated_csv = uc.validate(
    #             input_file  = results_one_engine,
    #             engine      = validation_engine,
    #         )
    #         # filtered_combined_results = uc.execute_misc_engine(
    #         #         input_file = validated_csv,
    #         #         engine='filter_csv',
    #         #     )

    #         validated_result_files.append(validated_csv)

    #     combined_results = uc.combine_search_results(
    #         input_files=validated_result_files,
    #         engine='combine_pep_1_0_0',
    #     )

    #     uc.params['csv_filter_rules'] = [
    #         ['combined PEP','lte', 0.01],
    #         ['Conflicting uparam', 'contains_not', 'enzyme'],
    #     ]
    #     filtered_combined_results = uc.execute_misc_engine(
    #         input_file = combined_results,
    #         engine='filter_csv',
    #     )
    #     all_result_files.append(filtered_combined_results)
    
    # results_all_files = uc.execute_misc_engine(
    #     input_file = all_result_files,
    #     engine='merge_csvs',
    #     merge_duplicates=True,
    # )

    # uc.params.update({
    #     'validation_score_field': 'combined PEP',
    #     'bigger_scores_better': False,
    #     'num_compared_psms': 10,
    #     'accept_conflicting_psms': False,
    #     'threshold_is_log10': True,
    #     'score_diff_threshold': 1,
    #     'psm_defining_colnames': [
    #         'Spectrum Title',
    #         'Sequence',
    #     ],
    # })

    # sanitized_combined_results = uc.execute_misc_engine(
    #     input_file = results_all_files,
    #     engine='sanitize_csv',
    # )

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(main.__doc__)
        exit()
    main(
        folder                = sys.argv[1],
        enzyme                = sys.argv[2],
        target_decoy_database = sys.argv[3],
    )
