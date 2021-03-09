#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os

offsets = {
    'Lana-SILAC-1raw':{
        'Lana-SILAC-10.mzML': -2,
        'Lana-SILAC-11.mzML': -2,
        'Lana-SILAC-12.mzML': -2,
        'Lana-SILAC-13.mzML': 0,
        'Lana-SILAC-1.mzML': 0,
        'Lana-SILAC-2.mzML': -6,
        'Lana-SILAC-2-1_1.mzML': 0,
        'Lana-SILAC-3.mzML': -8,
        'Lana-SILAC-4.mzML': -6,
        'Lana-SILAC-5.mzML': -4,
        'Lana-SILAC-6.mzML': -6,
        'Lana-SILAC-7.mzML': -6,
        'Lana-SILAC-8.mzML': -2,
        'Lana-SILAC-9.mzML': -6,
    },
    'Lana-SILAC-2raw':{
        'Lana-SILAC-2-10.mzML': 0,
        'Lana-SILAC-2-11.mzML': 0,
        'Lana-SILAC-2-12.mzML': 0,
        'Lana-SILAC-2-13.mzML': 0,
        'Lana-SILAC-2-1_2.mzML': 0,
        'Lana-SILAC-2-2-1.mzML': 0,
        'Lana-SILAC-2-2.mzML': -4,
        'Lana-SILAC-2-3.mzML': 2,
        'Lana-SILAC-2-4.mzML': 0,
        'Lana-SILAC-2-5.mzML': 2,
        'Lana-SILAC-2-6.mzML': 2,
        'Lana-SILAC-2-7.mzML': 4,
        'Lana-SILAC-2-8.mzML': 4,
        'Lana-SILAC-2-9.mzML': 4,
    },
    'Lana-SILAC-3raw':{
        'Lana-SILAC-3-10.mzML': 0,
        'Lana-SILAC-3-11.mzML': 2,
        'Lana-SILAC-3-12.mzML': 0,
        'Lana-SILAC-3-13.mzML': 2,
        'Lana-SILAC-3-1.mzML': 2,
        'Lana-SILAC-3-2-1.mzML': 6,
        'Lana-SILAC-3-2.mzML': 2,
        'Lana-SILAC-3-3.mzML': 6,
        'Lana-SILAC-3-4.mzML': 4,
        'Lana-SILAC-3-5.mzML': 4,
        'Lana-SILAC-3-6.mzML': 4,
        'Lana-SILAC-3-7.mzML': 4,
        'Lana-SILAC-3-8.mzML': 0,
        'Lana-SILAC-3-9.mzML': 0,
    },
    'Lana-SILAC-4raw':{
        'Lana-SILAC-4-10.mzML': 0,
        'Lana-SILAC-4-11.mzML': 0,
        'Lana-SILAC-4-12.mzML': 0,
        'Lana-SILAC-4-13.mzML': 0,
        'Lana-SILAC-4-1.mzML': 0,
        'Lana-SILAC-4-2-1.mzML': -2,
        'Lana-SILAC-4-2.mzML': -2,
        'Lana-SILAC-4-3.mzML': -2,
        'Lana-SILAC-4-4.mzML': -2,
        'Lana-SILAC-4-5.mzML': 0,
        'Lana-SILAC-4-6.mzML': 0,
        'Lana-SILAC-4-7.mzML': 0,
        'Lana-SILAC-4-8.mzML': 2,
        'Lana-SILAC-4-9.mzML': -2,
    },
}

def main(folder=None, enzyme=None, target_decoy_database=None):
    '''
    Workflow for the analysis a dataset with one run per sample.
    Usage:
        python <script_name.py> <folder_with_mzML> <enzyme> <path_to_database>
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
            exit()
    for sample in offset_files:
        if sample not in mzML_files:
            print('Sample in offset dict but mzML file NOT in folder: {}'.format(sample))
            exit()

    mass_spectrometer = 'QExactive+'
    search_engines  = [
       'xtandem_vengeance',
       'msfragger_2_3',
       'msgfplus_v2019_07_03',
    ]

    validation_engine = 'percolator_3_4_0'

    params = {
        'database'      : target_decoy_database,
        'enzyme'        : enzyme,
        'precursor_mass_tolerance_minus' : 10,
        'precursor_mass_tolerance_plus': 10,
        'frag_mass_tolerance' : 20,
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
    
    combined_pep_result_files = []
    for n, sample in enumerate(sorted(offsets.keys(), reverse=True)):
        validated_result_files = []
        for search_engine in search_engines:
            engine_results_validated = []
            for n, mod in enumerate(Hvo_Glyco):
                results_one_mod = []
                for spec_file in sorted(offsets[sample].keys()):
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
                            'C,fix,any,Methylthio',
                            'M,opt,any,Oxidation',
                            'P,opt,any,Label:13C(5)',
                            'K,opt,any,Label:13C(6)15N(2)',
                            'R,opt,any,Label:13C(6)',
                            '*,opt,Prot-N-term,Acetyl',
                        ]
                    else:
                        uc.params['modifications'] = [
                            'C,fix,any,Methylthio',
                            'M,opt,any,Oxidation',
                            'P,opt,any,Label:13C(5)',
                            'K,opt,any,Label:13C(6)15N(2)',
                            'R,opt,any,Label:13C(6)',
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

                    results_one_mod.append(unified_search_results)

                uc.params['prefix'] = sample
                merged_1engine_1mod_1sample = uc.execute_misc_engine(
                    input_file = results_one_mod,
                    engine='merge_csvs',
                    # merge_duplicates=True,
                )
                uc.params['prefix'] = ''
                # engine_results_unvalidated.append(merged_1engine_1mod_1sample)

                validated_csv = uc.validate(
                    input_file  = merged_1engine_1mod_1sample,
                    engine      = validation_engine,
                )
                engine_results_validated.append(validated_csv)

            merged_1engine_all_mods_validated = uc.execute_misc_engine(
                input_file = engine_results_validated,
                engine='merge_csvs',
                merge_duplicates=False,
            )
            validated_result_files.append(merged_1engine_all_mods_validated)

        uc.params['prefix'] = sample
        combined_pep_validated = uc.combine_search_results(
            input_files=validated_result_files,
            engine='combine_pep_1_0_0',
        )
        uc.params['prefix'] = ''
        uc.params['csv_filter_rules'] = [
            ['Is decoy', 'equals', 'false'],
            ['combined PEP','lte', 0.01],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ]
        filtered_validated_results = uc.execute_misc_engine(
            input_file = combined_pep_validated,
            engine='filter_csv',
        )
        combined_pep_result_files.append(filtered_validated_results)

        # uc.params['peptide_forest_initial_engine'] = 'msfragger_2_3'
        # uc.params['peptide_forest_file_params'] = {}
        # uc.params['prefix'] = 'peptide_forest_' + sample
        # peptide_forest_validated =  uc.validate(
        #     input_file=unvalidated_result_files,
        #     engine='peptide_forest',
        # )
        # uc.params['csv_filter_rules'] = [
        #     ['Is decoy', 'equals', 'false'],
        #     ['q-value_RF-reg','lte', 0.01],
        #     ['Conflicting uparam', 'contains_not', 'enzyme'],
        # ]
        # filtered_peptide_forest = uc.execute_misc_engine(
        #     input_file = peptide_forest_validated,
        #     engine='filter_csv',
        # )
        # peptide_forest_result_files.append(filtered_peptide_forest)

    uc.params['prefix'] = ''
    results_all_combined_pep = uc.execute_misc_engine(
        input_file = combined_pep_result_files,
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
        input_file = results_all_combined_pep,
        engine='sanitize_csv',
    )

    uc.params['prefix'] = 'Glyco_everywhere'
    uc.params['csv_filter_rules'] = [
        ['Modifications', 'contains', 'Hex'],
        # ['Sequence','contains_glycosite', 'N[^P][ST]']
    ]
    Glyco_filtered = uc.execute_misc_engine(
        input_file = sanitized_combined_results,
        engine='filter_csv',
    )

    uc.params['prefix'] = 'Glyco_glycosite'
    uc.params['csv_filter_rules'] = [
        ['Modifications', 'contains', 'Hex'],
        ['Sequence','contains_glycosite', 'N[^P][ST]']
    ]
    Glyco_filtered = uc.execute_misc_engine(
        input_file = sanitized_combined_results,
        engine='filter_csv',
    )
    uc.params['prefix'] = ''

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(main.__doc__)
        exit()
    main(
        folder                = sys.argv[1],
        enzyme                = sys.argv[2],
        target_decoy_database = sys.argv[3],
    )
