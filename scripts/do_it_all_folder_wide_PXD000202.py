#!/usr/bin/env python3
# encoding: utf-8
import ursgal
import sys
import glob
import os

offsets = {
    'SACO1':{
        '110224_OV3_P3_BH_HA_SACO_1_10.mzML':0 ,
        '110224_OV3_P3_BH_HA_SACO_1_1.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_1_2.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_1_3.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_1_4.mzML': 2,
        '110224_OV3_P3_BH_HA_SACO_1_5.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_1_6.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_1_7.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_1_8.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_1_9.mzML': 0,
    },
    'SACO2':{
        '110224_OV3_P3_BH_HA_SACO_2_10.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_2_1.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_2_2.mzML': -2,
        '110224_OV3_P3_BH_HA_SACO_2_3.mzML': 2,
        '110224_OV3_P3_BH_HA_SACO_2_4.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_2_5.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_2_6.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_2_7.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_2_8.mzML': 0,
        '110224_OV3_P3_BH_HA_SACO_2_9.mzML': 0,
    },
    'SA1':{
        '110317_OV3_P3_BH_HA_SA1_01.mzML': -2,
        '110317_OV3_P3_BH_HA_SA1_02.mzML': -2,
        '110317_OV3_P3_BH_HA_SA1_03.mzML': -2,
        '110317_OV3_P3_BH_HA_SA1_04.mzML': 0,
        '110317_OV3_P3_BH_HA_SA1_05_110318144238.mzML': 0,
        '110317_OV3_P3_BH_HA_SA1_06.mzML': 0,
        '110317_OV3_P3_BH_HA_SA1_07.mzML': 0,
        '110317_OV3_P3_BH_HA_SA1_08.mzML': 0,
        '110317_OV3_P3_BH_HA_SA1_09.mzML': 0,
        '110317_OV3_P3_BH_HA_SA1_10.mzML': 0,
    },
    'SA2':{
        '110318_OV3_P3_BH_HA_SA2_01.mzML': 0,
        '110318_OV3_P3_BH_HA_SA2_02.mzML': -2,
        '110318_OV3_P3_BH_HA_SA2_03.mzML': -2,
        '110318_OV3_P3_BH_HA_SA2_04.mzML': 0,
        '110318_OV3_P3_BH_HA_SA2_05.mzML': 0,
        '110318_OV3_P3_BH_HA_SA2_06.mzML': 0,
        '110318_OV3_P3_BH_HA_SA2_07.mzML': 0,
        '110318_OV3_P3_BH_HA_SA2_08.mzML': 0,
        '110318_OV3_P3_BH_HA_SA2_09_110320172044.mzML': 0,
        '110318_OV3_P3_BH_HA_SA2_10.mzML': 0,
    },
    'SA3':{
        '110318_OV3_P3_BH_HA_SA3_01.mzML': 0,
        '110318_OV3_P3_BH_HA_SA3_02.mzML': -2,
        '110318_OV3_P3_BH_HA_SA3_03.mzML': 0,
        '110318_OV3_P3_BH_HA_SA3_04.mzML': 0,
        '110318_OV3_P3_BH_HA_SA3_05.mzML': 0,
        '110318_OV3_P3_BH_HA_SA3_06.mzML': 0,
        '110318_OV3_P3_BH_HA_SA3_07.mzML': 0,
        '110318_OV3_P3_BH_HA_SA3_08.mzML': 0,
        '110318_OV3_P3_BH_HA_SA3_09.mzML': 0,
        '110318_OV3_P3_BH_HA_SA3_10.mzML': 0,
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
        'precursor_mass_tolerance_minus' : 5,
        'precursor_mass_tolerance_plus': 5,
        'frag_mass_tolerance' : 0.2,
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
