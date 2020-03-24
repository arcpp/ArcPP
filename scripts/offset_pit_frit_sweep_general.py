#!/usr/bin/env python3
# encoding: utf-8

import ursgal
import glob
import csv
import os
from collections import defaultdict as ddict
import sys
import re


def main(folder=None, enzyme=None, target_decoy_database=None, analyze_only=False):
    ''' 
    '''
    mass_spectrometer = 'QExactive+'
    search_engine = 'xtandem_vengeance',
    validation_engine = 'percolator_2_08'
    all_mods = [
        # 'C,fix,any,Carbamidomethyl',
        'M,opt,any,Oxidation',
        '*,opt,Prot-N-term,Acetyl',
    ]

    # Initializing the Ursgal UController class with
    # our specified modifications and mass spectrometer
    params = {
        'database'      : target_decoy_database,
        'modifications' : all_mods,
        'enzyme'        : enzyme,
        'csv_filter_rules' : [
            ['Is decoy', 'equals', 'false'],
            ['PEP', 'lte', 0.01],
            ['Conflicting uparam', 'contains_not', 'enzyme'],
        ],
        'precursor_mass_tolerance_minus' : 5,
        'precursor_mass_tolerance_plus': 5,
        'frag_mass_tolerance' : 0.4,
        'frag_mass_tolerance_unit': 'da',
        'rounded_mass_decimals': 2,
        '-xmx': '32g',
        'peptide_mapper_class_version' : 'UPeptideMapper_v4',
        'use_pyqms_for_mz_calculation' : False,
        'scan_skip_modulo_step': 2,
    }

    if analyze_only:
        mzML_basename_list = set()
        for mzML_path in glob.glob(os.path.join(folder, '*.idx.gz')):
            mzML_basename = os.path.basename(mzML_path)
            mzML_basename_list.add(mzML_basename)
    else:
        mzML_basename_list = search(
            input_folder=folder,
            params=params,
            search_engine=search_engine,
        )
    analyze(
        folder=folder,
        params=params,
        validation_engine=validation_engine,
        mzML_basename_list = mzML_basename_list,
    )

def search(input_folder=None, params=None, search_engine=None,):
    '''
    Does the parameter sweep on every nth (defined by scan_skip_modulo_step) 
    MS2 spectrum of the given file and engine.

    Sweeps over:

        * machine_offset_in_ppm from -10 to +10 ppm offset
        * precursor mass tolerance from 5 to 20 ppm
        * fragment mass tolerance from 5 to 40 ppm

    The search can be very time consuming (depending on your machine/cluster),
    therefore the analyze step can be performed separately by calling analyze()
    instead of search() when one has already performed the searches and wants
    to analyze the results.
    '''

    precursor_ion_tolerance_list = [5, 7.5, 10, 20]
    frag_ion_tolerance_list = [5, 7.5, 10, 20, 40]

    R = ursgal.UController(
        params=params
    )

    mzML_basename_list = []
    for mzML_path in glob.glob(os.path.join(input_folder, '*.idx.gz')):
        mzML_basename = os.path.basename(mzML_path)
        mzML_basename_list.append(mzML_basename)
        for ppm_offset in range(-10, 12, 2):
            R.params['machine_offset_in_ppm'] = ppm_offset
            R.params['prefix'] = 'ppm_offset_{0}'.format(
                int(ppm_offset)
            )
            mgf_file = R.convert(
                input_file=mzML_path,
                engine='mzml2mgf_2_0_0'
            )
            for precursor_ion_tolerane in precursor_ion_tolerance_list:
                for frag_ion_tolerance in frag_ion_tolerance_list:

                    new_prefix = '_pit_{0}_fit_{1}'.format(
                        precursor_ion_tolerane,
                        frag_ion_tolerance
                    )
                    R.params[
                        'precursor_mass_tolerance_minus'] = precursor_ion_tolerane
                    R.params[
                        'precursor_mass_tolerance_plus'] = precursor_ion_tolerane
                    R.params['frag_mass_tolerance'] = frag_ion_tolerance
                    R.params['prefix'] = new_prefix

                    search_result = R.search_mgf(
                        input_file = mgf_file,
                        engine = search_engine,                        )

                    converted_result = R.convert(
                        input_file=search_result,
                        guess_engine = True,
                    )
                    
                    try:
                        mapped_results = R.execute_misc_engine(
                            input_file=converted_result,
                            engine='upeptide_mapper',
                        )
                    except:
                        continue

                    unified_search_results = R.execute_misc_engine(
                        input_file = mapped_results,
                        engine='unify_csv'
                    )

    return mzML_basename_list


def analyze(folder=None, validation_engine=None, params=None, mzML_basename_list=None):
    '''
    Parses the result files form search and write a result .csv
    containing the number of identifications for each file and setting
    '''

    R = ursgal.UController(
        params=params
    )
    csv_collector = {}
    sample_offset_combos = []

    all_tested_offsets = [str(n) for n in range(-10, 12, 2)]

    for mzML_file in mzML_basename_list:
        for theo_offset in all_tested_offsets:
            sample_offset_combos.append((mzML_file, theo_offset))

    for csv_path in glob.glob(os.path.join('{0}'.format(folder), '*', '*_unified.csv')):
        print('found unified csv file:', csv_path)
        dirname = os.path.dirname(csv_path)
        splitted_basename = os.path.basename(csv_path).split('_')
        offset = splitted_basename[7]
        precursor_ion_tolerance = splitted_basename[2]
        frag_ion_tolerance = splitted_basename[4]
        prefix = '_'.join(splitted_basename[:8])
        mzML_basename = '_'.join(splitted_basename[8:-4]) + '.idx.gz'
        print('corresponding to mzML file:', mzML_basename)
        if mzML_basename not in mzML_basename_list:
            print(mzML_basename, 'not found in mzML_basename_list')

        R.params['machine_offset_in_ppm'] = offset
        R.params['precursor_mass_tolerance_minus'] = precursor_ion_tolerance
        R.params['precursor_mass_tolerance_plus'] = precursor_ion_tolerance
        R.params['frag_mass_tolerance'] = frag_ion_tolerance
        R.params['prefix'] = prefix

        validated_path = csv_path.replace(
            '_unified.csv',
            '_{0}_validated.csv'.format(validation_engine)
        )
        if os.path.exists(validated_path):
            csv_path = validated_path
            print('found validation file:', validated_path)
        else:
            try:
                csv_path = R.validate(
                    input_file=csv_path,
                    engine=validation_engine
                )
                print('validation successful:', csv_path)
            except:
                print('----- validation failed --------')
                continue

        pit_fit = (precursor_ion_tolerance, frag_ion_tolerance)

        if pit_fit not in csv_collector.keys():
            csv_collector[pit_fit] = ddict(set)

        csv_key = (sample, offset)

        print('Reading file: {0}'.format(csv_path))
        for line_dict in csv.DictReader(open(csv_path, 'r')):
            if line_dict['Is decoy'] == 'true':
                continue
            if float(line_dict['PEP']) <= 0.01:
                csv_collector[pit_fit][csv_key].add(
                    '{0}{1}'.format(
                        line_dict['Sequence'],
                        line_dict['Modifications']
                    )
                )
                csv_collector[pit_fit]['set_of_all'].add(
                    '{0}{1}'.format(
                        line_dict['Sequence'],
                        line_dict['Modifications']
                    )
                )

    fieldnames = [
        'Sample',
        'tested_ppm_offset',
        'peptide_count'
    ]

    for pit_fit in csv_collector.keys():
        outfile_name = 'data_ppm_sweep_precursor_mass_tolerance_{0}_fragment_mass_tolerance_{1}.csv'.format(*pit_fit)
        out_file_path = os.path.join(folder, outfile_name)
        with open(out_file_path, 'w') as io:
            csv_writer = csv.DictWriter(io, fieldnames)
            csv_writer.writeheader()

            # write missing values
            for sample_offset in sample_offset_combos:
                sample, ppm_offset = sample_offset
                if sample_offset not in csv_collector[pit_fit].keys():
                    dict_2_write = {
                        'Sample': sample,
                        'tested_ppm_offset': ppm_offset,
                        'peptide_count': 0,
                    }
                    csv_writer.writerow(dict_2_write)

            for csv_key, peptide_set in csv_collector[pit_fit].items():
                if csv_key != 'set_of_all':
                    sample, ppm_offset = csv_key
                    dict_2_write = {
                        'Sample': sample,
                        'tested_ppm_offset': ppm_offset,
                        'peptide_count': len(peptide_set),
                    }
                    csv_writer.writerow(dict_2_write)
                else:
                    print('>>>>>>>>>>>>>>>>>>>>>>')
                    print('pit: {0} ; fit: {1}'.format(*pit_fit))
                    print('peptide_count:', len(peptide_set))
                    print('<<<<<<<<<<<<<<<<<<<<<<')
    return

if __name__ == '__main__':
    if sys.argv[4] == 'True':
        analyze_only = True
    else:
        analyze_only = False
    main(
        folder=sys.argv[1],
        enzyme=sys.argv[2],
        target_decoy_database=sys.argv[3],
        analyze_only=analyze_only
    )
