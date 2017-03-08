# -*- coding: utf-8 -*-
"""
Created on Mon May 02 10:31:09 2016

@author: Claire
"""
from lxml import etree
import pandas as pd
import argparse

## This searches all of a downloaded HMDB xml file for each of your input masses

## things this code needs to do:
##1. Read the entire HMDB into a dictionary
#1.1 Fields of interest: HMDB ID, mass, name, synonyms, biofluid, biofunction
#2. Convert mz's to neutral masses, using golden dictionary
#2.1   probably should do both -H nd +Cl
#3. Search all masses against the masses in HMDB
#4. Write results to a file

def parse_HMDB(xml_file):
    """ Parses a cleaned up HMDB xml file into a dictionary.
    
    Dictionary contains: {HMDB_ID: {'neutral_mass': m, 'name': name, 
    'synonyms': [list of names], 'biofluid': [list of fluids], 
    'biofunction': [list of functions]}, ...}
    
    Input: xml_file: this file should be cleaned up using 
    remove_excess_xml_declarations.py
        Why? Because downloading all HMDB metabolites actually downloads
        45,000 separate xml files. Will concatenate into one file, but need to 
        remove new XML declarations.
    
    Output: dictionary
    """
    outdict = {}

    tree = etree.iterparse(xml_file, tag='metabolite')
    counter = 0
    
    for event, elem in tree:
        counter += 1
        if counter % 2000 == 0:
            print('Reading {}th metabolite from HMDB'.format(counter))
        tmpdict = {}
        tmpdict['name'] = elem.findtext('name')
        if elem.findtext('monisotopic_moleculate_weight'):
            tmpdict['neutral_mass'] = elem.findtext('monisotopic_moleculate_weight')
        else:
            tmpdict['neutral_mass'] = '0.0'
        tmpdict['synonyms'] = [i.text for i in elem.find('synonyms').findall('synonym')]
        tmpdict['biofluids'] = [i.text for i in elem.find('biofluid_locations').findall('biofluid')]
        tmpdict['biofunctions'] = [i.text for i in elem.find('ontology').find('biofunctions').findall('biofunction')]
        tmpdict['description'] = elem.findtext('description')
        # Get taxonomy info, substituents, and alternative_parents. Note - not all HMDB metabolites have this info
        try:
            tmpdict['kingdom'] = elem.find('taxonomy').findtext('kingdom')
        except:
            tmpdict['kingdom'] = ''
        try: 
            tmpdict['super_class'] = elem.find('taxonomy').findtext('super_class')
        except:
            tmpdict['super_class'] = ''
        try:
            tmpdict['class'] = elem.find('taxonomy').findtext('class')
        except:
            tmpdict['class'] = ''
        try:
            tmpdict['sub_class'] = elem.find('taxonomy').findtext('sub_class')
        except:
            tmpdict['sub_class'] = ''
        try:
            tmpdict['molecular_framework'] = elem.find('taxonomy').findtext('molecular_framework')
        except:
            tmpdict['molecular_framework'] = ''
        try:
            tmpdict['alternative_parents'] = [i.text for i in elem.find('taxonomy').find('alternative_parents').findall('alternative_parent')]
        except:
            tmpdict['alternative_parents'] = ['']
        try:
            tmpdict['substituents'] = [i.text for i in elem.find('taxonomy').find('substituents').findall('substituent')]        
        except:
            tmpdict['substituents'] = ['']       
        
        # Find diseases and pathways. Note that pathways also can have different IDs like KEGG
        try:
            tmpdict['diseases'] = [i.findtext('name') for i in elem.find('diseases').findall('disease')]
        except:
            tmpdict['diseases'] = ['']
        try:
            tmpdict['pathways'] = [i.findtext('name') for i in elem.find('pathways').findall('pathway')]
        except:
            tmpdict['pathways'] = ['']
            
        outdict[elem.findtext('accession')] = tmpdict
        elem.clear()
        
    return outdict

def extract_mzs(fname, colname=None, indexname=None, sep='\t', header=0):
    """ Extracts the mz's from a file containing mz's in a column
    named colname. Each mz is labeled with values in column labeled indexname.
    
    sep indicates the separator for columns in file fname. Default is tab-delimited.
    
    If header=None, file assumes no header row. Else it reads in the first row as the header.    
    
    If file only contains neutral masses (i.e. it's not a dataframe) and 
    colname = None, just reads in the file with each mz on a newline, and labels
    each mz with a counter.
    """
    ## Read in mz's from file
    if not colname:
        with open(fname, 'r') as f:
            mzs = f.readlines()
        mzs = [m.strip() for m in mzs]
    else:
        feattable = pd.read_csv(fname, sep=sep)
        mzs = list(feattable[colname])
    
    # Read in mz labels from file
    if not colname and not indexname:
        mznames = [str(i) for i in range(0, len(mzs))]
    else:
        feattable = pd.read_csv(fname, sep=sep, index_col=0)
        if not indexname:
            mznames = list(feattable.index)
        else:
            mznames = list(feattable(indexname))
    
    return mzs, mznames    
    
def calculate_neutral_masses(mzs):
    """ Calculates the neutral masses from a list of mz's provided.
    
    Currently returns two lists of neutral masses: one corresponding to the
    [M-H]- adduct and one to the [M+Cl]- adduct
    """
    
    ## Calculate neutral masses for both types of common neg ions.
    # Numbers include electron masses    
    # [M-H]- adducts
    minush = [float(i) + 1.007276 for i in mzs]
    # [M+Cl]- adducts    
    pluscl = [float(i) - 34.969402 for i in mzs]
    
    return minush, pluscl

def get_hmdb_hits(masses, hmdb_dict, ppm_tol):
    """ For each mass in masses, scan the given hmdb dictionary for any hits within
    ppm_tol tolerance.
    
    INPUTS
    masses = list of neutral masses
    
    hmdb_dict = dictionary created by running parse_HMDB() on an HMDB xml file. 
    
    ppm_tol = maximum error tolerated to consider mass a "hit", in ppm's
    
    OUTPUT
    allhits = dictionary containing all the hits for each mass
              dictionary contains: {mass: {hmdb_id: hmdb_dict}}
                         each hmdb_dict is as returned by parse_HMDB, and has as top key
                         the HMDB_ID. It also has one additional subkey, which is ppm.
    """
    allhits = {}
    count = 0
    print('Getting hits for {} unique neutral masses'.format(len(set(masses))))
    for m in set(masses):
        count += 1
        if count % 100 == 0:
            print('Getting hits for the {}th neutral mass'.format(count))
        mhits = [i for i in hmdb_dict if abs(float(hmdb_dict[i]['neutral_mass']) - float(m)) <= ppm_tol*float(m)/1e6]
        allhits[m] = {i: hmdb_dict[i] for i in mhits}
        # Add in ppm for each hit
        for i in allhits[m]:
            allhits[m][i]['ppm'] = str(abs(float(hmdb_dict[i]['neutral_mass']) - float(m))/float(m) * 1e6)
    
    return allhits

def write_hits_to_file(adduct_type, neutral_masses, mzs, mznames, allhits, fout, overwrite=False, describe=False):
    """ Write all of the hits for each mass for each adduct to fout.
    
    INPUTS
    adduct_type = list of adduct types of len(neutral_masses)
    neutral_masses = list of neutral masses
    mzs = list of mz's corresponding of len(neutral_masses). mzs[i] = neutral_masses[i] with adduct[i]
    allhits = dictionary of HMDB hits. Keys are all of the unique neutral masses
            allhits[neutral_masses[i]] = {hmdb_id: hmdb_dict[hmdb_id]}
    describe = whether to include the HMDB description in the output file. This feature is buggy
               and has weird formatting at times. Default is False
    overwrite = whether to overwrite the fout or append to existing file fout          
    
    OUTPUTS
    None.
    Writes tab-delimited file to fout with each hit for each mz on its own line.
    """
    
    if overwrite:
        readtype = 'w'
    else:
        readtype = 'a'
    
    with open(fout, readtype) as f:
        if overwrite:
            if describe:
                f.write('\t'.join(['featname', 'mz', 'adduct', 'neutral_mass', 
                                   'hmdb_id', 'hmdb_name', 'monoisotopic_mass', 
                                   'ppm', 'synonyms', 'biofluids', 'biofunctions', 
                                   'description', 'kingdom', 'super_class',
                                   'class', 'sub_class', 'molecular_framework', 
                                   'alternative_parents', 'substituents',
                                   'pathways', 'diseases']) + '\n')
            else:
                f.write('\t'.join(['featname', 'mz', 'adduct', 'neutral_mass', 
                                   'hmdb_id', 'hmdb_name', 'monoisotopic_mass', 
                                   'ppm', 'synonyms', 'biofluids', 'biofunctions',
                                   'kingdom', 'super_class',
                                   'class', 'sub_class', 'molecular_framework', 
                                   'alternative_parents', 'substituents',
                                   'pahways', 'diseases']) + '\n')
                                   
        for a, nm, mz, mzname in zip(adduct_type, neutral_masses, mzs, mznames):
            # If the neutral mass has an HMDB hit
            if allhits[nm]:
                for h_id in allhits[nm]:
                    # Assemble data from the dictionary
                    h_name = allhits[nm][h_id]['name'].encode('utf-8')
                    ppm = allhits[nm][h_id]['ppm']
                    mono_mass = allhits[nm][h_id]['neutral_mass']
                    syns = '; '.join(allhits[nm][h_id]['synonyms']).encode('utf-8')
                    fluids = '; '.join(allhits[nm][h_id]['biofluids']).encode('utf-8')
                    functions = '; '.join(allhits[nm][h_id]['biofunctions']).encode('utf-8')
                    kingdom = allhits[nm][h_id]['kingdom']
                    super_class = allhits[nm][h_id]['super_class']
                    cls = allhits[nm][h_id]['class']
                    sub_class = allhits[nm][h_id]['sub_class']
                    molec_framework = allhits[nm][h_id]['molecular_framework']
                    alt_parents = '; '.join(allhits[nm][h_id]['alternative_parents']).encode('utf-8')
                    substituents = '; '.join(allhits[nm][h_id]['substituents']).encode('utf-8')
                    pathways = '; '.join(allhits[nm][h_id]['pathways']).encode('utf-8')
                    diseases = '; '.join(allhits[nm][h_id]['diseases']).encode('utf-8')
                    
                    if describe:
                        description = allhits[nm][h_id]['description'].encode('utf-8')
                        f.write('\t'.join([str(i) for i in [mzname, mz, a, nm, 
                                                           h_id, h_name, mono_mass, 
                                                           ppm, syns, fluids, functions, 
                                                           description, kingdom,
                                                           super_class, cls,
                                                           sub_class, molec_framework,
                                                           alt_parents, 
                                                           substituents,
                                                           pathways, diseases]]) + '\n')
                    else:
                        f.write('\t'.join([str(i) for i in [mzname, mz, a, nm, 
                                                           h_id, h_name, mono_mass, 
                                                           ppm, syns, fluids, functions, 
                                                           kingdom,
                                                           super_class, cls,
                                                           sub_class, molec_framework,
                                                           alt_parents, 
                                                           substituents,
                                                           pathways, diseases]]) + '\n')
            else:
                f.write('\t'.join([str(i) for i in [mzname, mz, a, nm]]) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', help='clean hmdb xml file [required, default = data/hmdb_metabolites_clean.xml]', default='data/hmdb_metabolites_clean.xml', required=True)
    parser.add_argument('-f', help='input feature table. has column labeled mz and feature names in the first column [required, no default]', required=True)
    parser.add_argument('-s', help='feature table separator. default is tab-delimited', default='\t')
    parser.add_argument('-p', help='ppm tolerance. default is 5', default=5, type=int)
    parser.add_argument('-o', help='output file to write results to [required, no default]', required=True)
    parser.add_argument('-d', help='whether to include description in output (may mess up formatting)', default=False)
    args = parser.parse_args()
    
    hmdb_xml = args.x
    feattable = args.f
    sep = args.s
    ppm_tolerance = args.p
    outfile = args.o
    describe = args.d

    ## 1. Parse HMDB xml
    hmdb_dict = parse_HMDB(hmdb_xml)
#    try:
#        hmdb_dict = parse_HMDB(hmdb_xml)
#    except:
#        print('failed to parse hmdb')
        
    ## 2. Extract mz's from feature table file
    try:
        mzs, mznames = extract_mzs(feattable, colname='mz', indexname=None, sep=sep, header=0)
    except:
        print('failed to extract mzs from feature table')
        
    ## 3. Calculate neutral masses for each mz
    print('Calculating neutral masses')
    minush_nm, pluscl_nm = calculate_neutral_masses(mzs)
    
    ## 4. Get hits to HMDB, within a tolerance
    print('Getting hits to HMDB metabolites, [M-H]-')
    hits_h = get_hmdb_hits(minush_nm, hmdb_dict, ppm_tolerance)
    print('Getting hits to HMDB metabolites, [M+Cl]-')
    hits_cl = get_hmdb_hits(pluscl_nm, hmdb_dict, ppm_tolerance)
    
    ## 5. Write to file
    print('Writing file ' + outfile)
    # 5.1 [M-H]- adducts first
    adducts = ['[M-H]-']*len(minush_nm)
    newfile = True
    write_hits_to_file(adducts, minush_nm, mzs, mznames, hits_h, outfile, overwrite=newfile, describe=describe)
    
    # [M+Cl]- next
    adducts = ['[M+Cl]-']*len(pluscl_nm)
    newfile = False
    write_hits_to_file(adducts, pluscl_nm, mzs, mznames, hits_cl, outfile, overwrite=newfile, describe=describe)
