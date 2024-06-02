#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@Date : 2024/05/30
@Author : Sarah Eger

Full article on detecting sample swaps:
https://gatk.broadinstitute.org/hc/en-us/articles/360041696232-Detecting-sample-swaps-with-Picard-tools

`CrosscheckFingerprints` (picard): 
https://gatk.broadinstitute.org/hc/en-us/articles/21905112242075-CrosscheckFingerprints-Picard

STEP 2: `CrosscheckFingerprints` across multiome GEX (still need to add ATAC)

"""

import os
import glob
import numpy as np
import pandas as pd
from itertools import combinations

def picard_CrosscheckFingerprints(prefix, indir, outdir):
    # get fingerprints
    fingerprints = sorted(glob.glob(indir+'/*.fingerprint.vcf'))
    
    # make a file list
    output = os.path.join(outdir, prefix)
    df = pd.DataFrame(fingerprints, columns=['vcf'])
    df.to_csv(output+'.vcf_list', index=False, header=False, sep='\t')

    # run CrosscheckFingerprints
    cmd = ' '.join([picard,
                    'CrosscheckFingerprints',
                    '--NUM_THREADS %s' % n_threads,
                    '--CALCULATE_TUMOR_AWARE_RESULTS false',
                    '--LOD_THRESHOLD -5',
                    '--CROSSCHECK_BY SAMPLE',
                    '--HAPLOTYPE_MAP %s' % hapmap,
                    '--INPUT %s.vcf_list' % output,
                    '--MATRIX_OUTPUT %s.crosscheck_mtx' % output,
                    '--OUTPUT %s.crosscheck_metrics' % output])
    print(cmd)
    os.system(cmd)
    return output+'.crosscheck_metrics'

if __name__ ==  '__main__':
    # number of threads
    n_threads = 10

    PREFIX = 'COL_GEX_FTX'

    gatk = '/home/eger/software/gatk-4.5.0.0/gatk'
    picard = 'java -jar /home/eger/software/picard.jar'

    REF_DIR1 = '/home/eger/references/10xgenomics/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta'
    REF_DIR2 = '/home/eger/references/GATK'

    DIR1 = '/home/eger/projects/multiome/v1/1_cellranger'
    DIR2 = '/scratch/eger/projects/multiome/v1/4_Fingerprinting'
    SUBDIR2a = DIR2+'/reheadered_bams'
    SUBDIR2b = DIR2+'/ExtractFingerprint_vcfs'
    SUBDIR2c = DIR2+'/CrosscheckFingerprints_metrics'

    if not os.path.exists(SUBDIR2c):
        os.makedirs(SUBDIR2c)
    
    # hap map: from GATK bundle, but header replaced with cellranger bam header
    hapmap = os.path.join(REF_DIR2, 'Homo_sapiens_assembly38.haplotype_database.GEX.txt')

    # run cross check
    fname1 = picard_CrosscheckFingerprints(PREFIX, SUBDIR2b, SUBDIR2c)
    df1 = pd.read_csv(fname1, sep='\t', comment='#')
    
    # get the individual ID
    df1['LEFT_INDIV'] = df1['LEFT_SAMPLE'].apply(lambda x: x.split('_')[1])
    df1['RIGHT_INDIV'] = df1['RIGHT_SAMPLE'].apply(lambda x: x.split('_')[1])

    # remove the duplicate rows & the ones that compare the fingerprint to itself
    df1['sample_pair'] = df1.apply(lambda x: '|'.join(sorted([x['LEFT_SAMPLE'], x['RIGHT_SAMPLE']])), axis=1)
    df1 = df1.drop_duplicates(subset=['sample_pair'], keep='first')
    df1 = df1.drop('sample_pair', axis=1)
    df1 = df1[df1['RESULT'] != 'EXPECTED_MATCH']
    df1 = df1.sort_values(['LEFT_SAMPLE', 'RIGHT_SAMPLE'])

    # create a 'FINAL_RESULT' column based on individual and sample
    df1['FINAL_RESULT'] = np.nan
    df1.loc[(df1['LEFT_INDIV'] == df1['RIGHT_INDIV']) & (df1['RESULT'] == 'UNEXPECTED_MATCH'), 'FINAL_RESULT'] = 'BATCH_MATCH'
    df1.loc[(df1['LEFT_INDIV'] != df1['RIGHT_INDIV']) & (df1['RESULT'] == 'UNEXPECTED_MATCH'), 'FINAL_RESULT'] = 'UNEXPECTED_MATCH'
    df1.loc[(df1['LEFT_INDIV'] != df1['RIGHT_INDIV']) & (df1['RESULT'] == 'EXPECTED_MISMATCH'), 'FINAL_RESULT'] = 'EXPECTED_MISMATCH'
    df1.loc[(df1['LEFT_INDIV'] == df1['RIGHT_INDIV']) & (df1['RESULT'] == 'EXPECTED_MISMATCH'), 'FINAL_RESULT'] = 'UNEXPECTED_MISMATCH'
    print(df1['FINAL_RESULT'].value_counts())

    df1 = df1[['LEFT_INDIV', 'LEFT_SAMPLE', 'RIGHT_INDIV', 'RIGHT_SAMPLE', 
            'FINAL_RESULT', 'LOD_SCORE']]

    # find the samples that match across batches
    df2 = df1[df1['FINAL_RESULT'] == 'BATCH_MATCH'].sort_values('LOD_SCORE', ascending=False)

    # find the problems
    df3 = df1[df1['FINAL_RESULT'].str.contains('UNEXPECTED')].sort_values('LOD_SCORE', ascending=False)

    # save results
    batch_match = os.path.join(SUBDIR2c, PREFIX+'_batch_match.txt')
    mismatch = os.path.join(SUBDIR2c, PREFIX+'_unexpected.txt')
    df2.to_csv(batch_match, sep='\t', index=False)
    df3.to_csv(mismatch, sep='\t', index=False)