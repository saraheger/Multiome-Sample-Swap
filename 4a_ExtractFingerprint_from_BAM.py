#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@Date : 2024/05/30
@Author : Sarah Eger

Full article on detecting sample swaps:
https://gatk.broadinstitute.org/hc/en-us/articles/360041696232-Detecting-sample-swaps-with-Picard-tools

STEP 1: `ExtractFingerprint` (gatk) to extract fingerprint from a cellranger BAM
This makes the STEP 2 faster because the fingerprints are VCFs.

"""

import os
import glob
from multiprocessing import Pool

def reheader_addVN(bam_path, sample_alias):
    output = os.path.join(SUBDIR2a, sample_alias+'.reheader.bam')
    new_header = r'@HD\tVN:1.4\tSO:coordinate\0'

    cmd = (f'(printf "{new_header}" && tail -c +13 "{bam_path}") | '
           f'dd conv=notrunc of="{output}" bs=1 seek=8')

    return cmd

def samtools_reheader_addVN(bam_path, sample_alias):
    # create sample alias
    case = bam_path.split('/')[-3]
    lib_type = bam_path.split('/')[-1].split('_')[0].upper()
    batch = bam_path.split('/')[-4].split('Batch')[1]
    alias = 'B'+str(batch)+'_'+case+'_'+lib_type

    output = os.path.join(SUBDIR2a, sample_alias+'.reheader.bam')
    cmd = (f'samtools view -H "{bam_path}" | '
           f'sed "s/@HD\s*\(.*\)/@HD\tVN:1.4\tSO:coordinate/" | '
           f'samtools reheader - "{bam_path}" > "{output}"')
    os.system(cmd)

    # index the reheadered bam
    cmd = f'samtools index "{output}"'
    os.system(cmd)

def gatk_ExtractFingerprint(bam_path, sample_alias):
    output = os.path.join(SUBDIR2b, sample_alias+'.fingerprint.vcf')
    cmd = ' '.join([gatk,
                    'ExtractFingerprint',
                    '--SAMPLE_ALIAS %s' % sample_alias,
                    '--HAPLOTYPE_MAP %s' % hap_map,
                    '--REFERENCE_SEQUENCE %s' % ref_fasta,
                    '--INPUT %s' % bam_path,
                    '--OUTPUT %s' % output])
    return cmd

def run(cmd):
    print(cmd)
    os.system(cmd)

if __name__ ==  '__main__':
    # number of parallel tasks
    tasks = 4
    
    lib_type = 'GEX' # doesn't have to be reheadered like ATAC

    gatk = '/home/eger/software/gatk-4.5.0.0/gatk'
    picard = 'java -jar /home/eger/software/picard.jar'

    REF_DIR1 = '/home/eger/references/10xgenomics/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta'
    REF_DIR2 = '/home/eger/references/GATK'

    DIR1 = '/home/eger/projects/multiome/v1/1_cellranger'
    DIR2 = '/scratch/eger/projects/multiome/v1/4_Fingerprinting'
    SUBDIR2a = DIR2+'/reheadered_bams'
    SUBDIR2b = DIR2+'/ExtractFingerprint_vcfs'

    for dir in [DIR2, SUBDIR2a, SUBDIR2b]:
        if not os.path.exists(dir):
            os.makedirs(dir)

    # reference fasta used by cellranger
    ref_fasta = os.path.join(REF_DIR1, 'genome.fa')

    # hap map: from GATK bundle, but header replaced with cellranger bam header
    hapmap = os.path.join(REF_DIR2, 'Homo_sapiens_assembly38.haplotype_database.GEX.txt')

    # print and save start time
    start_time = os.system('date')
    print(start_time)

    # get lists of BAMs
    BAMs = sorted(glob.glob(DIR1+'/Batch*/*/outs/'+lib_type.lower()+'_possorted_bam.bam'))

    # loop through bams to create a list of commands
    cmds = []
    for bam in BAMs:
        # create sample alias
        case = bam.split('/')[-3]
        lib_type = bam.split('/')[-1].split('_')[0].upper()
        batch = bam.split('/')[-4].split('Batch')[1]
        alias = 'B'+str(batch)+'_'+case+'_'+lib_type
        
        # check that it doesn't already exist
        if not os.path.exists(SUBDIR2a+'/'+alias+'.fingerprint.vcf'):
            cmd = gatk_ExtractFingerprint(bam, alias, lib_type)
            cmds.append(cmd)

    pool = Pool(processes = tasks)
    pool.map(run, cmds)
    pool.close()
    pool.join()

    # print and save end time
    end_time = os.system('date')
    print(end_time)

    # calculate total time
    total_time = end_time - start_time
    print(total_time)
