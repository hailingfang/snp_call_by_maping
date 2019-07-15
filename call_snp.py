#! /usr/bin/env python3

import argparse
import os
import time
import shutil
import sys
import fbio.fparse

def get_args():
    args = argparse.ArgumentParser(description='This utility is used to get SNPs \
        consitance sequence.')
    args.add_argument('-ref_genome', type=str, required=True, help='the file name \
        of refrence genome. all reads or genome of target will mapped to this genome')
    args.add_argument('-target_reads_dir', type=str, default=None, help='target reads \
        file. reads file of one item needed put into a directory and put all this item file \
        into target_genome_reads directory. Utill now, only pair end reads are implementation.')
    args.add_argument('-target_assembly_dir', type=str, default=None, help='put all target \
        genome into a single directory.')
    args.add_argument('-outgroup_assembly_dir', type=str, default=None, help='genome \
        assembly of outgroup, multiple outgroup can be added, put all assembly into a \
        single directory.')
    args.add_argument('-outgroup_reads_dir', type=str, default=None, help='genome \
        reads of outgroup. multiple genome can be used. put reads file into a directory of \
        a single item and put all item into outgroup_genome_reads directory.')
    args.add_argument('-num_threads', type=int, default=1, help='number of threads this pipeline \
        used')
    args.add_argument('-snp_qual_cutoff', type=int, default=30, help='QUAL value cutoff of vcf result, default is 30.')
    args.add_argument('-snp_coverage_cutoff', type=int, default=20, help='coverage cutoff of snp site, default is 20.')
    args.add_argument('-snp_diversity_cutoff', type=float, default=0.5, help='snp diversity cutoff. reduce \
        single nucletide varition, default is 0.5.')
    args.add_argument('-clean_up', type=str, default='No', choices=['No', 'Yes'], help='clean up \
        all tmperary file when runing this pipeline. "No" was set as default.')
    args = args.parse_args()

    ref_genome, target_reads_dir, target_assembly_dir, outgroup_assembly_dir, outgroup_reads_dir, \
        num_threads, snp_qual_cutoff, snp_coverage_cutoff, snp_diversity_cutoff = \
        args.ref_genome, args.target_reads_dir, args.target_assembly_dir, \
        args.outgroup_assembly_dir, args.outgroup_reads_dir, args.num_threads, \
        args.snp_qual_cutoff, args.snp_coverage_cutoff, args.snp_diversity_cutoff

    if (not (target_reads_dir or target_assembly_dir)) or (target_reads_dir and target_assembly_dir):
        exit('you have to input one type target data at most and least.')

    return ref_genome, target_reads_dir, target_assembly_dir, outgroup_assembly_dir, outgroup_reads_dir, \
        num_threads, snp_qual_cutoff, snp_coverage_cutoff, snp_diversity_cutoff


def target_dt_info(target_reads, target_assembly):
    data_out = {}
    if target_reads:
        target_reads = os.path.abspath(target_reads)
        data_out['reads_target'] = {}
        for item in os.listdir(target_reads):
            data_out['reads_target'][item]=[]
            for ff in os.listdir(os.path.join(target_reads,item)):
                if ff.split('.')[-1] == 'fastq' or ff.split('.')[-1] == 'fq':
                    data_out['reads_target'][item].append(os.path.join(target_reads,item,ff))
    elif target_assembly:
        target_assembly = os.path.abspath(target_assembly)
        data_out['assembly_target'] = []
        for item in os.listdir(target_assembly):
            data_out['assembly_target'].append(os.path.join(target_assembly, item))

    return data_out


def outgroup_dt_info(outgroup_assembly_dir, outgroup_reads_dir):
    data_out = {}
    if outgroup_assembly_dir:
        data_out['assembly_outgroup'] = []
        outgroup_assembly_dir = os.path.abspath(outgroup_assembly_dir)
        for item in os.listdir(outgroup_assembly_dir):
            data_out['assembly_outgroup'].append(os.path.join(outgroup_assembly_dir, item))
    elif outgroup_reads_dir:
        data_out['reads_outgroup'] = {}
        outgroup_reads_dir = os.path.abspath(outgroup_reads_dir)
        for item in os.listdir(outgroup_reads_dir):
            data_out['reads_outgroup'][item] = []
            for ff in os.listdir(os.path.join(outgroup_reads_dir,item)):
                data_out['reads_outgroup'][item].append(os.path.join(outgroup_reads_dir, item, ff))

    return data_out


def index_ref_genome_BWA(ref_genome, out_prefix):
    '''
    using "bwa index" to index refrence genome for BWA mapping reads. command using like:
    `bwa index -a is -p ref_genome ref_genome`
    '''
    # Before perform mapping whit BWA, FM-index should been carry out.
    cmdname = 'bwa'
    methods = 'index'
    opt_p = '-p ' + out_prefix #prefix of results
    opt_a = '-a is' # Note, IS algoritmn suit for genome not longer than 2Gbp.
    para_genome = ref_genome
    cmd = ' '.join([cmdname, methods, opt_p, opt_a, para_genome])
    print(cmd)
    os.system(cmd)
    ref_genome_indexed_bwa = ref_genome
    #ref_genome_indexed is a refrence genome indexed prefix with path.
    return ref_genome_indexed_bwa


def map_reads_BWA(ref_genome_indexed_bwa, read1, read2, num_threads, out_file):
    '''
    Map pair end reads into ref_genome using BWA. command using like:
    'bwa mem -t num_threads ref_genome read1.fq read2.fq'
    '''
    cmdname = 'bwa'
    methods = 'mem'
    opt_t = '-t ' + str(num_threads)
    opt_o = '-o ' + out_file
    para_ref_genome_prefix = ref_genome_indexed_bwa
    para_read1 = read1
    para_read2 = read2
    cmd = ' '.join([cmdname, methods, opt_t, opt_o, para_ref_genome_prefix, \
        para_read1, para_read2])
    print(cmd)
    os.system(cmd)
    #bwa_map_res_path is bwa mem result file with path.
    return out_file


def index_ref_genome_samtools(ref_genome):
    '''
    using "samtools index" to index refrence genome for bcftools.
    `samtools faidx ref_genome`
    '''
    cmdname = 'samtools'
    methods = 'faidx'
    para_ref_genome = ref_genome
    cmd = ' '.join([cmdname, methods, para_ref_genome])
    print(cmd)
    os.system(cmd)
    ref_genome_indexed_samtools = ref_genome
    #ref_genome_indexed_samtools is prefix of refrence genome with path after indexed.
    return ref_genome_indexed_samtools


def sort_sam_samtools(sam_file, num_threads, out_file):
    '''
    Using samtools to sort BWA reads mapping result, The command using like:
    `samtools sort -O BAM -@ num_threads -o in.bam_sorted in.bam`
    '''
    cmdname = 'samtools'
    methods = 'sort'
    opt_O = '-O BAM'
    opt_threads = '-@ ' + str(num_threads)
    opt_o = '-o ' + out_file
    para_sam_file = sam_file
    cmd = ' '.join([cmdname, methods, opt_threads, opt_O, opt_o, para_sam_file])
    print(cmd)
    os.system(cmd)
    # sam file sorted result, it is a bam file with path.
    return out_file


def index_bam_samtools(sorted_bam_file, num_threads):
    '''
    Using "samtools index" to index bam file. Command using like:
    `samtools index in.bam`
    '''
    cmdname = 'samtools'
    methods = 'index'
    opt_threads = '-@ ' + str(num_threads)
    para_bam_file = sorted_bam_file
    cmd = ' '.join([cmdname, methods, opt_threads, para_bam_file])
    print(cmd)
    os.system(cmd)
    indexed_bam_samtools = para_bam_file + '.bai'
    return indexed_bam_samtools


def mpileup_bam_bcftools(ref_genome, bam_file_sorted_name, num_threads, out_file):
    '''
    Using "bcftools mpileup" to calculate base coverage on refrence genome. command
    using like:
    `bcftools mpileup -f ref_genome -o mpileup.vcf in.bam`
    '''
    cmdname = 'bcftools'
    methods = 'mpileup'
    opt_f = '-f ' + ref_genome
    opt_o = '-o ' + out_file
    opt_threads = '--threads ' + num_threads
    para_bam_file = bam_file_sorted_name
    cmd = ' '.join([cmdname, methods, opt_threads, opt_f, opt_o, para_bam_file])
    print(cmd)
    os.system(cmd)
    return out_file


def call_snp_bcftools(bcftools_mpileup_file, out_file):
    '''
    Using "bcftools call" to call SNPs. The command using like:
    `bcftools call -mv --ploidy 1 in1.bam ...`
    '''
    cmdname = 'bcftools'
    methods = 'call'
    opt_mv = '-mv'
    opt_ploidy = '--ploidy 1'
    opt_o = '-o ' + out_file
    para_bam_file = bcftools_mpileup_file
    cmd = ' '.join([cmdname, methods, opt_mv, opt_ploidy, opt_o, para_bam_file])
    print(cmd)
    os.system(cmd)
    return out_file


def add_ref_genome_data(results_container, ref_genome):
    return 0


def parse_target_reads_snp_data(mpileup_file, snp_file):
    def struc_mpileup_file(mpileup_file):
        data_out = {}
        a_line = [line.rstrip().split('\t') for line in open(mpileup_file) if line[0] != '#']
        for line in a_line:
            contig = line[0]
            position = int(line[1])
            ref_base = line[3]
            dp = line[7].split(';')[0]
            if dp == 'INDEL':
                continue
            dp = int(dp.split('=')[1])
            if contig not in data_out:
                data_out[contig] = {}
            data_out[contig][position] = [position, ref_base, dp]
        return data_out

    def struc_snp_file(snp_file):
        data_out = {}
        a_line = [line.rstrip().split('\t') for line in open(snp_file) if line[0] != '#']
        for line in a_line:
            contig = line[0]
            position = int(line[1])
            alt_base = line[4]
            qual = round(float(line[5]), 3)
            dp = line[7].split(';')[0]
            if dp == 'INDEL':
                continue
            if contig not in data_out:
                data_out[contig] = {}
            data_out[contig][position] = [position, alt_base, qual]
        return data_out

    def merge_mpileup_snp_dt(mpileup_dt, snp_dt):
        data_out = {}
        for contig in mpileup_dt:
            data_out[contig] = {}
            for position in mpileup_dt[contig]:
                dt1 = mpileup_dt[contig][position]
                if position in snp_dt[contig]:
                    dt2 = snp_dt[contig][position]
                    dt_new = dt1 + dt2[1:]
                else:
                    dt_new = dt1 + [None, None]
                data_out[contig][position] = dt_new
        return data_out

    mpileup_dt = struc_mpileup_file(mpileup_file)
    snp_dt = struc_snp_file(snp_file)
    mpileup_snp_dt = merge_mpileup_snp_dt(mpileup_dt, snp_dt)
    return mpileup_snp_dt




def main():
    ref_genome, target_reads_dir, target_assembly_dir, outgroup_assembly_dir, outgroup_reads_dir, \
        num_threads, snp_qual_cutoff, snp_coverage_cutoff, snp_diversity_cutoff = get_args()

    target_infor = target_dt_info(target_reads_dir, target_assembly_dir)
    print(target_infor)
    outgroup_infor = outgroup_dt_info(outgroup_assembly_dir, outgroup_reads_dir)
    print(outgroup_infor)

    res_dir = 'Result_' + time.strftime('%Y%m%d%H%M%S')
    tmp_dir = os.path.join(res_dir, 'tmp')
    ref_genome_dir = os.path.join(res_dir, 'ref_genome')
    snp_targe_dir = os.path.join(res_dir, 'snp_target')
    snp_outgroup_dir = os.path.join(res_dir, 'snp_outgroup')
    snp_all_dir = os.path.join(res_dir, 'snp_all')
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)
        os.mkdir(tmp_dir)
        os.mkdir(ref_genome_dir)
        os.mkdir(snp_targe_dir)
        os.mkdir(snp_outgroup_dir)
        os.mkdir(snp_all_dir)
    shutil.copy(ref_genome, ref_genome_dir)
    ref_genome = os.path.join(ref_genome_dir, os.path.basename(ref_genome))
    print(ref_genome)
    results_container = {'target_snp':{}, 'outgroup_snp':{}, 'ref_genome':{}}
    add_ref_genome_data(results_container, ref_genome)

    if list(target_infor.keys())[0] == 'reads_target':
        out_prefix = ref_genome
        ref_genome_indexed_bwa = index_ref_genome_BWA(ref_genome, out_prefix)
        ref_genome_indexed_samtools = index_ref_genome_samtools(ref_genome)
        for strain in target_infor['reads_target']:
            read1, read2 = target_infor['reads_target'][strain]
            out_file = os.path.join(tmp_dir, 'bwa_mem.tmp.sam')
            bwa_map_res_path = map_reads_BWA(ref_genome_indexed_bwa, read1, read2, num_threads, out_file)
            out_file = os.path.join(tmp_dir, 'samtools_sorted.bam')
            bam_file_sorted = sort_sam_samtools(bwa_map_res_path , num_threads, out_file)
            indexed_bam_samtools = index_bam_samtools(bam_file_sorted, num_threads)
            out_file = os.path.join(tmp_dir, 'mpileup_bcftools.vcf')
            mpileup_file = mpileup_bam_bcftools(ref_genome, bam_file_sorted, out_file)
            out_file = os.path.join(tmp_dir, 'bcf_call_snp.vcf')
            snp_file = call_snp_bcftools(mpileup_file, out_file)
            parse_target_reads_snp_data(results_container, mpileup_file, snp_file)

#    elif target_infor.keys()[0] == 'assembly_target':
#        print('utill now, only reads target data is supported.')

    return 0


if __name__ == '__main__':
    main()
