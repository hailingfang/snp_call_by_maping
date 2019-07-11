#! /usr/bin/env python3

import argparse
import os

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


def index_ref_genome_BWA(ref_genome):
    '''
    using "bwa index" to index refrence genome for BWA mapping reads. command using like:
    `bwa index -a is -p ref_genome ref_genome`
    '''
    # Before perform mapping whit BWA, FM-index should been carry out.
    cmdname = 'bwa'
    methods = 'index'
    opt_p = '-p ' + ref_genome
    opt_a = '-a is' # Note, IS algoritmn suit for genome not longer than 2Gbp.
    para_genome = ref_genome
    cmd = ' '.join([cmdname, methods, opt_p, opt_a, para_genome])
    print(cmd)
    os.system(cmd)
    ref_genome_indexed_bwa = ref_genome
    #ref_genome_indexed is a refrence genome indexed prefix with path.
    return ref_genome_indexed_bwa


def map_reads_BWA(ref_genome_indexed_bwa, read1, read2, num_threads):
    '''
    Map pair end reads into ref_genome using BWA. command using like:
    'bwa mem -t num_threads ref_genome read1.fq read2.fq'
    '''
    cmdname = 'bwa'
    methods = 'mem'
    opt_t = '-t ' + str(num_threads)
    bwa_map_res_path = os.path.join(os.path.dirname(read1), \
        os.path.basename(os.path.dirname(read1)) + '_bwa_mem_out.sam')
    opt_o = '-o ' + bwa_map_res_path
    para_ref_genome_prefix = ref_genome_indexed_bwa
    para_read1 = read1
    para_read2 = read2
    cmd = ' '.join([cmdname, methods, opt_t, opt_o, para_ref_genome_prefix, \
        para_read1, para_read2])
    print(cmd)
    os.system(cmd)
    #bwa_map_res_path is bwa mem result file with path.
    return bwa_map_res_path


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


def sort_sam_samtools(sam_file, num_threads):
    '''
    Using samtools to sort BWA reads mapping result, The command using like:
    `samtools sort -O BAM -@ num_threads -o in.bam_sorted in.bam`
    '''
    cmdname = 'samtools'
    methods = 'sort'
    opt_O = '-O BAM'
    opt_threads = '-@ ' + str(num_threads)
    bam_file_name = sam_file[:-4] + '_samtools_sorted.bam'
    opt_o = '-o ' + bam_file_name
    para_sam_file = sam_file
    cmd = ' '.join([cmdname, methods, opt_threads, opt_O, opt_o, para_sam_file])
    print(cmd)
    os.system(cmd)
    # sam file sorted result, it is a bam file with path.
    return bam_file_name


def index_bam_samtools():
    '''
    Using "samtools index" to index bam file. Command using like:
    `samtools index in.bam`
    '''
    pass



def call_snp_bcftools():
    '''
    Using "bcftools call" to call SNPs. The command using like:
    `bcftools mpilepy -f ref_genome in.bam`
    `bcftools call -mv -f ref_genome --ploidy 1 in1.bam ...`
    '''
    pass


def main():
    ref_genome, target_reads_dir, target_assembly_dir, outgroup_assembly_dir, outgroup_reads_dir, \
        num_threads, snp_qual_cutoff, snp_coverage_cutoff, snp_diversity_cutoff = get_args()
    target_infor = target_dt_info(target_reads_dir, target_assembly_dir)
    print(target_infor)
    outgroup_infor = outgroup_dt_info(outgroup_assembly_dir, outgroup_reads_dir)
    print(outgroup_infor)
    results_container = {'target_snp':{}, 'outgroup_snp':{}, 'ref_genome':{}}

    if list(target_infor.keys())[0] == 'reads_target':
        #ref_genome_indexed_bwa = index_ref_genome_BWA(ref_genome)
        #ref_genome_indexed_samtools = index_ref_genome_samtools(ref_genome)
        for strain in target_infor['reads_target']:
            read1, read2 = target_infor['reads_target'][strain]
            #bwa_map_res_path = map_reads_BWA(ref_genome_indexed_bwa, read1, read2, num_threads)
            #bam_file_sorted_name = sort_sam_samtools('/home/fanghl/git_work_space/snp_call_by_maping/testdata/target_genome_reads/CL100075791_L02_565/CL100075791_L02_565_bwa_mem_out.sam' , num_threads)

    elif target_infor.keys()[0] == 'assembly_target':
        print('utill now, only reads target data is supported.')

    return 0


if __name__ == '__main__':
    main()
