#! /usr/bin/env python3

import argparse
import os
import time
import shutil
import fbio.fparse
import pickle

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
    args.add_argument('-snp_qual_cutoff', type=int, default=40, help='QUAL value cutoff of vcf result, default is 40.')
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
    opt_threads = '--threads ' + str(num_threads)
    para_bam_file = bam_file_sorted_name
    cmd = ' '.join([cmdname, methods, opt_threads, opt_f, opt_o, para_bam_file])
    print(cmd)
    os.system(cmd)
    return out_file


def call_snp_bcftools(bcftools_mpileup_file, num_threads, out_file):
    '''
    Using "bcftools call" to call SNPs. The command using like:
    `bcftools call -mv --ploidy 1 in1.bam ...`
    '''
    cmdname = 'bcftools'
    methods = 'call'
    opt_mv = '-mv'
    opt_ploidy = '--ploidy 1'
    opt_num_threads = '--threads ' + str(num_threads)
    opt_o = '-o ' + out_file
    para_bam_file = bcftools_mpileup_file
    cmd = ' '.join([cmdname, methods, opt_mv, opt_ploidy, opt_num_threads, opt_o, para_bam_file])
    print(cmd)
    os.system(cmd)
    return out_file


def parse_ref_genome_data(ref_genome):
    dt_out = {}
    dt = fbio.fparse.Fasta_parse(ref_genome)
    dt.join_lines()
    dt = dt.data
    for head in dt:
        new_head = head.split()[0]
        dt_out[new_head] = dt[head]
    return dt_out


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


def merge_target_and_outgroup_dt():
    return 0


def transform_snp_dt(snp_dt_raw):
    dt_transformed = {}
    check_dt = {}
    for strain in snp_dt_raw:
        for contig in snp_dt_raw[strain]:
            if contig not in dt_transformed: dt_transformed[contig] = {}
            if contig not in check_dt: check_dt[contig] = []
            check_dt[contig].append(len(snp_dt_raw[strain][contig]))
            for position in snp_dt_raw[strain][contig]:
                if position not in dt_transformed[contig]: dt_transformed[contig][position] = {}
                dt_transformed[contig][position][strain] = snp_dt_raw[strain][contig][position]

    for contig in check_dt:
        if min(check_dt[contig]) != max(check_dt[contig]):
            raise Exception('opps')
    # dt_transformed: {contig:{position:{strain:[pos, ref, dp, alt, qual]}}}
    return dt_transformed


def filter_snp_dt(snp_dt_transed, dp_cutoff, qual_cutoff, diversity_cutoff):

    def filter(column, dp_cutoff, qual_cutoff, diversity_cutoff):
        for ele in column:
            if ele[2] < dp_cutoff:
                return 0
            if ele[3] and ele[4] < qual_cutoff:
                return 0
        seq_line = []
        for ele in column:
            if ele[3]:
                seq_line.append(ele[3])
            else:
                seq_line.append(ele[1])
        counter = {}
        for base in seq_line:
            if base not in counter: counter[base] = 0
            counter[base] += 1
        if len(counter) > 4:
            print('Worning, More than four kinds of base found in one position...')
        if len(counter) < 2:
            return 0
        tmp = []
        for base in counter:
            tmp.append(counter[base])
        tmp.sort()
        if sum(tmp[:-1])/sum(tmp) < diversity_cutoff:
            return 0
        else:
            return 1

    ok_position = {}
    snp_out = {}
    for contig in snp_dt_transed:
        if contig not in snp_out: snp_out[contig] = {}
        if contig not in ok_position: ok_position[contig] = []
        for position in snp_dt_transed[contig]:
            column = []
            for strain in snp_dt_transed[contig][position]:
                column.append(snp_dt_transed[contig][position][strain])
            judge_snp_bool = filter(column, dp_cutoff, qual_cutoff, diversity_cutoff)
            if judge_snp_bool:
                ok_position[contig].append(position)
                for strain in snp_dt_transed[contig][position]:
                    if strain not in snp_out[contig]:
                        snp_out[contig][strain] = {}
                        snp_out[contig][strain][position] = snp_dt_transed[contig][position][strain]
    return snp_out, ok_position


def print_res(snp_filtered, ok_position, ref_genome, out_dir):
    snp_res_detail = open(os.path.join(out_dir, 'snp_detail'), 'w')
    snp_seq = open(os.path.join(out_dir, 'snp_seq.fasta'), 'w')
    all_snp_seq = {'ref_pos':[], 'ref_base':[]}
    for contig in ok_position:
        print('>', contig, file=snp_filtered)
        positions = ok_position[contig]
        positions.sort()
        tmp_ref_pos = []
        tmp_ref_base = []
        for pos in positions:
            tmp_ref_pos.append(str(pos))
            tmp_ref_base.append(ref_genome[contig][pos - 1])
        all_snp_seq['ref_pos'] += tmp_ref_pos
        all_snp_seq['ref_base'] += tmp_ref_base
        print('ref_pos', ''.join(tmp_ref_pos), sep='\t', file=snp_res_detail)
        print('ref_base', ''.join(tmp_ref_base), sep='\t', file=snp_res_detail)

        for strain in snp_filtered[contig]:
            tmp_base = []
            for pos in positions:
                if snp_filtered[contig][strain][pos][3]:
                    tmp_base.append(snp_filtered[contig][strain][pos][3])
                else:
                    tmp_base.append(snp_filtered[contig][strain][pos][1])
            print(strain, ''.join(tmp_base), sep='\t', file=snp_res_detail)
            if strain not in all_snp_seq:
                all_snp_seq[strain] = []
            all_snp_seq[strain] += tmp_base

    print('>', 'ref_pos', file=snp_seq)
    print(''.join(all_snp_seq['ref_pos']), file=snp_seq)
    print('>', 'ref_base', file=snp_seq)
    print(''.join(all_snp_seq['ref_base']), file=snp_seq)
    all_snp_seq.pop('ref_pos')
    all_snp_seq.pop('ref_base')
    for strain in all_snp_seq:
        print('>', strain, file=snp_seq)
        print(''.join(all_snp_seq[strain]), fi松露le=snp_seq)
    return 0

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
    ref_genome_dt = parse_ref_genome_data(ref_genome)
    results_container['ref_genome'] = ref_genome_dt

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
            mpileup_file = mpileup_bam_bcftools(ref_genome, bam_file_sorted, num_threads, out_file)
            out_file = os.path.join(tmp_dir, 'bcf_call_snp.vcf')
            snp_file = call_snp_bcftools(mpileup_file, num_threads, out_file)
            mpileup_snp_dt = parse_target_reads_snp_data(mpileup_file, snp_file)
            # There maybe need to print some raw data for future analyse.
            #
            results_container['target_snp'][strain] = mpileup_snp_dt
            keeping_dir = os.path.join(snp_targe_dir, strain)
            os.mkdir(keeping_dir)
            shutil.move(mpileup_file, keeping_dir)
            shutil.move(snp_file, keeping_dir)
        tmp_bk = open('tmp.pickle', 'wb')
        pickle.dump(results_container['target_snp'], tmp_bk)
    elif list(target_infor.keys())[0] == 'assembly_target':
        print('utill now, only reads target data is supported.')

#    if list(outgroup_infor.keys())[0] == 'reads_outgroup':
#        pass
#    elif list(outgroup_infor.keys())[0] == 'assembly_outgroup':
#        pass
    snp_dt_raw = results_container['target_snp']
    snp_dt_transed = transform_snp_dt(snp_dt_raw)
    snp_filtered, ok_position = filter_snp_dt(snp_dt_transed, snp_coverage_cutoff, snp_qual_cutoff, snp_diversity_cutoff)
    print_res(snp_filtered, ok_position, results_container['ref_genome'], snp_all_dir)

    return 0


if __name__ == '__main__':
    main()
