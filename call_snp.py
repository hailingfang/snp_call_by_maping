#! /usr/bin/env python3

import argparse
import os

def get_args():
    args = argparse.ArgumentParser(description='This utility is used to get SNPs \
        consitance sequence.')
    args.add_argument('-ref_genome', type=str, required=True, help='the file name \
        of refrence genome. all reads or genome of target will mapped to this genome')
    args.add_argument('-target_genome_reads', type=str, default=None, help='target reads \
        file. reads file of one item needed put into a directory and put all this item file \
        into target_genome_reads directory. Utill now, only pair end reads are implementation.')
    args.add_argument('-target_genome_assembly', type=str, default=None, help='put all target \
        genome into a single directory.')
    args.add_argument('-outgroup_genome_assembly', type=str, default=None, help='genome \
        assembly of outgroup, multiple outgroup can be added, put all assembly into a \
        single directory.')
    args.add_argument('-outgroup_genome_reads', type=str, default=None, help='genome \
        reads of outgroup. multiple genome can be used. put reads file into a directory of \
        a single item and put all item into outgroup_genome_reads directory.')
    args.add_argument('-clean_up', type=str, default='No', choice=['No', 'Yes'], help='clean up \
        all tmperary file when runing this pipeline. "No" was set as default.')
    args = args.parse_args()
    ref_genome, target_reads_dir, target_assembly_dir, outgroup_assembly_dir, \
    outgroup_reads_dir = args.ref_genome, args.target_genome_reads, args.target_genome_assembly, \
    args.outgroup_genome_assembly, args.outgroup_genome_reads
    if (not (target_reads_dir or target_assembly_dir)) or (target_reads_dir and target_assembly_dir):
        exit('you have to input one type target data at most and least.')

    return ref_genome, target_reads_dir, target_assembly_dir, outgroup_assembly_dir, outgroup_reads_dir


def target_dt_info(target_reads, target_assembly):
    data_out = {}
    if target_reads:
        target_reads = os.path.abspath(target_reads)
        data_out['reads_target'] = {}
        for item in os.listdir(target_reads):
            data_out['reads_target'][item]=[]
            for ff in os.listdir(os.path.join(target_reads,item)):
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


def map_reads_BWA(ref_genome, read1, read2):
    # Before perform mapping whit BWA, FM-index should been carry out.
    def ref_genome_index(ref_genome):
        cmdname = 'bwa'
        methods = 'index'
        opt_p = '-p '+os.path.basename(ref_genome)
        opt_a = '-a is' # Note, IS algoritmn suit for genome not longer than 2Gbp.
        para_genome = ref_genome
        cmd = ' '.join([cmdname, methods, opt_p, opt_a, para_genome])
        print(cmd)
        # os.system(cmd)
        ref_genome_indexed = ref_genome
        return ref_genome_indexed

    def run_method_mem(ref_genome_indexed, read1, read2):
        cmdname = 'bwa'
        methods = 'mem'
        opt_t = '-t 8'
        ref_genome_prefix = ref_genome_indexed
        para_read1 = read1
        para_read2 = read2
        para_out = '>'+os.path.join(os.path.dirname(read1),'bwa_mem_out.sam')
        cmd = ' '.join([cmdname, methods, opt_t, ref_genome_prefix, para_read1, para_read2, para_out])
        print(cmd)
        #os.system(cmd)
        return 0

    return 0



def main():
    ref_genome, target_reads_dir, target_assembly_dir, outgroup_assembly_dir, outgroup_reads_dir \
        = get_args()
    target_infor = target_dt_info(target_reads_dir, target_assembly_dir)
    outgroup_infor = outgroup_dt_info(outgroup_assembly_dir, outgroup_reads_dir)


    return 0


if __name__ == '__main__':
    main()
