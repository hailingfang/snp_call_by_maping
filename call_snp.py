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
        into target_genome_reads directory.')
    args.add_argument('-target_genome_assembly', type=str, default=None, help='put all target \
        genome into a single directory.')
    args.add_argument('-outgroup_genome_assembly', type=str, default=None, help='genome \
        assembly of outgroup, multiple outgroup can be added, put all assembly into a \
        single directory.')
    args.add_argument('-outgroup_genome_reads', type=str, default=None, help='genome \
        reads of outgroup. multiple genome can be used. put reads file into a directory of \
        a single item and put all item into outgroup_genome_reads directory.')
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


def main():
    ref_genome, target_reads_dir, target_assembly_dir, outgroup_reads_dir, outgroup_assembly_dir \
        = get_args()

    return 0


if __name__ == '__main__':
    main()
