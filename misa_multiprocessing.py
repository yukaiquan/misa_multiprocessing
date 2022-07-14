'''
Author: yukaiquan
Date: 2022-06-24 08:07:38
Email: 1962568272@qq.com
LastEditor: yukaiquan
Description: Genome-wide multi-process chromosome-level SSR search.
FileName: misa_multiprocessing.py
'''
import os
import shutil
import sys
import time
from multiprocessing import Pool
import pandas as pd
from subprocess import Popen, PIPE
import pybedtools
from functools import partial
import argparse


def main():
    # 1. Get the genome file
    start_time = time.time()
    args = check_options(get_options())
    num_treads = int(args.threads)
    genome_file = args.genome
    misabin = './misa.pl'
    print('Genome file: %s' % genome_file)
    print('Misa bin: %s' % misabin)
    # get file path
    file_path = os.path.dirname(os.path.abspath(__file__))
    tmp_dir = file_path + '/tmp'
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.exists(genome_file + '.fai'):
        samtools = samtools_index(genome_file)
        if not samtools:
            print('Error: samtools index error')
            sys.exit(1)
    fai_file_name = genome_file + '.fai'
    if not os.path.exists(fai_file_name):
        print('Error: %s not exists' % fai_file_name)
        sys.exit(1)
    chr_length = check_chromosome(fai_file_name)
    print(chr_length)
    if len(chr_length) == 0:
        print('Error: %s not exists' % fai_file_name)
        sys.exit(1)
    val_chrname = []
    if args.chr == 'all':
        for i in chr_length.keys():
            val_chrname.append(i)
    else:
        val_chrname = args.chr.split(',')
    ##################################################多进程染色体拆分#############################################################
    # 2. Get the chr seq
    ex_chr_id = []
    for chr_name, length in chr_length.items():
        if chr_name not in val_chrname:
            continue
        print('%s length: %d' % (chr_name, length))
        ex_chr_id.append(chr_name)
    print('test:' + genome_file)
    print('test:' + str(ex_chr_id))
    pool = Pool(processes=num_treads)
    samtools_faidx_partial = partial(samtools_getseq, genome_file)
    pool.map(samtools_faidx_partial, ex_chr_id)
    pool.close()
    pool.join()
    #################################################first_misa############################################################
    print(chr_length)
    chr_names = []
    for chr_name, length in chr_length.items():
        chr_names.append(tmp_dir + '/' + chr_name + '.fasta')
    pool = Pool(processes=num_treads)
    pool.map(misa, chr_names)
    pool.close()
    pool.join()
    print('first misa pooldown!')
    print(chr_names)
    #################################################多进程提取bed to fasta####################################################
    pool = Pool(processes=num_treads)
    partial_func = partial(
        ex_bed_fasta, chr_length=chr_length)
    pool.map(partial_func, chr_names)
    pool.close()
    pool.join()
    second_fasta = []
    for chr_name in chr_names:
        misa_genome_second = chr_name + '_second.fasta'
        second_fasta.append(misa_genome_second)
    #################################################second_misa############################################################
    pool = Pool(processes=num_treads)
    pool.map(misa, second_fasta)
    pool.close()
    pool.join()
    print('seconde misa pooldown!')
    misa_result_df = pd.DataFrame()
    for misa_re in second_fasta:
        misa_df = pd.read_csv(misa_re + '.misa', sep='\t')
        misa_result_df = pd.concat(
            [misa_result_df, misa_df], ignore_index=True)
    merge_file(genome_file, tmp_dir)
    misa_result_df.to_csv(genome_file + '_second.tsv', sep='\t', index=False)
    for second_fa in second_fasta:
        mv_file(second_fa + '.statistics',
                './' + os.path.basename(second_fa) + '.statistics')
    shutil.rmtree(tmp_dir)
    stop_time = time.time()
    print('Total time: %s' % (stop_time - start_time))
    print('All Done!')


def ex_bed_fasta(chr_names, chr_length):
    bed_file_name = chr_names + '.bed'
    misa_name = chr_names + '.misa'
    misa_result = misa_result_to_bed(
        misa_name, bed_file_name, chr_length)
    if not misa_result:
        print('Error: %s not exists' % misa_name)
        sys.exit(1)
    bed_file = pybedtools.BedTool(bed_file_name)
    misa_genome_second = chr_names + '_second.fasta'
    bed_file.sequence(fi=chr_names,
                      fo=misa_genome_second)
    return True


def rm_file(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)
    return True


def mv_file(file_name, new_file_name):
    if os.path.exists(file_name):
        shutil.move(file_name, new_file_name)
    return True


def merge_file(genome_file, tmp_dir):
    cmd = 'cat ' + tmp_dir + '/' + '*_second.fasta > ' + \
        genome_file + '_second_all.fasta'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print('Error: %s' % err)
        return False
    return True


# def dict_getseq(genome_file, chr_dict, val_chrname, tmp_dir):
#     for chr_name, chr_length in chr_dict.items():
#         if chr_name not in val_chrname:
#             continue
#         print('%s length: %d' % (chr_name, chr_length))
#         samtools_get = samtools_getseq(genome_file, chr_name, tmp_dir)
#         if not samtools_get:
#             print('Error: %s not exists' % chr_name + '.fasta')
#             sys.exit(1)
#     return True


def samtools_getseq(genome_file, ex_chr_id):
    cmd = 'samtools faidx %s %s > %s' % (
        genome_file, str(ex_chr_id),  './tmp/' + str(ex_chr_id) + '.fasta')
    print(cmd)
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print('Error: %s' % err)
        return False
    samtools = samtools_index('./tmp/' + str(ex_chr_id) + '.fasta')
    if not samtools:
        print('Error: samtools index error')
        sys.exit(1)
    return True


def check_chromosome(fai_file_name):
    chr_length = {}
    with open(fai_file_name, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            seq_id = line[0]
            seq_length = int(line[1])
            if seq_id.lower() == 'chrun':
                print('Error: %s is unknown sequence' % seq_id)
                continue
            chr_length[seq_id] = seq_length
    return chr_length


def readfasta(lines):
    seq = []
    index = []
    seqplast = ""
    numlines = 0
    for i in lines:
        if ">" in i:
            index.append(i.replace("\n", "").replace(">", ""))
            seq.append(seqplast.replace("\n", ""))
            seqplast = ""
            numlines += 1
        else:
            seqplast = seqplast + i.replace("\n", "")
            numlines += 1
        if numlines == len(lines):
            seq.append(seqplast.replace("\n", ""))
    seq = seq[1:]
    return index, seq


# function of split table in txt
def str_split(lines):
    list2 = lines.split()
    return list2


def misa_result_to_bed(misa_result, bed_file, chr_length):
    with open(bed_file, 'w') as f:
        with open(misa_result, 'r') as o:
            for line in o:
                if line.startswith('ID'):
                    continue
                line = line.strip().split('\t')
                seq_id = str(line[0])
                start = int(line[5]) - 150
                if start < 0:
                    start = 0
                end = int(line[6]) + 150
                if end > chr_length[seq_id]:
                    end = chr_length[seq_id]
                f.write(seq_id + '\t' + str(start) + '\t' + str(end) + '\n')
    return True


def misa_result(genome_file):
    misa_name = genome_file + '.misa'
    statistics_name = genome_file + '.statistics'
    if not os.path.exists(misa_name):
        print('Error: %s not exists' % misa_name)
        return False
    if not os.path.exists(statistics_name):
        print('Error: %s not exists' % statistics_name)
        return False
    return True


def misa(genome_file):
    print(type(genome_file))
    misacmd = ' '.join(['perl ./misa.pl', genome_file])
    print(misacmd)
    runmisa = Popen(misacmd, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = runmisa.communicate()
    if runmisa.returncode != 0:
        print('Error: %s' % err)
        return False
    return True


def samtools_index(genome_file):
    samtools_index_cmd = ' '.join(
        ['samtools faidx', genome_file])
    print(samtools_index_cmd)
    run_samtools_index = Popen(
        samtools_index_cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = run_samtools_index.communicate()
    if run_samtools_index.returncode != 0:
        print('Error: %s' % err)
        return False
    return True


def check_options(parser):

    args = parser.parse_args()
    if not args.genome:

        print("Can not locate target genome file, please input a target genome file.\n")

        parser.print_help()

        sys.exit(1)

    print("#" * 40)

    print("genome file:",
          os.path.abspath(os.path.expanduser(args.genome)))
    print("chr_list:", args.chr)
    print("threads:", args.threads)

    print("#" * 40)

    return args


def get_options():

    parser = argparse.ArgumentParser(description="whole genome ssr find and primer design. Author:yukaiquan; Email:1962568272@qq.com", prog='kq_misa.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="Example:\n"
                                            "python kq_misa.py -g SFS.fasta -t 1\\ \n"
                                     )

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument(
        '-g',
        '--genome',
        dest='genome',
        help='Fatsa genome, should include all sequences from genome file.',
        type=str,
        required=True)
    parser.add_argument(
        '-c',
        '--chr',
        dest='chr',
        help='Chromosome name, should be in the format of chr1A,chr2A,chr3A...(default is all chromosomes).',
        default='all')
    parser.add_argument(
        '-t',
        '--threads',
        dest='threads',
        help='The number of threads (default is 1).',
        default=1,
    )

    return parser


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)
