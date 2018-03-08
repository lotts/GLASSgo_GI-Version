#!/usr/bin/env python3

import argparse
import re
import time
import os
import shlex
import subprocess
import sys
from subprocess import Popen, PIPE
from Bio import SeqIO
from Bio.Application import ApplicationError
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import pairwise2
from math import sqrt, pow, log10

"""
GLASSgo Version 1.5.0


MIT License

Copyright (c) 2017 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

def tmp_counter():
    if 'cnt' not in tmp_counter.__dict__:
        tmp_counter.cnt = 0
    tmp_counter.cnt += 1
    return tmp_counter.cnt


def compute_new_min_identity(in_seq_len):
    new_min_identity = 60 - pow((0.138 * sqrt(pow(in_seq_len, 0.95))), 2.8)
    if new_min_identity < 55:
        new_min_identity = 55
    elif new_min_identity > 62:
        new_min_identity = 62
    return round(new_min_identity, 3)


def analyze_input_sequence(in_seek_seq):
    f_fasta = open(in_seek_seq, 'rU')
    seek_id = ""
    seek_seq = ""
    seek_seq_len = 0
    count = 0
    for record in SeqIO.parse(f_fasta, "fasta"):
        seek_id = record.description
        seek_seq = str(record.seq)
        seek_seq_len = len(seek_seq)
        count += 1
    f_fasta.close()
    # to many sequences in the input file
    if count > 1:
        return "false"
    # create tmp_file with modified input seq
    file_name = './id.%s.tmp' % os.getpid()
    file_name = str(file_name) + str(tmp_counter())
    # check sequence (only 'A','T','C','G' are allowed)
    seek_seq = seek_seq.upper()
    seek_seq = seek_seq.replace("U", "T")
    regex = r'^[A,T,C,G]*$'
    if re.match(regex, seek_seq):
        # write to file
        f = open(file_name, 'w')
        f.write(">" + str(seek_id) + "\n")
        f.write(str(seek_seq) + "\n")
        f.close()
        return "true", file_name, seek_seq_len
    else:
        return "false", "", seek_seq_len


def plus_or_minus_strand(s_pos, e_pos):
    if s_pos < e_pos:
        ori = "'plus'"
    else:
        ori = "'minus'"
    return ori


def rotate_range(in_range1, in_range2):
    if int(in_range1) < int(in_range2):
        return int(in_range1), int(in_range2), "n"
    else:
        return int(in_range2), int(in_range1), "c"


def get_sequence_range(s_pos, e_pos):
    if s_pos < e_pos:
        out = str(s_pos) + "-" + str(e_pos)
        return out
    else:
        out = str(e_pos) + "-" + str(s_pos)
        return out


def check_seq_blast_length(s_pos, e_pos, seek_length):
    left = s_pos - 1
    right = seek_length - e_pos
    if right <= 0:
        right = 1
    if left <= 0:
        left = 1
    return left, right


def compress_blast_12cols(blast_in):
    id_index = dict()
    blast = blast_in.split("\n")
    for line in blast:
        line = line.rstrip()
        if line:
            line_arr = line.split("\t")
            if line_arr[1] in id_index:
                found = "false"
                (act_val1, act_val2, act_flag) = rotate_range(line_arr[4], line_arr[5])
                counter = 0
                for entry in id_index[line_arr[1]]:
                    (stored_val1, stored_val2, stored_flag) = rotate_range(entry[4], entry[5])
                    if (act_flag == stored_flag) and (stored_val1 <= act_val1 <= stored_val2 or \
                       stored_val1 <= act_val2 <= stored_val2 or \
                       act_val1 <= stored_val1 and stored_val2 <= act_val2):
                        # update values
                        if float(line_arr[6]) < float(entry[6]):
                            id_index[line_arr[1]][counter][0] = line_arr[0]
                            id_index[line_arr[1]][counter][1] = line_arr[1]
                            id_index[line_arr[1]][counter][2] = line_arr[2]
                            id_index[line_arr[1]][counter][3] = line_arr[3]
                            id_index[line_arr[1]][counter][4] = line_arr[4]
                            id_index[line_arr[1]][counter][5] = line_arr[5]
                            id_index[line_arr[1]][counter][6] = line_arr[6]
                            id_index[line_arr[1]][counter][7] = line_arr[7]
                        found = "true"
                    counter += 1
                if found == "false":
                    (id_index[line_arr[1]]).append(line_arr)
            else:
                id_index[line_arr[1]] = list()
                (id_index[line_arr[1]]).append(line_arr)
    # produce new cured BLAST out
    blast_out = ""
    for key in id_index:
        for entry in id_index[key]:
            blast_out += "\t".join(entry[0:]) + "\n"
    return blast_out


def extent_reads_from_blast_table(in_seek_seq, blast, db, e_value):
    # read files
    f_fasta = open(in_seek_seq, 'rU')
    for record in SeqIO.parse(f_fasta, "fasta"):
        seek_seq = record.seq
        seek_length = len(seek_seq)

    # compress blast result (remove and merge same results to speed up)
    blast_cured = compress_blast_12cols(blast)

    array = list()
    blast = blast_cured.split("\n")
    for line in blast:
        line = line.rstrip()
        if line:
            line = line.split("\t")
            # check e-value and compare it
            if float(line[6]) <= e_value and line[1].startswith("gi"):
                gi = -1
                gi = int(re.findall("gi\|(\d+)\|", line[1])[0])
                orientation = plus_or_minus_strand(int(line[4]), int(line[5]))
                if orientation == "'plus'":
                    (left, right) = check_seq_blast_length(int(line[2]), int(line[3]), seek_length)
                    s_pos_genome = int(line[4]) - left + 1
                    if s_pos_genome <= 0:
                        s_pos_genome = 1
                    e_pos_genome = int(line[5]) + right - 1
                else:
                    (left, right) = check_seq_blast_length(int(line[2]), int(line[3]), seek_length)
                    s_pos_genome = int(line[4]) + left - 1
                    e_pos_genome = int(line[5]) - right + 1
                    if e_pos_genome <= 0:
                        e_pos_genome = 1

                seq_range = get_sequence_range(int(s_pos_genome), int(e_pos_genome))
                act_args = shlex.split("blastdbcmd -db " + str(db) +
                                " -dbtype 'nucl' -entry " + str(gi) +
                                " -range " + str(seq_range) + " -strand " + str(orientation))

                p = subprocess.Popen(act_args, stdout=PIPE)
                (output, err) = p.communicate()
                p.wait()
                lines = output.decode('utf-8').split("\n")

                if len(lines) >= 2:
                    coord1 = re.findall(":[a-z]*(\d*)-", lines[0])[0]
                    coord2 = re.findall("\d*-(\d*)\s+", lines[0])[0]
                    id_taxID = str(lines[0]) + "-taxID:" + str(line[7])

                    if coord1 < coord2:
                        ori = "+"
                        # id, sequence, coord1, coord2, ori, identity
                        array.append([id_taxID, "".join(lines[1:]), coord1, coord2, ori, 0])
                    else:
                        ori = "-"
                        # id, sequence, coord1, coord2, ori, identity
                        array.append([id_taxID, "".join(lines[1:]), coord2, coord1, ori, 0])
    f_fasta.close()
    return array


def global_alignment(in_seek_seq, array_fasta):
    f_fasta = open(in_seek_seq, 'rU')
    seek_seq = ""
    for record in SeqIO.parse(f_fasta, "fasta"):
        seek_seq = record.seq
    f_fasta.close()

    for entry in array_fasta:
        aln_1, aln_2, score, begin, end = pairwise2.align.globalms(str(seek_seq), str(entry[1]), 2, -1, -.5, -.1)[0]
        matches = sum(aa1 == aa2 for aa1, aa2 in zip(str(aln_1), str(aln_2)))
        pct_identity = 100.0 * matches / len(str(aln_1))
        entry[5] = round(pct_identity, 2)
    return array_fasta


def analyze_sequences(fasta_hits_rated, min_identity, max_identity, protected_mode, range1, range2, range3):
    directly_stored_lum = list()
    directly_stored_high_sim = list()
    only_significant_fasta_hits = list()
    range1_arr = range1.split("-")
    range2_arr = range2.split("-")
    range3_arr = range3.split("-")
    count_range1 = 0
    count_range2 = 0
    count_range3 = 0

    for entry in sorted(fasta_hits_rated, key=lambda x: x[5]):
        if 100 >= entry[5] >= max_identity:
            directly_stored_high_sim.append(entry)
        if max_identity > entry[5] >= min_identity:
            directly_stored_lum.append(entry)
        if 100 > entry[5] >= min_identity:
            # protected mode
            if protected_mode == "on":
                # range 1
                if float(range1_arr[0]) > entry[5] >= float(range1_arr[1]) and count_range1 <= int(range1_arr[2]):
                    count_range1 += 1
                    only_significant_fasta_hits.append(entry)
                elif float(range2_arr[0]) > entry[5] >= float(range2_arr[1]) and count_range2 <= int(range2_arr[2]):
                    # range 2
                    count_range2 += 1
                    only_significant_fasta_hits.append(entry)
                elif float(range3_arr[0]) > entry[5] >= float(range3_arr[1]) and count_range3 <= int(range3_arr[2]):
                    # range 3
                    count_range3 += 1
                    only_significant_fasta_hits.append(entry)
            else:
                only_significant_fasta_hits.append(entry)

    return directly_stored_high_sim, directly_stored_lum, only_significant_fasta_hits


def blast_seed(in_seek_seq, db, e_value, gi_list, num_hits, n_threads):
    # compile the call
    format_str = "\"6 qseqid sseqid qstart qend sstart send evalue staxids\""
    if gi_list is None:
        blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn', query=in_seek_seq,
                                             db=db, evalue=e_value, outfmt=format_str, max_target_seqs=num_hits, num_threads=n_threads)
    else:
        blastn_cline = NcbiblastnCommandline(cmd='blastn', task='blastn', query=in_seek_seq,
                                             db=db, evalue=e_value, outfmt=format_str, gilist=gi_list, max_target_seqs=num_hits, num_threads=n_threads)
    # run call, but check for exceptions
    try:
        stdout, stderr = blastn_cline()
    except ApplicationError as blastError :
        # print error information
        sys.stderr.write( "\n\nERROR for blast call '" + str(blastn_cline) +
                          "':\n\treturn code = " + str(blastError.returncode) + "\n\terror cmd = " + str(blastError.cmd) +
                          "\n\terror stderr = " + str(blastError.stderr) + "\n\terror stdout = " + str(blastError.stdout) + "\n\n")
        # exit the tool
        sys.exit(-1)
    return stdout.rstrip()


def create_fasta_file(for_reblast):
    file_name = './id.%s.tmp' % os.getpid()
    file_name = str(file_name) + str(tmp_counter())
    fasta_w = open(file_name, "w")
    for entry in for_reblast:
        fasta_w.write(entry[0] + "\n" + entry[1] + "\n")
    fasta_w.close()
    return file_name


def create_input_for_londen(in_seek_seq, in_directly_stored, in_directly_stored_cure, trustable_hits, num_threads):
    compression_hash = dict()
    londen_scaling_mode = "relaxed"
    min_pc_allowed = 85
    max_trustable_seeds = 10
    thresh_perc = 60

    # sort output for reproducibility
    record = list()
    for key in in_directly_stored_cure:
        seq_pos_arr = in_directly_stored_cure[key].split(",")
        for i in seq_pos_arr:
            tmp_list = list()
            tmp_list.append(float(in_directly_stored[int(i)][5]))
            tmp_out_str = ">" + str(i) + "\n" + str(in_directly_stored[int(i)][1]) + "\n"
            tmp_list.append(tmp_out_str)
            record.append(tmp_list)
    # start sorting
    record_sorted = sorted(record, key=lambda x: x[0], reverse=False)

    # create input for clustal omega - londen hits
    file_name_in = "./IDo" + str(time.time()) + ".tmp"
    fasta_w = open(file_name_in, "w")
    counter_min2 = 0
    for i in range(0, len(record_sorted)):
        fasta_w.write(record_sorted[i][1])
        counter_min2 += 1
    fasta_w.close()

    if counter_min2 >= 3:
        # call clustal omega
        file_name_out = "./IDX" + str(time.time()) + ".tmp"
        clustalomega_cline = ClustalOmegaCommandline(infile=file_name_in, percentid=True,
                                                     outputorder="input-order", threads=num_threads,
                                                     max_hmm_iterations=-1,
                                                     distmat_full=True, distmat_out=file_name_out)
        # run call, but check for exceptions
        try:
            stdout, stderr = clustalomega_cline()
        except ApplicationError as error:
            sys.stderr.write("\n\nERROR for clustal omega call " + str(clustalomega_cline) + "\n\n")
            sys.exit(-1)

        handle = open(file_name_out, "r")
        matrix_size = 0
        row_lookup = list()
        matrix = list()
        count = 0
        for line in handle:
            line = line.rstrip()
            if count == 0:
                count += 1
                matrix_size = int(line)
            else:
                line = re.sub(' +',' ', line)
                line_arr = line.split(" ")
                row_lookup.append(line_arr[0])
                matrix.append(line_arr[1:])
        handle.close()
        os.system("rm " + str(file_name_out))

        # compress data
        # iterate through matrix
        forbidden_hash = dict()
        for i in range(0, matrix_size):
            nothing_clustered = "true"
            for j in range(i + 1, matrix_size):
                if i not in forbidden_hash and j not in forbidden_hash:
                    if float(matrix[i][j]) >= min_pc_allowed:
                        forbidden_hash[j] = 0
                        nothing_clustered = "false"
                        if row_lookup[i] in compression_hash:
                            compression_hash[row_lookup[i]].append(row_lookup[j])
                        else:
                            compression_hash[row_lookup[i]] = list()
                            compression_hash[row_lookup[i]].append(row_lookup[i])
                            compression_hash[row_lookup[i]].append(row_lookup[j])
            if nothing_clustered == "true" and i not in forbidden_hash:
                compression_hash[row_lookup[i]] = list()
                compression_hash[row_lookup[i]].append(row_lookup[i])
            forbidden_hash[i] = 0
    os.system("rm " + str(file_name_in))

    # store data to file - compressed Londen data + seed + max. 10 seqs of trustable reads (s0 (query),s1, ... , s9)
    # trustable_hits are needed!
    f_fasta = open(in_seek_seq, 'rU')
    seek_seq = ""
    for record in SeqIO.parse(f_fasta, "fasta"):
        seek_seq = record.seq
    f_fasta.close()

    file_name = "./IDFF" + str(time.time()) + ".tmp"
    fasta_w = open(file_name, "w")
    fasta_w.write(">s0" + "\n" + str(seek_seq) + "\n")

    # add up to x trustable seeds
    count = 0
    seq_duplication = dict()
    for i in range(0, len(trustable_hits)):
        if count == max_trustable_seeds:
            break
        if trustable_hits[i][5] < 100 and trustable_hits[i][1] not in seq_duplication:
            fasta_w.write(">s" + str(i + 1) + "\n" + str(trustable_hits[i][1]) + "\n")
            seq_duplication[trustable_hits[i][1]] = 0
            count += 1

    count_u60 = 0
    count_l60 = 0
    if counter_min2 >= 3:
        counter = 0
        for k in compression_hash:
            if in_directly_stored[int(k)][5] >= thresh_perc:
                count_u60 += 1
            else:
                count_l60 += 1
            counter += 1
            fasta_w.write(">" + str(k) + "\n" + str(in_directly_stored[int(k)][1]) + "\n")
        fasta_w.close()
    elif londen_scaling_mode == "relaxed":
        counter = 0
        for key in in_directly_stored_cure:
            seq_pos_arr = in_directly_stored_cure[key].split(",")
            for i in seq_pos_arr:
                if i in compression_hash:
                    compression_hash[i].append(i)
                else:
                    compression_hash[i] = list()
                    compression_hash[i].append(i)
                fasta_w.write(">" + str(i) + "\n" + str(in_directly_stored[int(i)][1]) + "\n")
                counter += 1
    else:
        counter = 0
        for key in in_directly_stored_cure:
            seq_pos_arr = in_directly_stored_cure[key].split(",")
            for i in seq_pos_arr:
                if in_directly_stored[int(i)][5] >= 56:
                    if i in compression_hash:
                        compression_hash[i].append(i)
                    else:
                        compression_hash[i] = list()
                        compression_hash[i].append(i)
                    fasta_w.write(">" + str(i) + "\n" + str(in_directly_stored[int(i)][1]) + "\n")
                    counter += 1
    if count_l60 == 0:
        londen_scaling_mode = 2
    else:
        londen_scaling_mode = (count_u60 / count_l60)
    return file_name, counter, compression_hash, londen_scaling_mode


def cure_results(in_directly_stored):
    id_lookup = dict()
    index = 0
    for entry in in_directly_stored:
        tmp_id = str((re.findall(">(.*?):", entry[0])[0])) + str(entry[4])
        if tmp_id in id_lookup:
            index_arr = id_lookup[tmp_id].split(",")
            found = 0
            pos = 0
            for i in index_arr:
                i = int(i)
                if int(in_directly_stored[i][2]) <= int(entry[2]) <= int(in_directly_stored[i][3]):
                    found = 1
                elif int(in_directly_stored[i][2]) <= int(entry[3]) <= int(in_directly_stored[i][3]):
                    found = 1
                elif int(entry[2]) <= int(in_directly_stored[i][2]) and int(in_directly_stored[i][3]) <= int(entry[3]):
                    found = 1

                if int(in_directly_stored[i][2]) >= int(entry[2]) >= int(in_directly_stored[i][3]):
                    found = 1
                elif int(in_directly_stored[i][2]) >= int(entry[3]) >= int(in_directly_stored[i][3]):
                    found = 1
                elif int(entry[2]) >= int(in_directly_stored[i][2]) and int(in_directly_stored[i][3]) >= int(entry[3]):
                    found = 1
                if found == 1 and float(in_directly_stored[i][5]) < float(entry[5]):
                    index_arr[pos] = str(index)
                    found = 2
                pos += 1
            if found == 0:
                index_arr.append(str(index))
                id_lookup[tmp_id] = ','.join(index_arr)
            elif found == 2:
                id_lookup[tmp_id] = ','.join(index_arr)
            #else:
            #    id_lookup[tmp_id] = str(index_arr[-1])
        else:
            id_lookup[tmp_id] = str(index)
        index += 1
    return id_lookup


def call_londen(input_fasta_file, in_directly_stored_lum, area_filter, cutting_method, londen_cl_ids_entries, londen_scaling_mode):
    # check if area filter is auto adjusted or not
    seq_len = 1000000
    if area_filter == -1:
        #for record in SeqIO.parse(query_seq, "fasta"):
        for record in SeqIO.parse(input_fasta_file, "fasta"):
            if len(record.seq) < seq_len:
                seq_len = len(record.seq)
        if londen_scaling_mode >= 1:
            cutting_method = 2
            area_filter = (0.01 * log10(seq_len)) + londen_scaling_mode
        else:
            cutting_method = 1
            area_filter = (0.29 * log10(seq_len)) + londen_scaling_mode

        if area_filter > 2.4:
            area_filter = 2.4

    scriptPath = os.path.dirname(os.path.realpath(__file__))
    act_args = shlex.split(str(scriptPath) + "/reqPackages/londen -f " + str(input_fasta_file) + " -t " + str(area_filter) + " -m " + str(cutting_method))
    p = subprocess.Popen(act_args, stdout=PIPE)
    (output, err) = p.communicate()
    p.wait()
    cluster_with_seed = output.decode('utf-8').split("\n")

    # read cluster with seed + trustable hits
    tmp_cluster_ids = list()
    final_record = list()
    collect_seqs = "false"
    trustable_seqs = "non"
    trustable_seqs_sec = 0
    for line in cluster_with_seed:
        if collect_seqs == "true" and line.startswith(">"):
            if line.startswith(">s"):
                trustable_seqs = "found"
            else:
                tmp_cluster_ids.append(line.split(">")[1])
                if in_directly_stored_lum[int(line.split(">")[1])][5] > 70:
                    trustable_seqs_sec += 1
        if line.startswith("CLUSTER"):
            if len(tmp_cluster_ids) > 0 and (trustable_seqs == "found" or trustable_seqs_sec >= 1):
                # recover all compressed sequences and store data to final_record
                for i in tmp_cluster_ids:
                    seq_id_list = londen_cl_ids_entries[i]
                    for ele in seq_id_list:
                        final_record.append(in_directly_stored_lum[int(ele)])
            tmp_cluster_ids = list()
            trustable_seqs = "non"
            trustable_seqs_sec = 0
            collect_seqs = "true"
    # evaluate final cluster
    if len(tmp_cluster_ids) > 0 and (trustable_seqs == "found" or trustable_seqs_sec >= 1):
        # recover all compressed sequences and store data to final_record
        for i in tmp_cluster_ids:
            seq_id_list = londen_cl_ids_entries[i]
            for ele in seq_id_list:
                final_record.append(in_directly_stored_lum[int(ele)])
    return final_record


def merge_preselection_and_londen(in_only_for_final_output, in_londen_result):
    merged_final_result = list()
    for i in in_only_for_final_output:
        merged_final_result.append(i)
    for i in in_londen_result:
        merged_final_result.append(i)
    return merged_final_result


def create_final_output(in_seek_seq, final_cure_ids, merged_final_result, out_path):
    f_fasta = open(in_seek_seq, 'rU')
    seek_id = ""
    seek_seq = ""
    for record in SeqIO.parse(f_fasta, "fasta"):
        seek_id = record.id
        seek_seq = record.seq
    f_fasta.close()
    if out_path == "stdout":
        print(">" + str(seek_id))
        print(str(seek_seq))
        for tmp_id in final_cure_ids:
            tmp_index_array = final_cure_ids[tmp_id].split(",")
            for index in tmp_index_array:
                # -p.c.VAL: -taxID:
                tmp_arr = str(merged_final_result[int(index)][0]).split("-taxID:")
                out_str = tmp_arr[0] + "-p.c.VAL:" + str(merged_final_result[int(index)][5]) + "%-taxID:" + tmp_arr[1]
                print(out_str)
                print(str(merged_final_result[int(index)][1]))
    else:
        # save to file
        f_out = open(out_path, 'w')
        f_out.write(">" + str(seek_id) + "\n")
        f_out.write(str(seek_seq) + "\n")
        for tmp_id in final_cure_ids:
            tmp_index_array = final_cure_ids[tmp_id].split(",")
            for index in tmp_index_array:
                # -p.c.VAL: -taxID:
                tmp_arr = str(merged_final_result[int(index)][0]).split("-taxID:")
                out_str = tmp_arr[0] + "-p.c.VAL:" + str(merged_final_result[int(index)][5]) + "%-taxID:" + tmp_arr[1]
                f_out.write(out_str + "\n")
                f_out.write(str(merged_final_result[int(index)][1]) + "\n")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--eValue", help="maximal excepted E-Value", type=float, default=1)
    parser.add_argument("-i", "--in_seek_seq", help="search sequence in fasta format", type=str)
    parser.add_argument("-o", "--output_name", help="output file", type=str, default="stdout")
    parser.add_argument("-g", "--gi_pos", help="allowed Gi list", type=str)
    parser.add_argument("-d", "--db", help="path to BLAST db", type=str)
    parser.add_argument("-t", "--num_threads", help="Number of threads", type=int, default=1)
    parser.add_argument("-r", "--blast_runs", help="Maximum number of BLAST runs", type=int, default=2)
    parser.add_argument("-a", "--area_filter", help="Filter distance for Londen ('-1' => auto adjusted)", type=float, default=-1)
    parser.add_argument("-p", "--min_identity", help="Lower limit for identity (0 - 100)%, default -1 (auto mode)", type=float, default=52)
    parser.add_argument("-m", "--max_identity", help="Upper limit for Londen (0 - 100)%, default 80%", type=float, default=80)
    parser.add_argument("-n", "--max_target_seqs", help="Maximum number of aligned sequences to keep, default 2000", type=int, default=2000)
    parser.add_argument("-s", "--secure_mode", help="The secure mode ensures that the number of new queries are "
                                                    "limited. (default \"on\", \"off\")", type=str, default="on")
    parser.add_argument("-l", "--londen_mode", help="Only adjustable, if --area_filter is manually set. "
                                                    "Turn off or select single or double cut method. turn off=0; Single=1; Double=2 (default 2)", type=int, default=2)
    parser.add_argument("-x", "--range1", help="Upper Range [start-end-number of hits]", type=str, default="100-80-10")
    parser.add_argument("-y", "--range2", help="Middle Range [start-end-number of hits]", type=str, default="80-70-45")
    parser.add_argument("-z", "--range3", help="Lower Range [start-end-number of hits]", type=str, default="70-65-45")
    args = parser.parse_args()

    # check paths
    if args.in_seek_seq is None:
        sys.stderr.write("ERROR: Please specify search sequence as input \"-i ./my.fasta\"\n")
        exit()
    if args.db is None:
        sys.stderr.write("ERROR: Please specify the path to the BLAST db \"-b ./my.db\"\n")
        exit()
    if args.blast_runs <= 0:
        sys.stderr.write("ERROR: Only values larger or equal than \"2\" are \"-r 2\"\n")
        exit()

    tmp_fasta_file_container = list()
    # analyze input sequence and modify the sequence (U->T, lower- to uppercase)
    value, input_filename, seq_len = analyze_input_sequence(args.in_seek_seq)
    tmp_fasta_file_container.append(input_filename)

    if value == "false":
        sys.stderr.write("ERROR: Please ensure the correct data format of the input file (DNA/RNA in FASTA format)\n")
        exit()

    # check if auto mode for min. identity is set or not
    # if args.min_identity == 0 -> auto mode else manuel mode
    if args.min_identity == -1:
        new_min_identity = compute_new_min_identity(seq_len)
    else:
        new_min_identity = args.min_identity

    # start algorithm
    directly_stored_londen = list()
    only_for_final_output = list()
    # initial run with standard parameters -> step 1 - 5
    blast_result = blast_seed(input_filename, args.db, args.eValue, args.gi_pos, args.max_target_seqs, args.num_threads)
    fasta_hits = extent_reads_from_blast_table(input_filename, blast_result, args.db, args.eValue)
    fasta_hits_rated = global_alignment(input_filename, fasta_hits)
    (directly_stored_high_sim, directly_stored_tmp, for_reblast) = analyze_sequences(fasta_hits_rated, new_min_identity,
                                                                                     args.max_identity, args.secure_mode,
                                                                                     args.range1, args.range2, args.range3)
    only_for_final_output += directly_stored_high_sim
    directly_stored_londen += directly_stored_tmp

    for i in range((args.blast_runs - 1)):
        # create new fasta file which contains all new queries
        tmp_fasta_file_name = create_fasta_file(for_reblast)
        tmp_fasta_file_container.append(tmp_fasta_file_name)
        # run with n-top-hits (matched best against seed)
        blast_result = blast_seed(tmp_fasta_file_name, args.db, args.eValue, args.gi_pos, args.max_target_seqs, args.num_threads)
        fasta_hits = extent_reads_from_blast_table(input_filename, blast_result, args.db, args.eValue)
        fasta_hits_rated = global_alignment(input_filename, fasta_hits)
        (directly_stored_high_sim, directly_stored_tmp, for_reblast) = analyze_sequences(fasta_hits_rated, new_min_identity,
                                                                                         args.max_identity, args.secure_mode,
                                                                                         args.range1, args.range2, args.range3)
        only_for_final_output += directly_stored_high_sim
        directly_stored_londen += directly_stored_tmp
    # cure directly_stored_londen results
    directly_stored_cure = cure_results(directly_stored_londen)

    num_seq = 0
    input_file_londen = ""
    londen_cl_ids_entries = dict()
    londen_scaling_mode = "relaxed"
    if args.londen_mode != 0 and len(directly_stored_cure) > 0:
        # create fasta file for londen
        (input_file_londen, num_seq, londen_cl_ids_entries, londen_scaling_mode) = create_input_for_londen(input_filename, directly_stored_londen,
                                                                       directly_stored_cure, only_for_final_output, args.num_threads)

    # call modified Londen only if more than one fasta seq in the input and Londen is activated (val 1 or 2)!
    if args.londen_mode != 0 and num_seq > 0:
        londen_result = call_londen(input_file_londen, directly_stored_londen, args.area_filter, args.londen_mode,
                                    londen_cl_ids_entries, londen_scaling_mode)
        # merge only_for_final_output and sequence_clusters (clustered file)
        merged_final_result = merge_preselection_and_londen(only_for_final_output, londen_result)
        final_cure_ids = cure_results(merged_final_result)
        # create final output with original identifiers
        create_final_output(input_filename, final_cure_ids, merged_final_result, args.output_name)
    else:
        if londen_scaling_mode == "relaxed" or args.londen_mode == 0:
            only_for_final_output += directly_stored_londen
        final_cure_ids = cure_results(only_for_final_output)
        create_final_output(input_filename, final_cure_ids, only_for_final_output, args.output_name)
    # clean up - delete tmp files
    if args.londen_mode != 0 and num_seq > 0:
        os.remove(input_file_londen)
    for tmp_files in tmp_fasta_file_container:
        os.remove(tmp_files)
