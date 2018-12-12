from __future__ import print_function

import glob
import os
from collections import defaultdict
from TCRSeqLib import VDJ_func

import six
from Bio import SeqIO


def makeOutputDir(output_dir_path):
    if not os.path.exists(output_dir_path):
        os.makedirs(output_dir_path)


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def get_filename_and_locus(name):
    name = name.split("_")
    cell = name[0]
    locus = "_".join(name[1:3])
    return ([cell, locus])


def sort_locus_names(dictionary_to_sort):
    for key, value in six.iteritems(dictionary_to_sort):
        sorted_value = sorted(value)
        dictionary_to_sort[key] = sorted_value
    return (dictionary_to_sort)


def load_IMGT_seqs(file):
    seqs = {}
    with open(file, 'rU') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seqs[record.id] = str(record.seq)
    return (seqs)


def parse_IgBLAST(fasta_file,receptor, loci, output_dir, cell_name, raw_seq_dir, species,
                  seq_method, max_junc_len=50, invariant_seqs=None):
    IMGT_seqs = dict()
    # expecting_D = dict()

    loci_for_segments = defaultdict(list)

    # for locus in loci:
    #    expecting_D[locus] = False
    for locus in loci:
        seq_files = glob.glob(os.path.join(raw_seq_dir,
                                           "{receptor}_{locus}_*.fa".format(
                                               receptor=receptor,
                                               locus=locus)))
        for f in seq_files:
            # if not f.endswith("_C.fa"):
            segment_name = os.path.splitext(os.path.split(f)[1])[0]
            IMGT_seqs[segment_name] = load_IMGT_seqs(f)
            # if segment_name.split("_")[2] == 'D':
            #    expecting_D[locus] = True
            loci_for_segments[segment_name.split("_")[2]].append(locus)

    # segment_names = ['TRAJ', 'TRAV', 'TRBD', 'TRBJ', 'TRBV']
    # IMGT_seqs = dict()
    # for segment in segment_names:
    #    IMGT_seqs[segment] = load_IMGT_seqs("{}/{}.fa".format(imgt_seq_location, segment))

    locus_names = ["_".join([receptor, x]) for x in loci]
    all_locus_data = defaultdict(dict)
    for locus in locus_names:
        file = "{output_dir}/{cell_name}_{locus}.IgBLASTOut".format(
                output_dir=output_dir, cell_name=cell_name,locus=locus)
        if os.path.isfile(file):
            igblast_result_chunks = split_igblast_file(file)
            for chunk in igblast_result_chunks:
                (query_name, chunk_details) = VDJ_func.process_chunk(chunk)
                all_locus_data[locus][query_name] = chunk_details
        else:
            all_locus_data[locus] = None

    cell = VDJ_func.find_possible_alignments(all_locus_data, locus_names, cell_name,
                                        IMGT_seqs, output_dir, species, seq_method,
                                        invariant_seqs, loci_for_segments, receptor,
                                        loci, max_junc_len,fasta_file)
    return (cell)


def split_igblast_file(filename):
    # code adapted from http://stackoverflow.com/questions/19575702/pythonhow-to-split-file-into-chunks-by-the-occurrence-of-the-header-word
    token = '# IGBLASTN'
    chunks = []
    current_chunk = []

    with open(filename) as fh:
        for line in fh:
            line = line.rstrip()

            if line.startswith(token) and current_chunk and not line.startswith(
                    "Total "):
                # if line starts with token and the current chunk is not empty
                chunks.append(current_chunk[:])  # add not empty chunk to chunks
                current_chunk = []  # make current chunk blank
            # just append a line to the current chunk on each iteration
            if not line.startswith("Total "):
                current_chunk.append(line)

        chunks.append(current_chunk)  # append the last chunk outside the loop
    return (chunks)


def check_binary(name, user_path=None):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    if user_path:
        if is_exe(user_path):
            return user_path
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, name)
            if is_exe(exe_file):
                return exe_file

    raise OSError(
        "Required binary not found: {name}. Please add to PATH or specify location in config file."
        .format(name=name))

