from __future__ import print_function

import matplotlib as mpl

mpl.use('pdf')
import argparse
import sys

from TCRSeqLib.tasks import Plate_Task, Cell_Task


def launch():
    parser = argparse.ArgumentParser(
        description='TCRseq.py: reconstruction of TCR sequences from single-cell RNAseq data',
        usage=''' TCRseq.py <mode> [<args>]

              Modes are :

              - plate: process regular plate fastq file - split them by cell barcodes
              - cell: assemble TCR sequences from single-cell RNA-sequencing reads

              use TCRseq.py <mode> -h for specific help
              ''')
    parser.add_argument('mode', metavar="<MODE>", help='TCRseq.py mode (plate, cell)',
                        choices=['plate', 'cell'])
    args = parser.parse_args(sys.argv[1:2])

    task_mapper = {
        'plate': Plate_Task,
        'cell': Cell_Task,
    }

    if args.mode not in task_mapper:
        print('Unrecognised mode')
        parser.print_help()
        exit(1)

    Task = task_mapper[args.mode]
    Task().run()
