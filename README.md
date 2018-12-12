# TCRseq
Reconstruction of T cell receptor sequences from single-cell RNA-seq data.

## Contents ##
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Usage](#using-TCRseq)
    - [*plate*](#plate)
    - [*cell*](#cell)


## Introduction
This tool was build to process single-cell RNA-seq data of TCRB gene and return an output with statistics about the most abundant VDJ sequence of each cell. 

The code is based on TraCer (https://github.com/Teichlab/tracer) architecture with adjustments to another version of input.


## Installation
TCRseq is written in Python and so can just be downloaded, and run with `python3.5 TCRseq.py`.
Download the latest version and accompanying files from https://github.com/DiklaGel/VDJ_pipeline. 
TCRseq relies on several additional tools and Python modules that you should install.

### Pre-requisites

#### Software 
1. [IgBLAST](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) - required for analysis of assembled contigs. (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/).
2. [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ) 

##### Installing IgBlast 
Downloading the executable files from `ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/<version_number>` is not sufficient for a working IgBlast installation. You must also download the `internal_data` directory (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data) and put it into the same directory as the igblast executable. This is also described in the igblast README file.

You should also ensure to set the `$IGDATA` environment variable to point to the location of the IgBlast executable. For example run `export IGDATA=/<path_to_igblast>/igblast/1.4.0/bin`.

TCRseq uses a configuration file to point it to the locations of files that it needs and a couple of other options.
An example configuration file is included in the repository - `TCRseq.conf`.

TCRseq looks for the configuration file, in descending order of priority, from the following sources:
1. The `-c` option used at run time for any of the TCRseq modes
2. The default global location of `~/.TCRseqrc`
3. If all else fails, the provided example `TCRseq.conf` in the directory is used

**Important:** If you  specify relative paths in the config file these will be used as relative to the main installation directory. For example, `resources/Mmus/igblast_dbs` will resolve to `/<wherever you installed TCRseq>/TCRseq/resources/Mmus/igblast_dbs`.

### External tool locations 
TCRseq will look in your system's `PATH` for external tools. You can override this behaviour by editing your `TCRseq.conf`.
Edit `TCRseq.conf` (or a copy) so that the paths within the `[tool_locations]` section point to the executables for all of the required tools.

	[tool_locations]
	#paths to tools used by TCRseq for alignment, quantitation, etc
	igblast_path = /path/to/igblastn

#### IgBLAST options 
##### Receptor type 
    igblast_seqtype = TCR

Type of sequence to be analysed. 

## Using TCRseq 
TCRseq has two modes: *plate* and *cell*

### *plate*:
##### Description:
Plate mode processes fastq files from single plate, split the reads by cell basrcodes and run TCRseq on each cell.
##### Input:
Fastq file for read1 (TCRB sequences) and read2 (cell_umi barcode) of a specific batch (plate)
##### Output:
CSV file with VDJ statistics about the plate's wells (only those who had mapped barcodes)
##### Main process:
In order to reduce the noise and running time, we need to drop out unwanted reads and split the process to "child" processes - one for each cell.
The filtering of noise will be done as follows:
1. First: filter by read2 (by kmers):
    1.1. Filter reads by their cell and UMI barcode (15-mers) frequencies: reads with 15-mer frequency less than f (by default = 0.96 ) percentile are excluded.
    1.2. Mapping of unmapped cell barcodes to wells by cell barcode similarity of 2 and umi barcode similarity of 1
    1.3. Filter reads by UMI barcode similarity:
        1.3.1. For each UMI barcode, define its consensus UMI: the most abundant UMI with hamming distance <= 2
        1.3.2. For each consensus UMI: keep only kmers from the most abundant well (kmers with less abundant cell barcodes are excluded)
2. Second: filter by read1 (by the gene sequence)
    2.1. The current  (and most cost effective) practice is choosing the reads with the most abundant (top 5) hyper variable regions (positions 80 to 130 in our reads) , So we might have several different reads from the same molecule of origin.

After a fasta file for each individual cell is created, we run TCRseq separately for each cell. 

##### Usage
    TCRseq.py plate [-h] [--ncores <CORES>] [--config_file <CONFIG_FILE>]
                 [--resume_with_existing_files] [--species {Mmus,Hsap}]
                 [--receptor_name RECEPTOR_NAME] [--loci [LOCI [LOCI ...]]]
                 [--full] [--filter FILTER]
                 <FASTQ1> <FASTQ2> <PLATE_NAME> <OUTPUT_DIR>


##### Positional arguments:
    <FASTQ1>              first fastq file - read1
    <FASTQ2>              second fastq file - read2
    <PLATE_NAME>          name of plate for file labels
    <OUTPUT_DIR>          directory for output as <output_dir>/<plate_name>

##### Optional arguments:
    -h, --help            show this help message and exit
    --ncores <CORES>, -p <CORES> number of processor cores to use (default: 1)
    --config_file <CONFIG_FILE>, -c <CONFIG_FILE> config file to use (default: ~/.TCRseqrc)
    --resume_with_existing_files, -r look for existing intermediate files and use those instead of starting from scratch (default: False)
    --species {Mmus,Hsap}, -s {Mmus,Hsap} Species to use for reconstruction (default: Hsap)
    --receptor_name RECEPTOR_NAME Name of receptor to reconstruct (default: TCR)
    --loci [LOCI [LOCI ...]]
                            Space-separated list of loci to reconstruct for
                            receptor (default: ['A', 'B'])
    --full                Continue the full process - after splitting to cells,
                            create new job for each cell (default: False)
    --filter FILTER       umis with more than filter reads (with respect to quantile) will be saved (default: 0.96)
 
### *cell*:
##### Description:
Reconstruct TCR sequences from RNAseq reads for a single cell. Running IgBlast and output statistics about the most abundant VDJ sequence
##### Input:
Fasta file with list of reads: for each molecule (consensus UMI barcode) there is one read (TCRB consensus sequences)
##### Output:
CSV file with one row of the VDJ statistics of this specific cell
##### Main process:
The code for this process is based on Tracer code 
1. Fasta file is used as input to IgBlast and the resulting output text file is processed with a custom parsing script (from Tracer).
2. Reads were classed as representing TCR sequences if they contained gene segments from the correct locus and if their reported V and J alignments had E-values below 5 × 10−3. 
3. Reads were determined to be productive, and CDR3 sequences were translated according to IgBlast output
4. From the set of productive reads, a csv file with IgBlast statistics is genereated


##### Usage
    TCRseq.py cell [-h] [--ncores <CORES>] [--config_file <CONFIG_FILE>]
                 [--resume_with_existing_files] [--species {Mmus,Hsap}]
                 [--receptor_name RECEPTOR_NAME] [--loci LOCI [LOCI ...]]
                 <FASTA> <CELL_NAME> <OUTPUT_DIR>

##### Positional arguments:
    <FASTA>               fasta file
    <CELL_NAME>           name of cell for file labels
    <OUTPUT_DIR>          directory for output as <output_dir>/<cell_name>

##### Optional arguments:
    -h, --help            show this help message and exit
    --ncores <CORES>, -p <CORES> 
        number of processor cores to use (default: 1)
    --config_file <CONFIG_FILE>, -c <CONFIG_FILE>
                        config file to use (default: ~/.TCRseqrc)
    --resume_with_existing_files, -r
                        look for existing intermediate files and use those instead of starting from scratch (default: False)
    --species {Mmus,Hsap}, -s {Mmus,Hsap}
                        Species to use for reconstruction (default: Hsap)
    --receptor_name RECEPTOR_NAME
                        Name of receptor to reconstruct (default: TCR)
    --loci LOCI [LOCI ...]
                        Space-separated list of loci to reconstruct for receptor (default: ['A', 'B'])


