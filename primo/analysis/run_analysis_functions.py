import os
import gzip
import shutil
import csv
from fpdf import FPDF
import pandas as pd
from tqdm import tqdm # for progress bar
import matplotlib.pyplot as plt
import random 
from Bio.SeqIO import parse
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def gz_extract(nanopore_data_path):
    """
    For a given directory, unzip all .gz files in its fastq_pass folder,
    save unzipped files in folder.
    """
    data_dir = '/fastq_pass'
    directory = nanopore_data_path + data_dir
    extension = ".gz"
    os.chdir(directory)
    # print(os.getcwd())
    for item in tqdm(os.listdir()): # loop through items in dir
        if item.endswith(extension): # check for ".gz" extension
            gz_name = os.path.abspath(item) # get full path of files
            file_name = (os.path.basename(gz_name)).rsplit('.',1)[0] #get file name for file within
            with gzip.open(gz_name,"rb") as f_in, open(file_name,"wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_name) # delete zipped file

def csv_to_dict(csv_file):
    """
    Given:
    A csv file with two columns
    Returns:
    A dictionary where the keys are the first column and the values are the second column
    """
    with open(csv_file, newline='') as cfile:
        reader = csv.reader(cfile)
        next(reader)
        results = dict(reader) # pull in each row as a key-val pair
    return results

def get_aligned_seqs(filename, query_dict, query_name, score_threshold, MAX_READ_LEN):
    """
    Given:
     fastq filename/path with filename
        ex. 'ALD310_pass_45751867_0.fastq'
     query
        ex. 'AAAGTGTCCGAACGTGTCAA'
     score_threshold
        ex. 100

    Returns:
     dataframe with the sequences in the file that align to the query with alignment score > threshold
        ex. | seq_ids | highest_score | seq    |
            |77c...9ab|     100       | ATCC...|
    """
    file = open(filename)
    query = Seq(query_dict[query_name])

    # records = parse(file, "fastq")

    scores = []
    aligned_seqs = []
    seq_ids = []

    for title, sequence, quality in FastqGeneralIterator(file):
        if len(sequence) > MAX_READ_LEN:
            continue
        # perform alignment
        alignments = pairwise2.align.localxs(sequence, query, -2, -1)
        for a in alignments:
            score = float(a[2])
            # record the sequence and sequence id if the score threshold is passed
            if score > score_threshold:
                scores.append(score)
                aligned_seqs.append(str(sequence))
                seq_ids.append(title.split(' ')[0])
                # only record the highest alignment score for that sequence
                if len(scores) > 1:
                    if seq_ids[-1] == seq_ids[-2]:
                        if scores[-1] > scores[-2]: # keep last info, delete second to last info
                            del scores[-2]
                            del seq_ids[-2]
                            del aligned_seqs[-2]
                        if scores[-2] > scores[-1]: # keep second to last info, delete last info
                            del scores[-1]
                            del seq_ids[-1]
                            del aligned_seqs[-1]
                        if scores[-2] == scores[-1]: # delete last info
                            del scores[-1]
                            del seq_ids[-1]
                            del aligned_seqs[-1]
    seq_alignment_dict = {'read_id': seq_ids,
                          'highest_score': scores,
                          'seq': aligned_seqs,
                         'query': [query_name]*len(scores)}
    sequence_alignment_df = pd.DataFrame(seq_alignment_dict)
    file.close()
    return sequence_alignment_df

def list_of_fastqs_to_analyze(prcnt_data):
    # determine which fastq files to analyze
    file_list = [i for i in os.listdir() if i.endswith('fastq')]
    num_files = int(prcnt_data/float(100)*len(file_list)) +  1
    files_to_analyze = random.sample(file_list, num_files)
    return(files_to_analyze)

def threaded_align_reads(target_dict, fastq_file, SCORE_THRESH, MAX_READ_LEN):
    for query_name in target_dict:
        alignment_df = get_aligned_seqs(fastq_file, target_dict, query_name, SCORE_THRESH, MAX_READ_LEN)
    alignment_df.to_csv(f"aligned_{fastq_file}.csv")

def align_reads(target_dict, data_path, SCORE_THRESH, MAX_READ_LEN, prcnt_data):
    """
    Given:
        A dictionary where the keys are target sequence IDs and values are DNA sequences,
        A string that is the path to the Nanopore data where the FASTQ_PASS dir is
    Returns:
        A pandas dataframe where each row is a read that aligned to 1+ of the target sequences.
        Each row contains the read ID, highest alignment score, read sequence, and the target 
        sequence ID that aligned to it.
    """
    # nanopore_data_path = data_path +  '/fastq_pass'
    
    # determine which fastq files to analyze
    file_list = [i for i in os.listdir() if i.endswith('fastq')]
    if prcnt_data != 100:
        num_files = int(prcnt_data/float(100)*len(file_list)) +  1
        files_to_analyze = random.sample(file_list, num_files)
    else:
        files_to_analyze = file_list
    
    
    # loop through all queries
    file_count = 0
    print('Each progress tick signals one query has been shown to all fastq files.')
    for query_name in tqdm(target_dict): #tqdm is the progress bar
        # loop through all fastq files
        for item in files_to_analyze:
            alignment_df = get_aligned_seqs(item, target_dict, query_name, SCORE_THRESH, MAX_READ_LEN)
            # add this fastq file's aligned sequences to the existing df
            if file_count == 0:
                query_df = alignment_df
            if file_count > 0:
                query_df = query_df.append(alignment_df)
            file_count+=1

    # make each rowname unique (important for plotting)
    query_df = query_df.reset_index()
    print(query_df.shape)
    print(query_df)
    return query_df

def remove_duplicate_reads(data):
    """
    Given:
        The output dataframe from align_reads
    Returns:
        The same dataframe, but all rows with a read_ID seen multiple times are now gone
    """
    # Get all reads that appear multiple times (meaning they aligned to multiple seqs)
    all_read_ids = data['read_id'].tolist()
    seen_reads = []
    duplicates = []
    for read_id in tqdm(all_read_ids):
        if read_id in seen_reads:
            duplicates.append(read_id)
        seen_reads.append(read_id)

    # Remove all reads that align to multiple sequences
    data = data[~data.read_id.isin(duplicates)]

    return data

def create_data_dir(date_label):
    """
    Assuming we're working in the fastq_pass dir:

    Checks to see if a directory for generated data has been made. If not, creates a dir.
    If a dir already exists, does nothing.
    Returns:
        A string that is the data directory name.
    """
    # save to figure directory (make dir if not there, overwrite fig in old dir if already there)
    dir_name = date_label + '_css_analysis'
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    return dir_name

def plot_target_freq(data, date_label):
    """
    Given:
        Dataframe from align_reads or remove_duplicate_reads
    Returns:
        A String that is the name of the saved bar plot of each sequence on the x axis, 
        with number of times it was seen on the y
    """
    data['query'].value_counts().plot(kind='bar')
    plt.ylabel('Number of Reads')
    plt.xlabel('Sequence ID')
    plt.title('Frequency of Each Sequence')
    plt.tight_layout()

    # checks to see if a directory to store data exists, makes one if not
    dir_name = create_data_dir(date_label)

    # save figure, return string that is the path to the figure, including figure
    fig_path_name = dir_name + '/' + 'alignment_freq_barplot.png'
    plt.savefig(fig_path_name, dpi = 500)
    plt.close()
    return fig_path_name

def get_read_lens():
    """
    Assuming you're starting in the fastq_pass dir,

    Returns:
        A list that has the length of all reads in all Fastq files in that dir
    """
    lengths = []
    for item in os.listdir():
        if item.endswith('fastq'):
            file = open(item)
            # get sequence data
            records = parse(file, "fastq")
            for record in records:
                seq = Seq(record.seq)
                lengths.append(len(seq))
    return(lengths)

def plot_len_distribution(lengths, date_label):
    """
    Given:
        A string that is the path to the FASTQ_PASS dir
    Returns:
        A String that is the name of the saved bar plot of each read's length on the x axis,
        with the number of reads that had that length on the y axis
    """
    plt.hist(lengths, log=True, bins=100)
    plt.xlabel('Read Length')
    plt.ylabel('Number of Reads')
    plt.title('Read Distribution of all FASTQ_PASS Data')
    plt.tight_layout()

    # checks to see if a directory to store data exists, makes one if not
    dir_name = create_data_dir(date_label)

    # save figure, return string that is the path to the figure, including figure
    fig_path_name = dir_name + '/' + 'all_read_lengths_barplot.png'
    plt.savefig(fig_path_name, dpi = 500)
    plt.close()
    return fig_path_name

def plot_len_distribution_zoomed(lengths, date_label, MAX_READ_LEN):
    """
    Given:
        A string that is the path to the FASTQ_PASS dir
    Returns:
        A String that is the name of the saved bar plot of each read's length on the x axis 
        (up to the read len cutoff as defined as a global variable), with the number of reads
        that had that length on the y axis
    """
    lengths = [i for i in lengths if i <= MAX_READ_LEN]
    plt.hist(lengths, log=True, bins=100)
    plt.xlabel('Read Length')
    plt.ylabel('Number of Reads')
    plt.title(f'Read Distribution of FASTQ_PASS Data Below Length {MAX_READ_LEN}')
    plt.tight_layout()

    # checks to see if a directory to store data exists, makes one if not
    dir_name = create_data_dir(date_label)

    # save figure, return string that is the path to the figure, including figure
    fig_path_name = dir_name + '/' + 'thresholded_read_lengths_barplot.png'
    plt.savefig(fig_path_name, dpi = 500)
    plt.close()
    return fig_path_name

def get_total_dir_reads(lengths):
    """
    Given:
        A string that is the path to the FASTQ_PASS dir
    Returns:
        A descriptive string of the total number of reads in the FASTQ_PASS dir
    """
    description = f'Total Number of Reads in FASTQ_PASS directory: {len(lengths)}'
    return(description)

def initialize_pdf_expt0(date_label):
    """
    Returns a PDF document we can write data to inside the output data dir.
    """
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size = 12)
    text = "This document contains the analysis of the baseline CRISPR-Cas9 Similarity "
    text +=  f"Search Project Experiment performed on {date_label}"
    pdf.multi_cell(200, 10, txt=text, align = 'C')

    return pdf

####################
####################
####################
# Functions below are just for css_run_analysis.py
####################
####################
####################


def initialize_pdf(date_label):
    """
    Returns a PDF document we can write data to inside the output data dir.
    """
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size = 20)
    text = f"This document contains the analysis of the CRISPR-Cas9 Similarity"
    text +=  f"Search Project Experiment performed on {date_label}"
    pdf.multi_cell(200, 10, txt=text, align = 'C')

    return pdf

def calculate_enrichment_scores(df):
    """
    Given:
        Pandas df of aligned data
    Returns:
        A dataframe where each target ('query') is a row (labeled in the first column), 
        and the second column is the proportion that query appeared
    """
    ratios = df.groupby('query')['query'].count()
    num_total_reads = sum(ratios)
    def ratio(x):
        return float(x)/float(num_total_reads)
    ratios = ratios.apply(ratio)
    return ratios

def plot_es_barplot(expt_ratios, initial_ratios, date_label):
    """
    Given:
        Two dataframes, each that lists all the sequence names and their ratios
    Returns:
        A barplot of each sequence's enrichment score
        A string that is the path to and name of the barplot
    """
    enrichment_score_df = expt_ratios / initial_ratios
    enrichment_score_df.plot(kind='bar', color="#6eb3e3")
    plt.xlabel('Sequence Name')
    plt.ylabel('Enrichment Score')
    plt.title('Enrichment Scores')
    plt.tight_layout()

    # checks to see if a directory to store data exists, makes one if not
    dir_name = create_data_dir(date_label)

    # save figure, return string that is the path to the figure, including figure
    fig_path_name = dir_name + '/' + 'enrichment_score_barplot.png'
    plt.savefig(fig_path_name, dpi = 500)
    plt.close()
    return fig_path_name