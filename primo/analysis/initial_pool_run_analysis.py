### Script Overview:
# This script analyzes the baseline pool's Oxford Nanopore data
# for the Cas9 Similarity Search project and outputs a .csv file
# with alignment information (read id, highest
# alignment score, read sequence, and the query (AKA aligned sequence))
# along with a PDF file the results of various potentially useful
# analyses, as well as a directory of all the figures generated for
# easy use elsewhere.

def main():
    # read in command line arguments
    SEQS_TO_ALIGN = sys.argv[1]
    INITIAL_POOL_DIR_PATH = sys.argv[2]
    date_label = INITIAL_POOL_DIR_PATH.split('/')[-1].split('_')[0]
    CURRENT_DIR = os.getcwd()

    # Check that these arguments seem correct:
    if '.csv' not in SEQS_TO_ALIGN:
        raise Exception('ARGUMENT ERROR: Your first argument must be a .csv file')
    if len(INITIAL_POOL_DIR_PATH.split('.')) > 1:
        raise Exception('ARGUMENT ERROR: The second argument should be a path that leads to FASTQ_PASS')
    print('Arguments (probably) successfully input')

    # Unzip all fastq files from FASTQ_PASS dir
    print('Unzipping Fastq files')
    gz_extract(INITIAL_POOL_DIR_PATH)

    # Read in all target seqs (seqs to align reads to) as dictionary
    print('Reading in target sequences')
    target_dict = csv_to_dict(SEQS_TO_ALIGN)

    # Read in FASTQ_PASS nanopore data and see which reads align to any of the seqs to align to
    # NOTE- any read aligning to more than seq will be discarded
    print('Aligning reads to target sequences')
    aligned_reads = align_reads(target_dict, INITIAL_POOL_DIR_PATH)
    aligned_reads_text = f'Number of aligned reads (with duplicate reads):{len(aligned_reads)}'

    # Remove all reads that aligned to multiple sequences
    print('Removing all reads with multiple alignments')
    aligned_reads = remove_duplicate_reads(aligned_reads)
    cleaned_aligned_reads_text = f'Number of aligned reads (no duplicate reads):{len(aligned_reads)}'

    # Write dataframe to csv
    print('Writing aligned read data to csv')
    csv_name = f'{create_data_dir(CURRENT_DIR, date_label)}/baseline_pool_data.csv'
    aligned_reads.to_csv(csv_name)

    print('Data is being analyzed')
    # Plot number of reads each sequence got
    barfig_name = plot_target_freq(aligned_reads, CURRENT_DIR, date_label)

    # Get a list of all the read lengths in the directory
    all_read_lengths_in_dir = get_read_lens(INITIAL_POOL_DIR_PATH)

    # Get a descriptive string of the total number of reads in the directory
    total_reads_str = get_total_dir_reads(all_read_lengths_in_dir)

    # Plot read length distribution for reads in FASTQ_PASS dir
    len_dist_name = plot_len_distribution(all_read_lengths_in_dir, CURRENT_DIR, date_label)
    # Plot read length distribution of reads length 0 to max_read_len
    len_dist_zoomed_name = plot_len_distribution_zoomed(all_read_lengths_in_dir, CURRENT_DIR, date_label)

    # Make PDF and write data to it
    print('Making PDF (this is a slow process, give it a minute or two)')
    pdf = initialize_pdf_expt0(date_label)
    text = total_reads_str + '\n' + aligned_reads_text + '\n' + cleaned_aligned_reads_text
    pdf.multi_cell(200,30, txt=text, align='L')

    default_image_width = 150
    default_image_height = 110
    pdf.image(barfig_name, w=default_image_width, h=default_image_height)
    pdf.image(len_dist_name, w=default_image_width, h=default_image_height)
    pdf.image(len_dist_zoomed_name, w=default_image_width, h=default_image_height)

    pdf.output(f'{create_data_dir(CURRENT_DIR, date_label)}/{date_label}_css_analysis_summary.pdf')


if __name__=="__main__":
    main()