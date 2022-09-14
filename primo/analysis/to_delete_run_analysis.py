### Script Overview:
# This script analyzes Oxford Nanopore data for the Cas9
# Similarity Search project and outputs a PDF file with
# the results of various potentially useful analyses, as
# well as a directory of all the figures generated for
# easy use elsewhere.


def main():
    # read in command line arguments
    SEQS_TO_ALIGN = sys.argv[1]
    INITIAL_POOL_INFO = sys.argv[2]
    EXPERIMENT_DIR = sys.argv[3]
    date_label = EXPERIMENT_DIR.split('/')[-1].split('_')[0]
    CURRENT_DIR = os.getcwd()

    # Check that these arguments seem correct:
    if '.csv' not in SEQS_TO_ALIGN:
        raise Exception('ARGUMENT ERROR: The first argument must be a .csv file')
    if INITIAL_POOL_INFO.split('/')[-1] != 'baseline_pool_data.csv':
        raise Exception('ARGUMENT ERROR: The second argument must be the .csv output of css_pool_0_analysis.py')
    if '.' in EXPERIMENT_DIR:
        raise Exception('ARGUMENT ERROR: The third argument should be a path that leads to FASTQ_PASS')
    print('Arguments (probably) successfully input')

    # read in baseline pool information
    unmodified_df = pd.read_csv(INITIAL_POOL_INFO)
    
    # Unzip all fastq files from FASTQ_PASS dir
    print('Unzipping Fastq files')
    gz_extract(EXPERIMENT_DIR)

    # Read in all target seqs (seqs to align reads to) as dictionary
    print('Reading in target sequences')
    target_dict = csv_to_dict(SEQS_TO_ALIGN)

    # Read in FASTQ_PASS nanopore data and see which reads align to any of the seqs to align to
    # NOTE- any read aligning to more than seq will eventually be discarded
    print('Aligning reads to target sequences')
    aligned_reads = align_reads(target_dict, EXPERIMENT_DIR)
    aligned_reads_text = f'Number of aligned reads (with duplicate alignments):{len(aligned_reads)}'

    # Remove all reads that aligned to multiple sequences
    print('Removing all reads with multiple alignments')
    aligned_reads = remove_duplicate_reads(aligned_reads)
    cleaned_aligned_reads_text = f'Number of aligned reads (no duplicate alignments):{len(aligned_reads)}'

    print('Data is being analyzed')
    # Plot raw number of reads each sequence got
    barfig_name = plot_target_freq(aligned_reads, CURRENT_DIR, date_label)

    # Get a list of all the read lengths in the directory
    all_read_lengths_in_dir = get_read_lens(EXPERIMENT_DIR)

    # Get a descriptive string of the total number of reads in the directory
    total_reads_str = get_total_dir_reads(all_read_lengths_in_dir)

    # Plot read length distribution for reads in FASTQ_PASS dir
    len_dist_name = plot_len_distribution(all_read_lengths_in_dir, CURRENT_DIR, date_label)

    # Plot read length distribution of reads length 0 to max_read_len
    len_dist_zoomed_name = plot_len_distribution_zoomed(all_read_lengths_in_dir, CURRENT_DIR, date_label)

    # Calculate enrichment score and plot those results
    initial_ratios = calculate_enrichment_scores(unmodified_df)
    expt_ratios = calculate_enrichment_scores(aligned_reads)
    enrichment_score_barfig_name = plot_es_barplot(expt_ratios, initial_ratios, CURRENT_DIR, date_label)

    # Make PDF and write data to it
    print('Making PDF (this is a slow process, give it a minute or two)')
    pdf = initialize_pdf(date_label)
    text = total_reads_str + '\n' + aligned_reads_text + '\n' + cleaned_aligned_reads_text
    pdf.multi_cell(200,30, txt=text, align='L')

    default_image_width = 150
    default_image_height = 110
    pdf.image(barfig_name, w=default_image_width, h=default_image_height)
    pdf.image(enrichment_score_barfig_name, w=default_image_width, h=default_image_height)
    pdf.image(len_dist_name, w=default_image_width, h=default_image_height)
    pdf.image(len_dist_zoomed_name, w=default_image_width, h=default_image_height)

    pdf.output(f'{create_data_dir(CURRENT_DIR, date_label)}/{date_label}_css_analysis_summary.pdf')


if __name__=="__main__":
    main()