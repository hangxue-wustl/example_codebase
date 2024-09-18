import subprocess
import os
import os.path
import click
from Bio import AlignIO, SeqIO
import pandas as pd
import pandas.errors
from seqclass import SeqObject

def prcyan(text):
    click.echo(click.style(text, fg='cyan'))


def prgre(text):
    click.echo(click.style(text, fg='green'))
    

def separate_sequences(input_file, output_dir, continue_analysis=False):
    """
    Separates input file into single separate FASTA files and creates objects for each input sequence
    """
    os.makedirs(output_dir, exist_ok=True)
    seq_list = []

    if not continue_analysis:
        print(
            "TE Trimmer is modifying sequence names; any occurrence of '/', '-', ':', '...', '|' and empty spaces before '#' "
            "will be converted to '_'.\n"
            "You can find the original and modified names in the 'Sequence_name_mapping.txt' file in the output directory.\n")
        # Initialize the name mapping file
        name_mapping_file = os.path.join(os.path.dirname(output_dir), "Sequence_name_mapping.txt")

        detected_pound = False
        with open(input_file, 'r') as fasta_file, open(name_mapping_file, 'w') as mapping_file:
            # Write header to the mapping file
            mapping_file.write("original_input_seq_name\tTEtrimmer_modified_seq_name\n")
            id_list = []
            # Required to add suffix 'fasta', as this pattern will be used for file deletion later
            for record in SeqIO.parse(fasta_file, 'fasta'):
                # Check if '#' is in seq.id. If '#' is found, the string before '#' becomes the seq_name and the string
                # after '#' is the seq_TE_type
                if len(record.id.split("#")) > 1:
                    detected_pound = True
                    sanitized_id = record.id.split("#")[0].replace('/', '_').replace(' ', '_').\
                        replace('-', '_').replace(':', '_').replace('...', '_').replace('|', '_')
                    te_type = record.id.split("#")[1]

                    # Normally SeqIO.parse only takes content before " " as record.id. Separate with " " to make
                    # the code more reliable
                    te_type = te_type.split(" ")[0]

                else:
                    sanitized_id = record.id.replace('/', '_').replace(' ', '_').replace('-', '_')\
                        .replace(':', '_').replace('...', '_')
                    te_type = "Unknown"
                    # modify header to add #Unknown 
                    record.id = f"{record.id}#{te_type}"
                    record.description = record.id

                # double check if sanitized_id is unique. If not, modify sanitized_id
                if sanitized_id not in id_list:
                    id_list.append(sanitized_id)
                else:
                    # print(f"Duplicated seq_name {sanitized_id} during separate_sequences.")
                    id_list.append(sanitized_id)
                    count = id_list.count(sanitized_id)
                    sanitized_id = f"{sanitized_id}_n{count}"

                # Write original and modified names to the mapping file
                mapping_file.write(f"{record.id}\t{sanitized_id}\n")

                # Define output file name
                output_filename = os.path.join(output_dir, f"{sanitized_id}.fasta")
                seq_obj = SeqObject(str(sanitized_id), str(output_filename), len(record.seq), te_type)

                # Store all input file information (object) to seq_list
                seq_list.append(seq_obj)

                # Convert sequence name to sanitized_id
                record.id = sanitized_id

                # Write single FASTA file using sanitized name
                # the record.id is now same as sanitized_id, which only contains string before #
                # the record.description = f"{record.id}#{te_type}"
                with open(output_filename, 'w') as output_file:
                    SeqIO.write(record, output_file, 'fasta')

            if detected_pound:
                print("TEtrimmer detected instances of '#' in your input FASTA sequence headers. The string before "
                      "'#' is denoted as the seq_name, and the string after '#' is denoted as the TE type.\n")
        print("Finish to generate single sequence files.\n")

    elif continue_analysis:
        # If continue_analysis is 'True', generate seq_list based on single FASTA files
        for filename in os.listdir(output_dir):
            file = os.path.join(output_dir, filename)
            with open(file, 'r') as fasta_file:
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    # Get sanitized_id from single FASTA file name
                    sanitized_id = os.path.splitext(filename)[0]

                    # single FASTA file name is the same as sanitized_id and record.id
                    te_type = record.description.split("#")[-1]
                    te_type = te_type.split(" ")[0]
                    seq_obj = SeqObject(str(sanitized_id), str(file), len(record.seq), te_type)
                    seq_list.append(seq_obj)
        print("\nFinished to read single sequence files generated by previous analysis.\n")
    
    single_fasta_n = len(seq_list)
    
    return seq_list, single_fasta_n


def repeatmasker_classification(final_unknown_con_file, final_classified_con_file, classification_dir, num_threads, progress_file,
                                final_con_file, proof_curation_dir, perfect_proof, good_proof, intermediate_proof,
                                need_check_proof, low_copy_dir, hmm, hmm_dir):
    if os.path.exists(final_unknown_con_file) and os.path.exists(final_classified_con_file):
        temp_repeatmasker_dir = os.path.join(classification_dir, "temp_repeatmasker_classification")
        reclassified_recording_path = os.path.join(temp_repeatmasker_dir, "Reclassified_recoring.txt")
        os.makedirs(temp_repeatmasker_dir, exist_ok=True)
        classification_out = repeatmasker(final_unknown_con_file, final_classified_con_file,
                                          temp_repeatmasker_dir, thread=num_threads, classify=True)

        if classification_out:
            repeatmasker_out = os.path.join(temp_repeatmasker_dir,
                                            "temp_TEtrimmer_unknown_consensus.fasta.out")
            reclassified_dict = repeatmasker_output_classify(repeatmasker_out, progress_file,
                                                             min_iden=70, min_len=80, min_cov=0.5)
            if reclassified_dict:
                click.echo(
                    f"\n{len(reclassified_dict)} TE elements were re-classified by the "
                    f"final classification module.")

                # Update final consensus file
                rename_cons_file(final_con_file, reclassified_dict)
                rename_files_based_on_dict(proof_curation_dir, reclassified_dict)
                rename_files_based_on_dict(perfect_proof, reclassified_dict)
                rename_files_based_on_dict(good_proof, reclassified_dict)
                rename_files_based_on_dict(intermediate_proof, reclassified_dict)
                rename_files_based_on_dict(need_check_proof, reclassified_dict)
                rename_files_based_on_dict(low_copy_dir, reclassified_dict, seq_name=True)
                if hmm:
                    rename_files_based_on_dict(hmm_dir, reclassified_dict)

                # Write reclassified ID into a file
                # Open the file for writing
                with open(reclassified_recording_path, 'w') as file:
                    # Iterate through the dictionary and write each key-value pair to the file
                    for key, value in reclassified_dict.items():
                        file.write(f'{key}\t{value}\n')

            else:
                click.echo("0 TE elements were re-classified by the final classification module.")

    else:
        prcyan(f"\nThe final classification module failed.")
        prgre("\nThis does not affect the final TE consensus sequences You can choose to ignore this error.\n")


def repeatmasker(genome_file, library_file, output_dir, thread=1, classify=False):
    """
    Run RepeatMasker with the provided parameters.
    """

    # Construct the RepeatMasker command
    if classify:
        command = ["RepeatMasker",
                   genome_file,
                   "-lib", library_file,
                   "-s",  # Slow search; 0-5% more sensitive, 2-3 times slower than default
                   "-dir", output_dir,
                   "-pa", str(thread)
                   ]
    else:
        command = ["RepeatMasker",
                   genome_file,
                   "-lib", library_file,
                   "-pa", str(thread),
                   "-dir", output_dir,
                   "-s",  # Slow search; 0-5% more sensitive, 2-3 times slower than default
                   "-gff",  # Creates an additional Gene Feature Finding format output
                   "-xm",  # Creates an additional output file in cross_match format (for parsing)
                   "-a",  # Writes alignments in .align output file
                   ]
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return True

    except FileNotFoundError:
        prcyan("'RepeatMasker' command not found. Please ensure 'RepeatMasker' is installed correctly.")
        raise Exception

    except subprocess.CalledProcessError as e:
        if classify:
            prcyan(f"\nRepeatMasker failed during final classification step with error code {e.returncode}")
            prcyan(f"\n{e.stdout}")
            prcyan(f"\n{e.stderr}\n")
            prgre("This will not affect the final result. Only the classification of TE may not be correct.")
            raise Exception
        else:
            prcyan(f"\nRepeatMasker failed during final whole-genome TE annotation with error code {e.returncode}")
            prcyan(f"\n{e.stdout}")
            prcyan(f"\n{e.stderr}\n")
            prgre("This does not affect the final TE consensus library. You can perform the final genome-wide TE"
                  " annotation by yourself with RepeatMasker.")
            raise Exception


def repeatmasker_output_classify(repeatmasker_out, progress_file, min_iden=70, min_len=80, min_cov=0.8):
    # Read RepeatMasker output file (.out) into a DataFrame
    # The regex '\s+' matches one or more whitespace characters
    # error_bad_lines=False to skip errors
    try:
        df = pd.read_csv(repeatmasker_out, sep=r'\s+', header=None, skiprows=3, usecols=range(15))
    except pandas.errors.EmptyDataError:
        return False
    except pd.errors.ParserError:
        df = pd.read_csv(repeatmasker_out, sep=r'\s+', header=None, skiprows=3, error_bad_lines=False, usecols=range(15))

    # Rename columns for easier referencing
    df.columns = [
        'score', 'perc_div', 'perc_del', 'perc_ins', 'query_name',
        'query_start', 'query_end', 'query_left', 'strand',
        'repeat_name', 'repeat_class', 'repeat_start',
        'repeat_end', 'repeat_left', 'ID'
    ]

    # Filter rows based on query identity
    df = df[df['perc_div'] <= (100 - min_iden)]

    # Calculate coverage length and add to a new column cov_len
    df['cov_len'] = abs(df['query_start'] - df['query_end'])

    # Select dataframe columns
    df_filter = df[["query_name", "repeat_name", "repeat_class", "cov_len"]]

    # Group by columns and calculate sum of 'cov_len'
    grouped_df = df_filter.groupby(['query_name', 'repeat_name', 'repeat_class'])['cov_len'].sum().reset_index()
    grouped_df_filter = grouped_df[grouped_df['cov_len'] >= min_len]

    # Group by 'repeat_name' and get the index of the row with the maximum 'cov_len'
    idx = grouped_df_filter.groupby('query_name')['cov_len'].idxmax()

    # Use these indices to filter the DataFrame
    max_cov_len_df = grouped_df_filter.loc[idx]

    # Convert max_cov_len_df to a dictionary for faster look-up
    max_cov_dict = max_cov_len_df.set_index('query_name')['cov_len'].to_dict()

    # Read progress file
    progress_df = pd.read_csv(progress_file)

    # Dictionary to store consensus_name and its reclassified_type
    reclassified_dict = {}

    # Iterate over each row in progress_df
    for index, row in progress_df.iterrows():
        consensus_name = row['consensus_name']
        cons_length = row['cons_length']
        cons_type = row['reclassified_type']

        # Check if consensus_name exists in max_cov_dict and compute the ratio
        if consensus_name in max_cov_dict:
            cov_len = max_cov_dict[consensus_name]
            ratio = int(cov_len) / int(cons_length)

            # Check if ratio meets the threshold
            if ratio >= min_cov and "Unknown" in cons_type:
                # Find the corresponding repeat_class for this consensus_name
                repeat_class = max_cov_len_df[max_cov_len_df['query_name'] == consensus_name]['repeat_class'].iloc[0]

                # Modify the reclassified_type column in progress_df
                progress_df.at[index, 'reclassified_type'] = repeat_class

                # Update the dictionary with the new reclassified_type
                reclassified_dict[consensus_name] = repeat_class

    # Save the modified progress_df back to the original file
    progress_df.to_csv(progress_file, index=False, na_rep='NaN')

    return reclassified_dict


    
def classify_single(consensus_fasta):
    """
    Run RepeatClassifier with the provided parameters.
    """

    # Store the current working directory
    original_dir = os.getcwd()

    # Change the working directory to the directory of the consensus_fasta
    os.chdir(os.path.dirname(consensus_fasta))

    # Define RepeatClassifier command, the output file will be stored in the same directory as consensus_fasta
    command = ["RepeatClassifier", "-consensi", consensus_fasta]

    try:
        # Run RepeatClassifier using subprocess
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("'RepeatClassifier' command not found. Please ensure 'RepeatModeler' is installed correctly.")
        return False

    except subprocess.CalledProcessError as e:
        prcyan(f"RepeatClassifier failed for {os.path.basename(consensus_fasta)} with error code {e.returncode}")
        prcyan(f"\n{e.stdout}")
        prcyan(f"\n{e.stderr}\n")
        prgre("This only affects classification but not the consensus sequence. "
              "You can run 'RepeatClassifier -consensi <your_consensus_file>' manually.")
        # Change the working directory back to the original directory
        os.chdir(original_dir)
        return False

    # Change the working directory back to the original directory
    os.chdir(original_dir)

    classified_file = f'{consensus_fasta}.classified'

    # Get the first record of file with classified consensus sequences
    record = next(SeqIO.parse(classified_file, "fasta"))

    # seq_name = record.id.split("#")[0]
    seq_TE_type = record.id.split("#")[-1]

    return seq_TE_type

def long_bed(input_file, output_dir):
    # calculate alignment length in column 6
    df = pd.read_csv(input_file, sep='\t', header=None)

    """
    # Define a threshold for outlier removal
    threshold = 0.5

    # Identify the top 10 sequence lengths
    top_10_lengths = df.nlargest(10, 6)[6]

    # Calculate the mean and standard deviation of the top 10 lengths
    mean_top_10_lengths = top_10_lengths.mean()
    std_top_10_lengths = top_10_lengths.std()

    # Identify values outside the threshold
    df['difference_from_mean'] = df[6] - mean_top_10_lengths

    # ~ is used to negate the condition, thus keeping rows where the condition is not satisfied.
    filtered_df = df[~((abs(df['difference_from_mean']) > threshold * std_top_10_lengths)
                       & (df['difference_from_mean'] < 0))]
    """
    filtered_df = df

    # Conditionally update values in column 1 and column 2 for right extension
    # 1 adjusted to avoid error message "Error: malformed BED entry at line 91. Start was greater than end"
    reset_right_filtered_df = filtered_df.copy()
    reset_right_filtered_df.loc[reset_right_filtered_df[5] == '+', 1] = reset_right_filtered_df.loc[reset_right_filtered_df[5] == '+', 2] - 1
    reset_right_filtered_df.loc[reset_right_filtered_df[5] == '-', 2] = reset_right_filtered_df.loc[reset_right_filtered_df[5] == '-', 1] + 1

    # Conditionally update values in column 1 and column 2 for left extension
    reset_left_filtered_df = filtered_df.copy()
    reset_left_filtered_df.loc[reset_left_filtered_df[5] == '+', 2] = reset_left_filtered_df.loc[reset_left_filtered_df[5] == '+', 1] + 1
    reset_left_filtered_df.loc[reset_left_filtered_df[5] == '-', 1] = reset_left_filtered_df.loc[reset_left_filtered_df[5] == '-', 2] - 1

    # write new BED file
    reset_right_long_bed = os.path.join(output_dir, f"{os.path.basename(input_file)}_reset_right_long.bed")
    reset_left_long_bed = os.path.join(output_dir, f"{os.path.basename(input_file)}_reset_left_long.bed")
    reset_right_filtered_df.to_csv(reset_right_long_bed, sep='\t', index=False, header=None)
    reset_left_filtered_df.to_csv(reset_left_long_bed, sep='\t', index=False, header=None)

    return reset_left_long_bed, reset_right_long_bed


def extend_end(max_extension, ex_step_size, end, input_file, genome_file, output_dir, crop_end_gap_thr,
               crop_end_gap_win, ext_threshold, define_boundary_win, bed_dic):
    # end: left or right extension
    if_ex = True
    ex_total = 0

    # The following code calculates the majority of the alignment length and finds the longest alignment for extension
    # while ignoring the shorter ones.
    # long_bed function won't do length selection
    # It generates the proper bed file for left and right extension.
    reset_left_long_bed, reset_right_long_bed = long_bed(input_file, output_dir)

    # 100 bp are added to avoid a case where the boundary is at the edge. Therefore, the alignment length is actually
    # ex_step_size + 100 bp
    adjust = 100
    ite = 1
    reset_left_long_bed_copy = reset_left_long_bed
    reset_right_long_bed_copy = reset_right_long_bed
    while if_ex and (ex_total < max_extension):
        # BEDtools will make sure that the extension will not excess the maximum length of the chromosome
        # track extend total
        ex_total += ex_step_size
        
        if end == "left":
            left_ex = ex_total
            right_ex = ex_step_size - ex_total + adjust
            bed_fasta, bed_out_flank_file = extract_fasta(reset_left_long_bed, genome_file, output_dir,
                                                          left_ex, right_ex, nameonly=True)
        if end == "right":
            left_ex = ex_step_size - ex_total + adjust
            right_ex = ex_total
            bed_fasta, bed_out_flank_file = extract_fasta(reset_right_long_bed, genome_file, output_dir,
                                                          left_ex, right_ex, nameonly=True)
        
        # align_sequences() will return extended MSA in absolute file
        bed_fasta_mafft_with_gap = align_sequences(bed_fasta, output_dir)

        if not os.path.isfile(bed_fasta_mafft_with_gap):
            click.echo(f"{input_file} encountered a problem during the MAFFT extension step.")
            return False

        # Remove nucleotide whose proportion is smaller than threshold
        bed_fasta_mafft_with_gap_column_clean_object = CleanAndSelectColumn(bed_fasta_mafft_with_gap, threshold=0.08)
        bed_fasta_mafft_with_gap_column_clean = bed_fasta_mafft_with_gap_column_clean_object.clean_column(output_dir)

        # Remove gaps with similarity check
        bed_fasta_mafft = remove_gaps_with_similarity_check(bed_fasta_mafft_with_gap_column_clean, output_dir,
                                                            gap_threshold=0.6, simi_check_gap_thre=0.4,
                                                            similarity_threshold=0.7,
                                                            min_nucleotide=5
                                                           )

        # bed_fasta_mafft will be false if the MSA column number is smaller than 50 after removing the gaps
        if not bed_fasta_mafft:
            break
        bed_fasta_mafft_object = CropEndByGap(bed_fasta_mafft, gap_threshold=crop_end_gap_thr,
                                              window_size=crop_end_gap_win)
        large_crop_ids, remain_n, remain_ids = bed_fasta_mafft_object.find_large_crops()
        cropped_alignment = bed_fasta_mafft_object.crop_alignment()
        bed_fasta_mafft_cop_end_gap = bed_fasta_mafft_object.write_to_file(output_dir, cropped_alignment)

        # ext_threshold refer to option --ext_thr
        # max_X is float number, means the proportion of X
        bed_boundary = DefineBoundary(bed_fasta_mafft_cop_end_gap, threshold=ext_threshold,
                                      check_window=define_boundary_win, max_X=0.3, if_con_generater=False,
                                      extension_stop_num=300)
        # Read bed_out_flank_file
        bed_out_flank_file_df = pd.read_csv(bed_out_flank_file, sep='\t', header=None)

        # Update bed_dic
        for i in bed_dic:
            id_list = large_crop_ids + remain_ids
            if i in id_list:
                matching_rows = bed_out_flank_file_df[bed_out_flank_file_df[3] == i]
                if end == "left":
                    if bed_dic[i][2] == "+":
                        bed_dic[i][0] = matching_rows.iloc[0, 1]
                    elif bed_dic[i][2] == "-":
                        bed_dic[i][1] = matching_rows.iloc[0, 2]
                if end == "right":
                    if bed_dic[i][2] == "+":
                        bed_dic[i][1] = matching_rows.iloc[0, 2]
                    elif bed_dic[i][2] == "-":
                        bed_dic[i][0] = matching_rows.iloc[0, 1]

        # If the remaining sequences that need further extension are fewer than 5, stop the extension
        if not bed_boundary.if_continue or remain_n < 5:
            break
        elif bed_boundary.if_continue:
            # Subset BED file to only keep sequences that need further extension
            if end == "left":
                if_ex = bed_boundary.left_ext
                df = pd.read_csv(reset_left_long_bed, sep='\t', header=None)

                # Filter out rows where the fourth column is in large_crop_ids
                df_filtered = df[~df[3].isin(large_crop_ids)]
                reset_left_long_bed = reset_left_long_bed_copy + f"_{ite}"

                # Write the filtered DataFrame back to reset_left_long_bed
                df_filtered.to_csv(reset_left_long_bed, sep='\t', index=False, header=None)

            if end == "right":
                if_ex = bed_boundary.right_ext
                df = pd.read_csv(reset_right_long_bed, sep='\t', header=None)

                # Filter out rows where the fourth column is in large_crop_ids
                df_filtered = df[~df[3].isin(large_crop_ids)]
                reset_right_long_bed = reset_right_long_bed_copy + f"_{ite}"

                # Write the filtered DataFrame into reset_left_long_bed
                df_filtered.to_csv(reset_right_long_bed, sep='\t', index=False, header=None)
            ite += 1

    return bed_dic, ex_total

