import gzip
import argparse
from itertools import product
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob
import shutil
import os

def avg_quality_score(sequence):
    '''
    Takes in a sequence of nucleotides, and returns the average quality score
    '''
    sum_quality_score = 0
    # for each character, get phred score using ord(), and adds it to cumulative sum variable
    for char in sequence:
        phred_score = ord(char)-33
        sum_quality_score += phred_score

    # gets avg by dividing sum by length of sequence
    avg_quality_score = sum_quality_score / len(sequence)
    return(avg_quality_score)

def wc_pairs():
    '''
    Returns a dictionary containing Watson-Crick base pairs
    '''
    wc_dict = {'A':'T', 'C':'G', 'T':'A', 'G':'C','N':'N'}
    return(wc_dict)

def reverse_complement(sequence):
    '''
    Takes in a sequence of nucleotides, and returns the reverse complement
    '''
    wc_dict = wc_pairs()
    final_str = ''

    #reverses string
    rev_sequence = sequence[::-1]

    #uses dicts to get watson-crick pair complement
    for char in rev_sequence:
        final_str += wc_dict[char]

    return(final_str)

def write_file_dict(index_list):
    '''
    Creates and returns a dictionary which contains information on index and associated output file
    '''

    fw_file_list = []
    rv_file_list = []

    # make files for fw and rv
    for index in index_list:
        fw_file_list.append(open('fw_' + index + '.fastq','a'))
        rv_file_list.append(open('rv_' + index + '.fastq','a'))

    # creates tuple with (forward file wrapper, reverse file wrapper) for each fw/rv combination
    file_list_tup = tuple(zip(fw_file_list,rv_file_list))

    # creates dictionary where the key is the barcode, value is tuple with (forward file open wrapper, reverse file open wrapper)
    output_file_dict = dict(zip(index_list,file_list_tup))

    return(output_file_dict)

def make_heatmap_matrix(index_list):

    '''
    Makes initial matrix that will be used to keep track of values plotted in the final heatmap.
    '''

    # creates dictionary with each possible index combination as the key, and 0 as the value.
    pairs_dict = dict.fromkeys(list(product(index_list,index_list)),0)

    # creates list from the dict values
    heatmap_vals = list(pairs_dict.values())

    # creates list of lists, with each list representing a barcode
    # each element within a sublist represents one of the 24 barcodes
    # think of it as a 2d matrix, with the same x and y labels
    heatmap_vals = [heatmap_vals[i:i+24] for i in range(0, len(heatmap_vals), 24)]

    return(pairs_dict, heatmap_vals)

def plot_heatmap(index_list, heatmap_vals):
    '''
    Takes in a list of indexes and the associated counts for each index pair, and plots a heatmap.
    '''
    # gets all records by individual barcode
    sum_records = [sum(sublist) for sublist in heatmap_vals]

    # divides each column by the total WITHIN EACH INDIVIDUAL COLUMN
    # this is different than what is being reported at the end of the script
    # i believe this is more useful to show, because it visualizes index hopping irrespective of how many times the index appears in the entire sample
    # each column and row within the heatmap should add up to 100%.
    # note that the colorbar is on a range from 0-1%, in order to make the "high" index hopped combinations pop-out on the heatmap (otherwise, those colors would be drowned out by the correctly matched pairs, the diagonals)
    heatmap_vals = [[100*j/sum_records[i] for j in sublist] for i, sublist in enumerate(heatmap_vals)]

    x = np.array((heatmap_vals))
    fig, ax = plt.subplots(figsize=(15,15))
    sns.heatmap(x, square=True, ax=ax, annot=True, cmap="Blues", xticklabels=index_list, yticklabels=index_list, vmin=0, vmax=1, cbar_kws={'label': '%'}, annot_kws={"size": 10})
    plt.yticks(rotation=0,fontsize=16);
    plt.xticks(fontsize=12);
    plt.tight_layout()
    plt.title('Heatmap of Index Frequency Within Each Index')
    plt.xlabel('P5 index')
    plt.ylabel('P7 index')
    fig.savefig('heatmap.png')

def gzip_files():
    '''
    This function will get a list of files with a .fastq file extension within the current directory, and gzip them. Then, it deletes the uncompressed fastq file.
    '''
    # get list of files with .fastq extension
    file_list = glob.glob('*.fastq')
    for file in file_list:
        with open(file, 'rb') as f_in, gzip.open(file+'.gz', 'wb') as f_out:
            # copy from .fastq to .fastq.gz
            shutil.copyfileobj(f_in, f_out)
        # remove original, uncompressed fastq
        os.remove(file)

def process_files(filenames,threshold=30):
    '''
    This function will process the FASTQ files.
    It will check for correct indexes, index hopping, and bad indexes
    It will also write to specific files based on the above
    '''

    seq_1 = filenames[0]
    seq_2 = filenames[3]
    index_1 = filenames[1]
    index_2 = filenames[2]

    # set counters to count number of occurrences for each situation
    index_hop_counter = 0
    correct_read_counter = 0
    und_read_counter = 0

    # get valid indexes into list
    index_list = ['GTAGCGTA','CGATCGAT','GATCAAGG','AACAGCGA','TAGCCATG','CGGTAATC','CTCTGGAT','TACCGGAT','CTAGCTCA','CACTTCAC','GCTACTCT','ACGATCAG','TATGGCAC','TGTTCCGT','GTCCTAAG','TCGACAAG','TCTTCGAC','ATCATGCG','ATCGTGGT','TCGAGAGT','TCGGATTC','GATCTTGC','AGAGTCCA','AGGATAGC']

    # get output file dictionary from function
    output_file_dict = write_file_dict(index_list)
    pairs_dict, heatmap_vals = make_heatmap_matrix(index_list)

    # open hopped files to write to
    fw_hop_output = open('fw_hopped.fastq','a')
    rv_hop_output = open('rv_hopped.fastq','a')
    # open undetermined files to write to
    fw_und_output = open('fw_undetermined.fastq','a')
    rv_und_output = open('rv_undetermined.fastq','a')

    # store amount of records read (to report later)
    record_counter = 0

    # open the four sequence files
    with gzip.open(seq_1,'r') as s1, gzip.open(seq_2,'r') as s2, gzip.open(index_1,'r') as i1, gzip.open(index_2,'r') as i2:
        print('Started demultiplexing...')
        while True:
            seq_1_list = []
            seq_2_list = []
            index_1_list = []
            index_2_list = []
            # make 4 lists for each file to store current record (should be length of 4 each time per list)
            # append four lines of record to the lists above using readline 4 times (or for/range)
            for i in range(0,4):
                seq_1_list.append(s1.readline().decode('UTF-8').strip())
                seq_2_list.append(s2.readline().decode('UTF-8').strip())
                index_1_list.append(i1.readline().decode('UTF-8').strip())
                index_2_list.append(i2.readline().decode('UTF-8').strip())


            # if first line of seq1 is empty string, means at end of file
            # break out of loop
            if seq_1_list[0] == '':
                break

            record_counter += 1

            # add that combined string to header of each sequence (1st line of seq1/seq2 += concatenated index string above)
            seq_1_list[0] += ' ' + index_1_list[1] + '-' + index_2_list[1]
            seq_2_list[0] += ' ' + index_1_list[1] + '-' + index_2_list[1]
            # add that combined string to header of each sequence (1st line of seq1/seq2 += concatenated index string above)
            combined_index = index_1_list[1] + index_2_list[1]

            # if the indexes are not in the index list, or the quality score of any of the biological/index sequences are less than the threshold (default is 30), then store the record in an undetermined file.
            if (index_1_list[1] not in index_list) or (reverse_complement(index_2_list[1]) not in index_list) or (avg_quality_score(seq_1_list[1]) < threshold) or (avg_quality_score(seq_2_list[1]) < threshold) or (avg_quality_score(index_1_list[1]) < threshold) or (avg_quality_score(index_2_list[1]) < threshold):

                und_read_counter += 1
                fw_und_output.write('\n'.join(seq_1_list) + '\n')
                rv_und_output.write('\n'.join(seq_2_list) + '\n')
                continue

            # if index_1 same as reverse_complement(index_2), this means NO index hopping
            # write record to barcode-specific output file
            if index_1_list[1] == reverse_complement(index_2_list[1]):
                correct_read_counter += 1

                fw_output_file = output_file_dict[index_1_list[1]][0]
                rv_output_file = output_file_dict[index_1_list[1]][1]

                # write R1 list to file, where filename is based on output file dictionary (key is index_1)
                # write R4 list to file, where filename is based on output file dictionary (key is index_1), except has reverse designation in filename
                fw_output_file.write('\n'.join(seq_1_list) + '\n')
                rv_output_file.write('\n'.join(seq_2_list) + '\n')
            # otherwise, index hopping has occurred
            else:
                index_hop_counter += 1

                # write R1 list to a index hopping designated file (forward strand)
                # write R4 list to index hopping designated file (reverse strand)
                fw_hop_output.write('\n'.join(seq_1_list) + '\n')
                rv_hop_output.write('\n'.join(seq_2_list) + '\n')

            # find current p5 barcode position within the index list
            coord_1 = index_list.index(index_1_list[1])
            # find current p7 barcode position within the index list
            coord_2 = index_list.index(reverse_complement(index_2_list[1]))

            # the list of lists with the given coordinates by 1
            heatmap_vals[coord_1][coord_2] += 1

    fw_hop_output.close()
    rv_hop_output.close()
    fw_und_output.close()
    rv_und_output.close()

    # close all correctly paired index read files
    for key in output_file_dict:
        # gets both fw and rv file wrappers
        values = output_file_dict[key]
        for file_wrapper in values:
            # close both of them
            file_wrapper.close()

    # calculates percent of reads which were correct, undetermined, and index hopped.
    percent_correct = str(100*round(correct_read_counter/record_counter,3))
    percent_index_hop = str(100*round(index_hop_counter/record_counter,3))
    percent_und = str(100*round(und_read_counter/record_counter,3))

    # saves heatmap figure
    plot_heatmap(index_list, heatmap_vals)

    # gets heatmap vals as a percentage in order to report
    heatmap_vals_freq = [[j/record_counter for j in sublist] for sublist in heatmap_vals]

    # prints statistics about the demultiplexing process to a text file
    f = open('output_statistics.txt','w')
    print('Finished!')
    print('Total number of reads processed: ' + str(record_counter),file=f)
    print(percent_correct + '% of reads have correctly paired indexes and are above the quality threshold. (' + str(correct_read_counter) + ')',file=f)
    print(percent_index_hop + '% of reads are index hopped. (' + str(index_hop_counter) + ')',file=f)
    print(percent_und + '% of reads are undetermined or below the quality threshold. (' + str(und_read_counter) + ')',file=f)
    print('---------------------------------------',file=f)
    print('***BREAKDOWN OF CORRECTLY MATCHED READS***',file=f)
    print('Index | Count | Percent of total',file=f)
    # [i][i] ensures only matched indexes will be printed
    # the coordinates of the matched indexes are [1][1], [2][2], [3][3] ... think along the diagonal
    for i, barcode in enumerate(index_list):
        print(barcode + ' | ' + str(heatmap_vals[i][i]) + ' | ' + str(100*heatmap_vals_freq[i][i]) + '%',file=f)
    print('---------------------------------------',file=f)
    print('***BREAKDOWN OF INDEX HOPPED READS***',file=f)
    for i, p5_barcode in enumerate(index_list):
        for j, p7_barcode in enumerate(index_list):
            # does not print matched indexes for this section
            if i != j:
                print(p5_barcode+'-'+p7_barcode + ' | ' + str(heatmap_vals[i][j]) + ' | ' + str(100*heatmap_vals_freq[i][j]) + '%',file=f)
    print('Statistics written to output_statistics.txt.')
    f.close()

# this is a common line in python scripts
# this evals to TRUE if the file is called directly, FALSE if it is imported
if __name__ == "__main__":
    # add arguments for argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-r1", "--r1", help="Enter R1 (read 1) filename.")
    parser.add_argument("-r2", "--r2", help="Enter R2 (index 1) filename.")
    parser.add_argument("-r3", "--r3", help="Enter R3 (index 2) filename.")
    parser.add_argument("-r4", "--r4", help="Enter R4 (read 2) filename.")
    args = parser.parse_args()

    # place user's arguments into list
    file_list = [args.r1, args.r2, args.r3, args.r4]

    # threshold refers to phred score
    process_files(file_list, threshold=30)
    gzip_files()
