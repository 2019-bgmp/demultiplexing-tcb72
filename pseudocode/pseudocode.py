# DEFINING THE PROBLEM: For Illumina sequencing, it is possible to process multiple samples at the same time. This is done by attaching unique sequences called “barcodes” near the end of the sequence prior to entering the flow cell. This is known as “multiplexing.” Then, after sequencing is complete, we can use programming to separate out the barcodes into independent files. This is known as “demultiplexing.” However, it is possible that the indexes get mixed up when they enter the flow cell, perhaps due to incorrect barcode addition along with extra primers sitting on the flow cell. This is a common phenomena known as “index hopping.” It is also possible that the indexes cannot be resolved, leading to “N” base pairs, or the barcodes do not match any of the barcodes we put into the flow cell. As a result, when demultiplexing, we need to separate the correct index pairs from the index hopped pairs (because we ultimately don’t know which index was hopped) and their attached sequenced. Furthermore, we need to separate unresolved indexes and their accompanied sequences to different files.

# OUTPUT FILES: 24 correct index pairs forward strand files, 24 correct index pairs reverse strand files, 1 index hopped forward strand file, 1 index hopped reverse strand file, 1 “undetermined” forward strand file, and 1 “undetermined” reverse strand file, making for 52 total output files. Furthermore, it would be great information to know the percentage of correct index pairs, percentage of index hopped pairs, and percentage of “undetermined” pairs. In addition, to see if there’s a pattern in index hopping.

# UNIT TESTS: Four files were created which contain an example of each possible outcome above. In each file, the first record is an index hopped read, the second is a correctly paired read, and the third is an “undetermined” read. The output will be 6 files: 1 correct index pairs forward strand files, 1 correct index pairs reverse strand file, 1 index hopped forward strand file, 1 index hopped reverse strand file, 1 undetermined forward strand file, and 1 undetermined reverse strand file. They are named as follows: R1_short_test.fastq, R2_short_test.fastq, R3_short_test.fastq, and R4_short_test.fastq. After running the code, there should be a total of 6 output files: three categories (correct, index hop, undetermined) and its a paired end read (all file categories get a end1 and end2 read), so 3*2 = 6 output files.

# PSEUDOCODE:
def reverse_complement(sequence):
    '''
    Takes in a sequence of nucleotides, and returns the reverse complement
    '''
    # reverse the string
    # take the complement
    # return the complement

    #test input: ATCGCCA
    #test output: TGGCGAT
    pass

def write_file_dict():
    '''
    Creates and returns a dictionary which contains information on index and associated output file
    '''
    # write dict that contains the indexes given in assignment as key, output file to write to as value

    # don't need test input, this function is mainly used for organization purposes
    pass

def process_files(filenames):
    '''
    This function will process the FASTQ files.
    It will check for correct indexes, index hopping, and bad indexes
    It will also write to specific files based on the above, and print number of times they appear

    '''
    # seq_1 = filenames[0]
    # index_1 = filenames[1]
    # etc...
    # set counters to count number of occurrences for each situation
    # get output file dictionary from function
    # with open seq1, index_1, index_2, seq_2:
        #while True
            # make 4 lists for each file to store current record (should be length of 4 each time per list)
            # append four lines of record to the lists above using readline 4 times (or for/range)
            # if first line of seq1 is empty string, means at end of file
                # break out of loop

            # if below quality score threshold, or index1/index2 not in dict:
                # this means bad indexes in some way, shape, or form
                # increment specific counter
                # get combined indexes by concatenating strings (2nd line from index1 list, 2nd line from index2 list)
                # add that combined string to header of each sequence (1st line of seq1/seq2 += concatenated index string above)
                # write R1 to a bad index file (forward strand)
                # write R4 to a bad index file (reverse strand)
                # continue (no point going to next if)

            # if index_1 same as reverse_complement(index_2)
                # this means NO index hopping
                # increment specific counter
                # get combined indexes by concatenating strings (2nd line from index1 list, 2nd line from index2 list)
                # add that combined string to header of each sequence (1st line of seq1/seq2 += concatenated index string above)
                # write R1 list to file, where filename is based on output file dictionary (key is index_1)
                # write R4 list to file, where filename is based on output file dictionary (key is index)1, except has reverse designation in filename
            # else
                # this means index hopping
                # increment specific counter
                # get combined indexes by concatenating strings (2nd line from index1 list, 2nd line from index2 list)
                # add that combined string to header of each sequence (1st line of seq1/seq2 += concatenated index string above)
                # write R1 list to a index hopping designated file (forward strand)
                # write R4 list to index hopping designated file (reverse strand)
    #print counters

    pass

if __name__ == "__main__":
    # call process_files with list of fastq filenames
