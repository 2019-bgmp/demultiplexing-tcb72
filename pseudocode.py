

def reverse_complement(sequence):
    '''
    Takes in a sequence of nucleotides, and returns the reverse complement
    '''
    # reverse the string
    # take the complement
    # return the complement
    pass

def write_file_dict():
    '''
    Creates and returns a dictionary which contains information on index and associated output file
    '''
    # write dict that contains the indexes given in assignment as key, output file to write to as value


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
