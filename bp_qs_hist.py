import gzip
import matplotlib.pyplot as plt

def convert_phred(letter):
    return(ord(letter) - 33)

def get_avg_quality_scores(filename):
    # counter for number of records (should be total lines divided by 4)
    num_records = 0
    # use gzip package to avoid having to gunzip the file (which is massive)
    with gzip.open(filename,'r') as f:
        for index, line in enumerate(f):
            # print every 1 million lines for debug purposes
            if index % 1000000 == 0:
                print(filename + ': ' + str(index), flush=True)
            # only on the first time it hits a quality score line
            # set empty array equal to that line
            if index == 3:
                line_for_count = line.decode('UTF-8').strip()
                line_length = len(line_for_count)
                arr_sum = [0 for i in range(line_length)]
            if index % 4 == 3:
                # need to decode bc gzip returns binary line
                line = line.decode('UTF-8').strip()
                num_records+=1
                for char_index, char in enumerate(line):
                    # sequentially sums quality scores by base pair position
                    arr_sum[char_index] += convert_phred(char)

    # takes sum and divides by total records to get avg
    arr_avg = [i/num_records for i in arr_sum]
    # gets R1, R2, R3, or R4.
    run_num = filename.strip().split('_')[4]
    return(arr_avg, run_num)

def plot_quality_scores(avg_arr, run_num):
    # larger figsize
    plt.rcParams['figure.figsize'] = [20, 10]
    # larger font
    plt.rcParams.update({'font.size': 22})
    plt.bar(range(len(avg_arr)), avg_arr)
    plt.xlabel('Base pair position')
    plt.ylabel('Average quality score (Phred + 33)')
    if len(avg_arr) != 8:
        plt.title(str(run_num) + ' (sequence file)' + ': Average quality score by base pair position')
    else:
        plt.title(str(run_num) + ' (index file)' + ': Average quality score by base pair position')
    plt.savefig(str(run_num) + 'q_score_dist.png')
    plt.close()


if __name__ == "__main__":
    index_files = ['/gpfs/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz','/gpfs/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz']
    sequence_files = ['/gpfs/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz','/gpfs/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz']
    for filename in index_files+sequence_files:
        arr_avg, run_num = get_avg_quality_scores(filename)
        plot_quality_scores(arr_avg, run_num)
