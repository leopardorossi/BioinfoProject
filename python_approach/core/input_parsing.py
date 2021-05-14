from .counting import Preprocessing, KmerData


def read_kmers_counting(input_file):
    # Define the array that will contain the counting
    counting = []

    # Open the kmer counting file
    with open(input_file, "r") as f:
        # Loop over the lines to retrieve the counting
        for line in f:
            # Split the current line in order to retrieve the kmer and its counting
            kmer, freq = line.split("\t")
            # Produce an integer number that identifies the class of the kmer 
            family = Preprocessing.map_first_last_bases(kmer[0], kmer[len(kmer) - 1])
            # Build an item to store all the information together
            freq = int(freq.rstrip("\n"))
            counting.append(KmerData(kmer, freq, family))
    
    return counting
