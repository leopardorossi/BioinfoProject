from .counting import Counter, KmerData

def read_kmers_countings(inputFile):
    # Define the array that will contain the countings
    countings = []

    # Open the kmer counting file
    with open(inputFile, "r") as f:
        # Loop over the lines to retrieve the countings
        for line in f:
            # Split the current line in order to retrieve the kmer and its counting
            kmer, freq = line.split("\t")
            # Produce an integer number that identifies the class of the kmer 
            family = Counter.mapFirstLastBases(kmer[0], kmer[len(kmer)-1])
            # Build an item to store all the information together
            freq = int(freq.rstrip("\n"))
            countings.append(KmerData(kmer, freq, family))
    
    return countings
