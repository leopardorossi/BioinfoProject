from python_approach.core import *

import sys

if __name__ == '__main__':

    # Parse command line arguments
    params = {}
    for i in range(len(sys.argv)-1):
        arg = sys.argv[i]
        i += 1
        if arg == "-m":
            params[arg] = Preprocessing.SortingMethods(int(sys.argv[i]))
        else:
            params[arg] = sys.argv[i]

    seed = ""
    try:
        seed = generate_seed_from_pattern(params["-p"], int(params["-k"]))
    except Exception as e:
        print(e)
        exit(1)

    # Read input file
    kmer_data, families_sizes = read_kmers_counting(params["-i"])

    # Define the counting strategy according to the given parameter
    cStrategy = None
    if params["-m"] == Preprocessing.SortingMethods.FAMILY:
        kmer_data = Preprocessing.order_kmers(Preprocessing.SortingMethods.FAMILY, kmer_data)
        cStrategy = FamilyCountingStrategy(seed, kmer_data, families_sizes, params["-o"])
    else:
        cStrategy = LexicoGraphicalCountingStrategy(seed, kmer_data, params["-o"])

    # Use the strategy to count space seeds
    cStrategy.execute()
