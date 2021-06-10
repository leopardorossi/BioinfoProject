from python_approach.core import *

import sys
import time
import tracemalloc

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
        print(seed)
    except Exception as e:
        print(e)
        exit(1)

    time_start = time.perf_counter()

    # Read input file
    kmer_data, families_sizes = read_kmers_counting(params["-i"], seed)

    # Define the counting strategy according to the given parameter
    tracemalloc.start()
    cStrategy = None
    if params["-m"] == Preprocessing.SortingMethods.FAMILY:
        Preprocessing.order_kmers(Preprocessing.SortingMethods.FAMILY, kmer_data)
        cStrategy = FamilyCountingStrategy(seed, kmer_data, families_sizes, params["-o"])
    else:
        Preprocessing.order_kmers(Preprocessing.SortingMethods.LEXICOGRAPHICAL, kmer_data)
        cStrategy = LexicoGraphicalCountingStrategy(seed, kmer_data, params["-o"])

    # Use the strategy to count space seeds
    counting = cStrategy.execute()

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    time_elapsed = (time.perf_counter() - time_start)

    print("%5.1f secs" % time_elapsed)
    print(f"Current memory usage is {current / 10 ** 6}MB; Peak was {peak / 10 ** 6}MB\n")

    store_output_file(counting, params["-o"], cStrategy.description())
    generate_stats_file(params["-o"], cStrategy.description())

