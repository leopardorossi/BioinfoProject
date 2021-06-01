from .counting import Preprocessing, KmerData
from functools import reduce

import csv


def read_kmers_counting(input_file, seed):
    # Define the array that will contain the counting
    counting = []
    # Define an array that contains the sizes of each family
    families_sizes = []
    for i in range(0, 16):
        families_sizes.append(0)

    # Get the indexes of the ones inside the seed
    ones_positions = [i for i, c in enumerate(seed) if c == "1"]

    # Open the kmer counting file
    with open(input_file, "r") as f:
        # Loop over the lines to retrieve the counting
        for line in f:
            # Split the current line in order to retrieve the kmer and its counting
            kmer, freq = line.split("\t")
            # Produce an integer number that identifies the family of the kmer
            family = Preprocessing.map_first_last_bases(kmer[0], kmer[len(kmer) - 1])
            # Increment the size of the family
            families_sizes[family-1] += 1
            # Build an item to store all the information together
            freq = int(freq.rstrip("\n"))
            # Apply space seed
            masked_kmer = apply_space_seed(ones_positions, kmer)

            counting.append(KmerData(masked_kmer, freq, family))

    return counting, families_sizes


def apply_space_seed(ones_positions, kmer):
    masked_kmer = ""
    for index in ones_positions:
        masked_kmer += kmer[index]

    return masked_kmer


def generate_seed_from_pattern(pattern, kmer_length):
    # Split the pattern into the different groups
    groups = pattern.split("-")
    pattern_length = reduce(lambda a, b: int(a)+int(b), groups)

    if pattern_length > kmer_length:
        raise Exception("The pattern length must be less or equal to the kmer length")

    # Generate the seed
    pattern = ""
    for i, group in enumerate(groups):
        if (i+1) % 2 == 0:
            zeros = "0" * int(group)
            pattern += zeros
        else:
            ones = "1" * int(group)
            pattern += ones

    # Check if the seed is palindrome
    reverse = pattern[::-1]
    if pattern != reverse:
        raise Exception("The seed must be palindrome")

    return pattern


def generate_stats_file(output_path, description):
    stats = {}
    with open(f"{output_path}/result_{description}.txt", "r") as f:
        for line in f:
            kmer, freq = line.split("\t")
            freq = int(freq)
            if freq in stats.keys():
                stats[freq] += 1
            else:
                stats[freq] = 1

    with open(f"{output_path}/profile_{description}.csv", "w") as f:
        f_writer = csv.writer(f, delimiter=";")
        f_writer.writerow(["Frequency", "Count"])
        for key in stats.keys():
            f_writer.writerow([key, stats[key]])


def store_output_file(counting, output_path, description):
    with open(f"{output_path}/result_{description}.txt", "w") as f:
        for key in counting.keys():
            f.write(f"{key}\t{counting[key]}\n")




