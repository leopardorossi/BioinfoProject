from python_approach.core import *

if __name__ == '__main__':
    
    kmerData = read_kmers_counting("resources/out.txt")
    kmerData = Preprocessing.order_kmers(Preprocessing.SortingMethods.FAMILY, kmerData)

    cStrategy = LexicoGraphicalCountingStrategy("1001001001001001001001001001", kmerData)
    cStrategy.execute()
