# BioinfoProject

## Approach to follow
We split our implementation in two main phases:
1. Naive implementation of space seeds counting (Python)
2. Distributed space seeds counting using MapReduce approach (Java)

### Naive implementation
1. Produce k-mer count file using Squeakr tool [V]
2. Order k-mers (different tests can be done) [V]
3. Generate space-seeds
   - Random generation of space-seeds (palindromi e non)
4. For each generated space-seed:
   - Apply the seed to the k-mers [V]
   - Put the result into an hash-table to count it
5. Produce the output file with (space-seed, counting)

**Inputs:**
- k-mer counter file (-i)
- space-seed (-p)
- method (-m)
- output file path (-o)
- k-mer length (-k)

**To-do:**
- What to do when the input sequence is larger? Scalability
- Make program independent from the input type

### MapReduce implementation
1. Produce k-mer count file using Squeakr tool
2. Generate space-seeds
3. Map phase:
   - Emit pairs (first-last chars, sequence)
4. Reduce phase:
   - Put the input pairs into a has-table and count them
5. Aggregate final countings