# Space seeds post-processing counter
This repo contains the implementation of a method to apply as a post-processing operation of a k-mer counter, to count spaced k-mers.

## How to use
1. Clone the repository.
2. Open a terminal window and move inside the project folder.
3. Execute the following command

```
./__main__.py -i {input_path} -m {method} -o {output_path} -p {space_seed_pattern} -k {k_size}
```
The parameters to provide are the following:
* -i: input_path = the path to the input file (k-mer counter output)
* -m: method = method to use. It is possible to choose between 1 (Lexicographical sorting) and 2 (Family sorting)
* -o: output_path = the path where to store the output files
* -k: k_size = the length of the k-mers
* -p: space_seed_pattern = the pattern of the space seed, where the odd characters are 1 and the even ones are 0. The pattern must be **palindrome**

*Example*
```
./__main__.py -i resources/inputFile -m 1 -o resources -p 3-3-2-12-2-3-3 -k 28
```
In this example the seed that corresponds to the pattern is: 1110001100000000000011000111

## Team
Alice Codogno
Leonardo Rossi
