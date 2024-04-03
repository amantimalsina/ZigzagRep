# Computing Zigzag Persistence with Representatives using wires
This project implements the computation of zigzag persistence for filtrations with representatives using wires [0].

## Description
The implementation focuses on efficient computation of zigzag persistence and representatives using C++. It parses a given filtration (a sequence of simplices along with their addition or removal operations) and computes the zigzag persistence intervals and their representatives.

## Getting Started

### Dependencies

- C++ compiler (C++11 or higher)
- [Boost](https://www.boost.org/) library for dynamic bitset operations.

Ensure Boost is installed and properly configured in your system.

### Building

Clone this repository and compile the source code using a C++ compiler supporting C++11 or higher. Here's an example using `g++`:

```bash
git clone https://github.com/yourusername/zzrep.git
cd zzrep
g++ -std=c++11 -I /path/to/boost_1_XX_0 main.cpp -o zzrep
```

Replace `/path/to/boost_1_XX_0` with the actual path to your Boost installation if it's not in the standard include path.

### Running the Program
To run the program, you need to provide an input file describing the filtration. The input file should list each simplex, prefixed by an operation (i for insertion, d for deletion), followed by the vertex indices:
```bash
git clone https://github.com/amantimalsina/zzrep.git
cd zzrep
g++ -std=c++11 -I /path/to/boost_1_XX_0 main.cpp -o zzrep
```

The program outputs the persistence intervals and representatives to files named after the input file but with `_pers_` and `_id_to_simp_` suffixes in the `../outputs/ directory`.


### Authors
Aman Timalsina — [@amantimalsina](https://github.com/amantimalsina)

### License 
This project is licensed under the MIT License - see the LICENSE file for details.
