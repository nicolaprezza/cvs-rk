Critical Variables Selection via Rabin-Karp hashing (cvs-rk)
==============

### Brief

Given an input text file containing a matrix A of ASCII characters of size M x L (see file example_input.txt), find a length-L bit-vector B maximizing the counts entropy $H_K(A')$ of the rows of the reduced matrix $A' = A[,B]$, where:

- $A[,B]$ is the matrix obtained from A by removing all columns in positions i such that B[i] = false
- $H_K(A')$ = $-\sum_{k>0}\frac{k\cdot m_k}{M}\log(\frac{k\cdot m_k}{M})$, with $m_k$ being the number of $distinct$ rows (strings) of $A'$ that appear $k$ times in $A'$

### Download

> git clone https://github.com/nicolaprezza/cvs-rk

### Compile

The library has been tested under linux using gcc 4.9.2/5.2.1.

firstly, create and enter a build/ directory

> mkdir build; cd build

Use cmake to generate the makefile (default build type is release):

> cmake ..

Alternatively, to build in debug mode:

> cmake -DCMAKE_BUILD_TYPE=Debug ..

Finally, build the executables:

> make

The above command creates the executables in the build directory.

### Run

from the build/ directory, execute

> ./cvs-rk ../example_input.txt 5 2

The above command executes 2 times the local search strategy on the matrix 'example_input.txt', always keeping the number n of selected columns fixed at n=5.

For more details, run

> ./cvs-rk
