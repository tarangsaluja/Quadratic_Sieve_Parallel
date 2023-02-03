# Quadratic_Sieve_Parallel

Final Project with Tai Thongthai for Parallel and Distributed Systems course. 

Project compared different parallelizations for particular steps of the Quadratic Sieve attack algorithm. 

Code for implementation and final paper are included in this repository. 

Write "make clean" and then "make" to clear old and build new executables.

Sequential:
First run ./step1 N to run the first part of the algorithm. This should generate text files. 
Without clearing those textfiles, run ./step2 N to run the second part of the algorithm. 
Then the three resulting textfiles have the numbers which been factorized and the c
orresponding matrices (one with the powers of the primes for each factorization and another with those powers taken modulo 2).

Parallel versions:
Same instructions and output, except use the command line arguments for numbers of processors and hostfile to run mpi.
Some hostfiles are included. 
Hostfiles should be modified accordingly after running check.py.

Running and Checking: 
Sample values of N can be found in the run.sh files. 
Correctness can be checked by running the check.py script (Python3 check.py) to see that a row 
in the matrix really does correspond to the power of each of the primes in the factor 
base in the factorization of the corresponding smooth number.
