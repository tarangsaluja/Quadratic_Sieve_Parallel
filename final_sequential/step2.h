/*

Tai Thongthai and Tarang Saluja

*/

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <gmp.h>
#include <stdarg.h>
#include <obstack.h>

#include "step1.h"

using namespace std;

#pragma once

/*
struct holding polynomials, an mpz_t types 
*/
struct polynomial_element {
    mpz_t poly;
};

/*
Description: load in the factor base elements and the corresponding values
            of a and b which correspond to the smallest values a and b such
            that a^2 - T is 0 mod p and b^2 - T is 0 mod p.
Params: (1) mpz_t number N (product of two primes);
        (2) integer factor base size (fbs)
Return: A pointer to prime_element array. The struct prime_element holds ints,
        prime p, and their solutions a and b.
*/
prime_element * load(mpz_t N, int fbs);


/*
Description: Generates our sieving interval given a size and starting bound.
Params: (1) mpz_t number N (product of two primes);
        (2) int size_SSI (interval size);
        (3) mpz_t number T (interval starting bound starting bound).
Return: A pointer to the polynomial_element array. The polynomial_element
        struct holds an mpz_t type (we named "poly"). This is done
        due to weird memory constraints of mpz_t types
*/
polynomial_element * generate_sieving_interval(mpz_t N, int size_SSI, mpz_t T);



/*
Description: step where we repeatedly divide our sieving interval using our
              prime factor base until we have the required number of relations
              and writing the relevant data to file.
Params: (1) Pointer to prime_element array (the factorbase);
        (2) mpz_t N;
        (3) int factorbasesize (fbs);
        (4) int polynomial element array size(pes)
Return: Nothing
Output [All Text Files]: (1) Power_Matrix.txt (matrix representation of
                              prime powers of found relations);
                          (2) Bit_Matrix.txt (Power_Matrix mod 2);
                          (3) Smooth_Numbers.txt (Numbers we have found
                          to be smooth)
*/
void sieving_step(prime_element *FB, mpz_t N, int fbs, int pes);


/*
Description: Helper function used in sieving_step. Here we divide
            elements in the sieving interval by primes when
            applicable
Params: (1) polynomial_element array pointer SI;
        (2) int double pointer to 2D-power_storage array;
        (3) int pointer counter (keeps track of number of relations)
        (4) int sieving interval size (size_SI)
        (5) int factorbase size (size_FB)
        (6) int smallest index for which this prime divides
            this sieving element
        (7) int prime (prime of interest)
        (8) int i (index of the prime of interest)
Return [Explicit]: Nothing
Return [Implicit]: [Note: We write in place using passed pointers]:
                    (1) polynomial_element array pointer SI;
                    (2) int double pointer to 2D-power_storage array;
                    (3) int pointer counter (keeping track of number of
                    relations)
*/
void prime_divide(polynomial_element* SI, int** power_storage, int* counter, 
                    int size_SI,  int size_FB, int smallest, int prime, int i);


/*
Description: Finds the smallest index in the current sieving interval
              where the value stored is a multiple of our prime
              of interest, p.
Params: (1) int size of interval (size_SI)
        (2) mpz_t solution (a) to tonelli shanks (Q(x) is divisble
            by p when evaluated at numbers of form a+kp and b+kp)
        (3) mpz_t prime of interest p
        (4) mpz_t min: the smallest index in the sieving intereval
            where the element is divisible by p
        (5) mpz_t T, where T is the input to our polynomial
            function Q(T) = T^2 - N
        (6) mpz_t idx:
Return: unsigned int min: the smallest index in the sieving intereval
            where the element is divisible by p
*/                    
int prime_find_min(int size_SI, mpz_t a, mpz_t p, mpz_t min, mpz_t T);


/*
Description: writing to file smooth numbers, power matrix, and bit matrix.
Params: (1) ofstream reference smooth_numstore
        (2) ofstream reference power_matrix_store
        (3) ofstream reference bit_matrix_store
        (4) int size_SSI: sieving sub interval size (currently sieving
            interval)
        (5) int size_FB: size of factor base
        (6) int double pointer to 2d array of power_storage
        (7) pointer to polynomial_element array SI_SAVE
Return: Nothing
Output [All Text Files]: (1) Power_Matrix.txt (matrix representation of
                              prime powers of found relations);
                          (2) Bit_Matrix.txt (Power_Matrix mod 2);
                          (3) Smooth_Numbers.txt (Numbers we have found
                          to be smooth)
*/
void write_to_file(ofstream& smooth_num_store, ofstream& power_matrix_store, 
    ofstream& bit_matrix_store, int size_SSI, int size_FB, 
    int** power_storage, polynomial_element* SI_SAVE);