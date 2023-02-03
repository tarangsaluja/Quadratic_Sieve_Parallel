/*
Tai Thongthai and Tarang Saluja


*/
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <gmp.h>
#include <stdarg.h>
#include <obstack.h>


using namespace std;

#pragma once


/*
struct holding prime p, and solutions, a and b, to the polynomial X^2 - N, for a 
given p.
*/
struct prime_element {
    int p, a, b;
};

/*
Function which performs the sieve of erathosenes on the interval provided
by [0, l]

Inputs: (1) int l -> upper bound for the kth prime
        (2) array of prime_elements
        (3) int size of factorbase (unitialized)

returns factorbase size. 
*/
int getprimes(int l, mpz_t N, prime_element * primes, int fbs);

/*Shank-Tonellis algorithm, taken verbatim from the Bytopia MPQS
implementation, we did not commment it is very complicated*/
void shanktonellis(mpz_t N, prime_element *prime);

/*
Helper function for Tonellis Shanks. Taken verbatim from
Bytopia MPQS implementation. We are unable to provide
more comments, as we do not fully grasp the nature 
of this function
*/
int mpz_sqrtm(mpz_ptr rop, mpz_t a, mpz_t q);

