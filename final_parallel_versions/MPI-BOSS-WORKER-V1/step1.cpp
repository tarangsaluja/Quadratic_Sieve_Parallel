/*
Tai Thongthai and Tarang Saluja

This program generates a factor base, based on our key N, and writes it to
file.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <gmp.h>
#include <stdarg.h>
#include <obstack.h>

#include "step1.h"

using namespace std;

int main(int argc, char *argv[]){

    if (argc != 2) {
      cout << "Wrong amount of arguments. Please enter desired N" << endl;
      exit(1);
    }

    char* pq = argv[1];
    cout << "Step1 now running for N=" << pq << endl;

    //Setting some value of N = pq for prime p and prime q.
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, pq, 10);

    //calculate number of digits
    size_t j = mpz_sizeinbase (N, 10);
    int size = static_cast<int>(j);

    //Size of Factorbase to be allocated depending on the number of digits
    int fbs, l;

    if (size < 10) {
      fbs = 150;
    } else if (size < 20) {
      fbs = 300;
    } else if (size < 30) {
      fbs = 600;
    } else if (size < 40){
      fbs = 1200;
    } else if (size < 50){
      fbs = 2400;
    } else {
      fbs = 3200;
    }

    //approximation for upper bound on the fbs-th largest prime and allocation
    //of array to store factor base
    l = 2*fbs*log(2*fbs);
    prime_element * primes = (prime_element *)calloc(fbs, sizeof(prime_element));

    //gets all primes less than l = 2*fbs*log(2*fbs) and the amount
    fbs = getprimes(l, N, primes, fbs);



    //Write primes and solutions a, b suh that a^2 = b^2 = N (mod p) into a
    // text file; also write in nuber of primes
    ofstream fb;
    fb.open ("factorbase.txt");
    for (int i = 0; i < fbs; i++){
      fb << primes[i].p <<" " << primes[i].a << " " << primes[i].b << endl;
    }
    fb.close();


    fb.open ("fb_size.txt");
    fb << fbs << endl;
    fb.close();
    cout << "Step1 done running for N=" << pq << endl;

    return 0;
}


/*
Function which performs the sieve of erathosenes on the interval provided
by [0, l]

Inputs: (1) int l -> upper bound for the kth prime
        (2) array of prime_elements
        (3) int size of factorbase (unitialized)

returns factorbase size. 
*/
int getprimes(int l, mpz_t N, prime_element * primes, int fbs){

    //initialize index and mpz ttype to store prime in question
    int idx = 0;
    mpz_t pp;
    mpz_init(pp);

    //An array of booleans indicating if the val is prime or not
    //initialize as True
    bool *truth = new bool[l+1];
    for (int i = 0; i < l+1; i++){
        truth[i] = true;
    }
    //make 0 and 1 index false
    truth[0] = false;
    truth[1] = false;
    //check all i less than sqrt(l), when something is prime, make
    //all multiples not prime (starting from squared, as others arleady covered)
    for (int i = 0; i*i < l; i++){
        if (truth[i]){
            for (int j=i*i; j<=l; j=j+i){
                truth[j] = false;
            }
        }
    }
    //Check if N is QR mod p with legendre
    for (int i = 0; i<l+1; i++){
        if (truth[i]){
            mpz_set_ui(pp, i);

            //do it for 2
            if (i == 2){
              primes[idx].p = i;
              primes[idx].a = 1;
              primes[idx].b = 1;
              idx++;
            }
            //if qr then shank tonelli until fbs size isn't exceeded
            if (mpz_legendre(N, pp) != -1 && idx < fbs){
              primes[idx].p = i;
              shanktonellis(N, &primes[idx]);
              idx++;
            }
        }
    }
    return idx;
}


/*Shank-Tonellis algorithm, taken verbatim from the Bytopia MPQS
implementation, we did not commment it is very complicated*/
void shanktonellis(mpz_t N, prime_element *prime){

 int res_int;
 mpz_t b, q, res;
 mpz_init(b);
 mpz_init(q);
 mpz_init(res);
 mpz_set(b, N);
 mpz_set_ui(q, prime->p);


 mpz_sqrtm(res, b, q);

 res_int = mpz_get_ui(res);
 prime->a = res_int;
 prime->b = prime->p - prime->a;

 return;

}

/*
Helper function for Tonellis Shanks. Taken verbatim from
Bytopia MPQS implementation. We are unable to provide
more comments, as we do not fully grasp the nature 
of this function
*/
int mpz_sqrtm(mpz_ptr rop, mpz_t a, mpz_t q)
{
  mpz_t g, temp, t, gInv, qDiv, h, b;
  int i, s, e, y;

  if (mpz_legendre(a,q) == -1) {
    return 0;
  }

  mpz_init(g); mpz_init(temp);

  mpz_sub_ui(temp, q, 1);


  while (mpz_legendre( g, q ) != -1)
    {
      mpz_random(g, 2);
      mpz_mod(g, g, temp);
      mpz_add_ui(g, g, 1);
    }


  mpz_init_set(t, q);
  mpz_sub_ui(t, t, 1);
  s = mpz_scan1(t, 0);
  mpz_tdiv_q_2exp(t, t, s);

  e = 0;

  mpz_init(gInv);
  if (!mpz_invert(gInv, g, q))
    return 0;

  mpz_init(qDiv);
  mpz_init(h);
  for (i = 2; i <= s; i++)
    {
      mpz_powm_ui(temp, gInv, e, q);
      mpz_mul(h, a, temp);
      mpz_sub_ui(temp, q, 1);
      mpz_tdiv_q_2exp(qDiv, temp, i);
      mpz_powm(temp, h, qDiv, q);
      if (mpz_cmp_ui(temp, 1 ) != 0)
{
 y = 1 << (i - 1);
 e += y;
}
    }

  mpz_powm_ui(temp, gInv, e, q);
  mpz_mul(h, a, temp);

  mpz_init(b);
  mpz_add_ui(t, t, 1);
  mpz_tdiv_q_2exp(t, t, 1);
  mpz_powm(h, h, t, q);
  mpz_powm_ui(g, g, (int)e/2, q);
  mpz_mul(b, g, h);
  mpz_mod(b, b, q);
  mpz_set(rop, b);

  mpz_clear(g); mpz_clear(temp); mpz_clear(t); mpz_clear(gInv);
  mpz_clear(qDiv); mpz_clear(h); mpz_clear(b);
  return 1;
}
