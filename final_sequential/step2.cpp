/*

Tai Thongthai and Tarang Saluja

This program takes in our key N and the factorbase file as input, a 
generates an exponential matrix, a bit matrix, and finds
smooth numbers, and writes them to file.

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

#include "step2.h"


int main(int argc, char *argv[]){

    if (argc != 2) {
      cout << "Wrong amount of arguments. Please enter desired N" << endl;
      exit(1);
    }

    char* pq = argv[1];

    //Set value of N
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, pq, 10);

    size_t j = mpz_sizeinbase (N, 10);
    int size = static_cast<int>(j);

    //sieving chunk size
    int pes = 16000;
    int fbs;
    int relation_count = 0;

    //grab size of factor base from file
    string line;
    ifstream myfile ("fb_size.txt");
    if (myfile.is_open()){
      getline(myfile, line);
      char str_array[line.length()];
      strcpy(str_array, line.c_str());
      fbs = atoi(str_array);
    }

    //fill in factor base
    prime_element * FB = load(N, fbs);


    //Write complete columns as rows into a text file
    sieving_step(FB, N, fbs, pes);



    return 0;
}

/*
Description: load in the factor base elements and the corresponding values
            of a and b which correspond to the smallest values a and b such
            that a^2 - T is 0 mod p and b^2 - T is 0 mod p.
Params: (1) mpz_t number N (product of two primes);
        (2) integer factor base size (fbs)
Return: A pointer to prime_element array. The struct prime_element holds ints,
        prime p, and their solutions a and b.
*/
prime_element * load(mpz_t N, int fbs){

    prime_element * FB = (prime_element *)calloc(fbs, sizeof(prime_element));
    char* number;
    int item;


    string line;
    ifstream myfile ("factorbase.txt");
    if (myfile.is_open()){
        int idx = 0;
        //go through each line and get the required data
        while ( getline (myfile,line) ){
          char str_array[line.length()];
          strcpy(str_array, line.c_str());
          //Use strtok to get one nuber at a time
          number = strtok(str_array, " ");
          item = 0;
          while (number != NULL){
            if (item == 0){
                FB[idx].p = atoi(number);
                item++;
            } else if (item == 1) {
                FB[idx].a = atoi(number);
                item++;
            } else if (item == 2) {
                FB[idx].b = atoi(number);
                item = 0;
            }
            number = strtok(NULL, " ");
        }
          //iterate index, so that it is stored in the correct place
          idx++;
        }
        myfile.close();
    }

    else cout << "Unable to open file";

    return FB;
}

/*
Description: Generates our sieving interval given a size and starting bound.
Params: (1) mpz_t number N (product of two primes);
        (2) int size_SSI (interval size);
        (3) mpz_t number T (interval starting bound starting bound).
Return: A pointer to the polynomial_element array. The polynomial_element
        struct holds an mpz_t type (we named "poly"). This is done
        due to weird memory constraints of mpz_t types
*/
polynomial_element * generate_sieving_interval(mpz_t N, int size_SSI, mpz_t T){

    int size = size_SSI;
    polynomial_element * SI = new polynomial_element[size];

    mpz_t Tsq, res;
    mpz_init(Tsq);
    mpz_init(res);


    //Evaluate for chunk size amount of values and add to array.
    for (int i = 0; i < size_SSI; i++){
        mpz_pow_ui(Tsq, T, 2);
        mpz_sub(res, Tsq, N);

        mpz_init(SI[i].poly);
        mpz_set(SI[i].poly, res);

        mpz_add_ui(T, T, 1);
    }
    return SI;
}



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
void sieving_step(prime_element *FB, mpz_t N, int fbs, int pes){

  //intialize values and recompute value for T
  mpz_t T, T_hold, p, a, b, min1, min2;
  int size_FB = fbs;
  int size_SI = pes;
  int size_SSI = 32000;
  int power;

  unsigned long int init1, init2, step;
  string temp;

  int** power_storage;
  polynomial_element * SI;
  polynomial_element * SI_SAVE;

  mpz_init(min1);
  mpz_init(min2);
  mpz_init(T_hold);
  mpz_init_set_ui(T, 1);
  mpz_root(T, N, 2); // T = sqrt(N)
  mpz_add_ui(T, T, 1); //Buffer T by one to ensure non
  mpz_set(T_hold, T);
  temp = mpz_get_str(NULL, 10, T);

  //counter for total number of relations
  int counter = 0;

  //create text files that we will write to
  ofstream smooth_num_store;
  smooth_num_store.open ("Smooth_Num.txt");

  ofstream power_matrix_store;
  power_matrix_store.open("Power_Matrix.txt");

  ofstream bit_matrix_store;
  bit_matrix_store.open("Bit_Matrix.txt");


  while (counter < size_FB + 10){
    temp = mpz_get_str(NULL, 10, T);

    //initialize matrix which stores max power of a prime which divides
    //any given polynomial evaluation
    power_storage = new int*[size_FB+1];
    for (int i = 0; i < size_FB+1; i++){
      power_storage[i] = new int[size_SSI];
      for (int j = 0; j< size_SSI; j++){
        power_storage[i][j] = 0;
      }
    }

    // Create two copies of sieving interval and maintain value of T
    SI = generate_sieving_interval(N, size_SSI, T);
    mpz_set(T, T_hold);
    SI_SAVE = generate_sieving_interval(N, size_SSI, T);
    mpz_set(T, T_hold);

    //iteration through primes for sieving
    for (int i = 0; i < size_FB; i++){
      //convert p, a, b to mpz types
      mpz_init_set_ui(p, FB[i].p);
      mpz_init_set_ui(a, FB[i].a);
      mpz_init_set_ui(b, FB[i].b);

      //find smallest indices such that the polynomial evaluation at
      //that index is divisble by p in chunk
      init1 = prime_find_min(size_SI, a, p, min1, T_hold);
      init2 = prime_find_min(size_SI, b, p, min2, T_hold);
      // string temp;
      temp = mpz_get_str(NULL, 10, T_hold);
      step = mpz_get_ui (p);

      //go ahead and do all of the divisions
      prime_divide(SI, power_storage, &counter, size_SSI, size_FB, init1,
        step, i);
      prime_divide(SI, power_storage, &counter, size_SSI, size_FB, init2,
        step, i);


      //if there are enough relations, it is time to return
      if (counter >= size_FB + 10){
        break;
      }
    }

    //write the relations to the file
    write_to_file(smooth_num_store, power_matrix_store, bit_matrix_store,
                  size_SSI, size_FB, power_storage, SI_SAVE);

    //clean memory and update smallest value of next chunk
    delete[] SI;
    delete[] SI_SAVE;

    mpz_add_ui(T, T, size_SSI);
    mpz_set(T_hold, T);

    for (int i = 0; i < size_FB + 1; i++){
      delete[] power_storage[i];
    }
    delete[] power_storage;

  }

  if (counter < size_FB + 1){
    // cout << "Not enough relations found" << endl;
    exit(0);
  }
  //return power_storage;

  smooth_num_store.close();
  power_matrix_store.close();
  bit_matrix_store.close();


}

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
  int size_SI, int size_FB, int smallest, int prime, int i){
    int power = 0;
    int step = prime;
    int q;

    //iterate through all numbers which will be divisble by prime
    for (int j = smallest; j < size_SI; j = j + step){

      //keep dividing until 1
      q = mpz_cmp_ui(SI[j].poly, 1);
      if (q != 0){
        power = 0;
        q = mpz_divisible_ui_p(SI[j].poly, step);
        while (q != 0){
          mpz_divexact_ui(SI[j].poly, SI[j].poly, step);
          power += 1;
          q = mpz_divisible_ui_p(SI[j].poly, step);
        }
        //store the power
        power_storage[i][j] += power;

        q = mpz_cmp_ui(SI[j].poly, 1);
        if (q == 0){
          *counter += 1; //iterate counter if it has now been reduced to 1
          power_storage[size_FB][j] = 1;
          if (*counter >= size_FB + 10){ //return if enough relations found
            return;
          }
        }
      }
    }
}


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
int prime_find_min(int size_SI, mpz_t a, mpz_t p, mpz_t min, mpz_t T){

  mpz_t r;
  mpz_t idx;
  mpz_init(r);
  mpz_init(idx);
  int q;

  for (unsigned long int j = 0; j < size_SI; j++){
      //for each prime, figure out the smallest polynomial expressed as
      // (a+pk)^2 - N
      mpz_set_ui(idx, j);
      mpz_add(idx, idx, T);
      mpz_sub(idx, idx, a);
      //check if i congruent to a - T mod p
      mpz_mod(r, idx, p);
      q = mpz_cmp_ui(r, 0);
      if (q == 0){
        mpz_set_ui(min, j);
        break;
      }
    }
    return mpz_get_ui (min);
}

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
  int** power_storage, polynomial_element* SI_SAVE){

  string temp;
  int val;

  for (int j = 0; j < size_SSI; j++){
    if (power_storage[size_FB][j] == 1){
      temp = mpz_get_str(NULL, 10, SI_SAVE[j].poly);
      smooth_num_store << temp << endl;
      for (int i = 0; i < size_FB; i++){
        power_matrix_store << power_storage[i][j];
        val = power_storage[i][j];
        val = power_storage[i][j] % 2;
        bit_matrix_store << val;
      }
      power_matrix_store << endl;
      bit_matrix_store << endl;
    }
  }
}
