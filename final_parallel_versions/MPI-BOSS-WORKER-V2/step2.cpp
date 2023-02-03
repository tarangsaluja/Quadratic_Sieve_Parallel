/*

Tai Thongthai and Tarang Saluja

This program takes in our factorbase file from step1, and
generates the exponential matrix, bit matrix, finds the smooth_numbers,
and writes them all to file.

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
#include <mpi.h>

#include </usr/include/eigen3/Eigen/Dense>

#include "step2.h"


int main(int argc, char *argv[]){

    if (argc != 2) {
      cout << "Wrong amount of arguments. Please enter desired N" << endl;
      exit(1);
    }

    char* pq = argv[1];

    int block_size = 0;
    unsigned int rank = 0;
    unsigned long num_proc = 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, (int*) &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, (int*) &rank);
    //Set value of N
    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, pq, 10);

    size_t j = mpz_sizeinbase (N, 10);
    int size = static_cast<int>(j);

    //current size of sieving interval
    block_size = 16000;

    int relation_count = 0;

    //grab size of factor base from file
    int fbs = 0;
    string line;

    ifstream myfile ("fb_size.txt");
    if (myfile.is_open()){
      getline(myfile, line);
      char str_array[line.length()];
      strcpy(str_array, line.c_str());
      fbs = atoi(str_array);
    } else {
      cout << "FB text file problem" << endl;
    }


    //fill in factor base
    prime_element * FB = load(N, fbs);

    //Write complete columns as rows into a text file
    sieving_step(FB, N, fbs, rank, status, block_size, num_proc);

    MPI_Finalize();
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

    string line;
    ifstream myfile ("factorbase.txt");
    if (myfile.is_open()){
        int idx = 0;
        //go through each line and get the required data
        while ( getline (myfile,line) ){
          char str_array[line.length()];
          strcpy(str_array, line.c_str());
          //Use strtok to get one nuber at a time
          char* number = strtok(str_array, " ");
          int item = 0;
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
Description: step where we repeatedly divide our sieving interval using our
              prime factor base until we have the required number of relations
              and writing the relevant data to file. Rank 1 is the boss node,
              andall other ranks are the worker nodes. Only the workers sieve,
              theboss compiles results
Params: (1) Pointer to prime_element array (the factorbase);
        (2) mpz_t N;
        (3) int factorbasesize (fbs);
        (4) int rank of the current node
        (5) MPI_Status object, status
        (6) int block_size: the chunk size that we generate our temporary s
        ieving  sub intervals on
        (7) int num_proc: the number of processes allocated
Return: Nothing
Output [All Text Files]: (1) Power_Matrix.txt (matrix representation of
                              prime powers of found relations);
                          (2) Bit_Matrix.txt (Power_Matrix mod 2);
                          (3) Smooth_Numbers.txt (Numbers we have found
                          to be smooth)
*/
void sieving_step(prime_element *FB, mpz_t N, int fbs, int rank,
  MPI_Status status, int block_size, int num_proc){

  //intialize values and recompute value for T
  mpz_t T, T_hold, poly;
  int size_FB = fbs;
  int continue_sieving = 1;

  //Find smallest value of T such that T^2 - N >= 0
  mpz_t Tsq, res;
  mpz_init(Tsq);
  mpz_init(res);
  mpz_init(T_hold);
  mpz_init(poly);

  mpz_init_set_ui(T, 1);
  mpz_root(T, N, 2); // T = sqrt(N)
  mpz_add_ui(T, T, 1); //Buffer T by one to ensure non negativity
  mpz_set(T_hold, T);

  int received_processes = 0;
  int tot_relations = 0;
  int max_relations = (size_FB + 10)/(num_proc - 1) + 1;

  //boss process work
  if (rank == 0){
    ofstream smooth_num_file;
    smooth_num_file.open ("Smooth_Num.txt");
    ofstream expo_matrix_file;
    expo_matrix_file.open ("Expo_Matrix.txt");
    ofstream bit_matrix_file;
    bit_matrix_file.open("Bit_Matrix.txt");

    //keep performing master task until received from all processes
    while (received_processes < num_proc - 1){

      master_unpack_save(size_FB, expo_matrix_file, bit_matrix_file,
        smooth_num_file);

      received_processes += 1;
    }
    smooth_num_file.close();
    expo_matrix_file.close();
    bit_matrix_file.close();

 } else{ //worker block

   string* all_smooth = new string[max_relations];
   int** all_relations;
   int** relations;
   string* smooth_nums;
   int** power_storage;

   //create relatins array in advance, knowing how many there are
   all_relations = alloc_2d_int(max_relations, size_FB);

   int block_offset = (rank - 1)*block_size;
   mpz_add_ui(T, T, block_offset);
   mpz_set(T_hold, T);

   int counter = 0;

   //keeping sieving through chunks until enough relations are found
   while (tot_relations < max_relations){

     int relations_amt = 0;
     power_storage = new int*[size_FB+1];
     for (int i = 0; i < size_FB+1; i++){
       power_storage[i] = new int[block_size];
       for (int j = 0; j< block_size; j++){
         power_storage[i][j] = 0;
       }
     }

     polynomial_element * SI = generate_sieving_interval(N, block_size, T);
     mpz_set(T, T_hold);
     polynomial_element * SISAVE = generate_sieving_interval(N, block_size, T);
     mpz_set(T, T_hold);

     worker_sieves(power_storage, &counter, block_size, N, T, size_FB, FB,
       &relations_amt, rank, SI, max_relations);

     //if relations retrieved, then write them into the 2D array of relations
     // and 1D array of strings which will be sent later.
     if (relations_amt > 0){
       smooth_nums = new string[relations_amt];
       relations = alloc_2d_int(relations_amt, size_FB);

       reduce_and_transpose(smooth_nums, relations, power_storage, block_size,
         size_FB, SISAVE);

       update_total(all_smooth, all_relations, max_relations, &tot_relations,
         size_FB, relations, smooth_nums, relations_amt);

       delete [] smooth_nums;

       free(relations[0]);
       free(relations);
     }

     //update T and clean memory
      int offset = block_size * (num_proc - 1);
      mpz_add_ui(T, T, offset);
      mpz_set(T_hold, T);

      for (int i = 0; i < size_FB+1; i ++){
        delete [] power_storage[i];
      }
      delete [] power_storage;

      delete [] SI;
      delete [] SISAVE;

   }

   //worker sender function, once enough relations are found
   worker_pack_send(tot_relations, size_FB, all_relations, all_smooth);

   delete [] all_smooth;
   free(all_relations[0]);
   free(all_relations);

  }
}

/*
Description: Receives expo matrix and smooth nums (creates bit matrix
from expo matrix)
              from workers and saves them to file.
               Unpacks smooth nums before saving to file.
Params:
        (1) int factor base size
        (2) ofstream reference to expo_matrix_file
        (3) ofstream reference to bit_matrix_file
        (4) ofstream reference to smooth_num_file
Return: none
Output [All Text Files]: (1) Power_Matrix.txt (matrix representation of
                              prime powers of found relations);
                          (2) Bit_Matrix.txt (Power_Matrix mod 2);
                          (3) Smooth_Numbers.txt (Numbers we have found
                          to be smooth)
*/
void master_unpack_save(int size_FB, ofstream& expo_matrix_file,
  ofstream& bit_matrix_file, ofstream& smooth_num_file){
  int location, new_relations, size, bit_val, packed_str_length;
  int* smooth_nums_storage;
  int **relations_storage;
  char* str;
  MPI_Status status;

  //Receive number of relations
  MPI_Recv(&new_relations, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
    &status);

  location = status.MPI_SOURCE;
  smooth_nums_storage = new int[new_relations];
  relations_storage = alloc_2d_int(new_relations, size_FB);

  size = new_relations * size_FB;
  //Receive relations
  MPI_Recv(&relations_storage[0][0], size, MPI_INT, location, 0, MPI_COMM_WORLD,
    &status);

  for (int i = 0; i < new_relations; i++){
    for (int j = 0; j < size_FB; j++){
      expo_matrix_file << relations_storage[i][j];
      bit_val = relations_storage[i][j] % 2;
      bit_matrix_file << bit_val;
    }
    expo_matrix_file << endl;
    bit_matrix_file << endl;
  }

  //Receive length char array with packed smooth numbers
  MPI_Recv(&packed_str_length, 1, MPI_INT, location, 0, MPI_COMM_WORLD,
     &status);

  str = new char[packed_str_length];

  //Receive the array
  MPI_Recv(&str[0], packed_str_length, MPI_CHAR, location, 0, MPI_COMM_WORLD,
    &status);

  for (int i = 0; i < packed_str_length; i++){
    if (str[i] != '|' && str[i] != '\0') {
      smooth_num_file << str[i];
    } else if (str[i] == '|') {
      smooth_num_file << endl;
    }
  }

  delete [] smooth_nums_storage;
  delete [] str;
  //free(relations_storage[0]);
  free(relations_storage);
}

/*
Description: Helper function to sieving step: the bulk of the work
              of the worker's sieving task is done here.
Params: (1) int double pointer to power_storage 2D array
        (2) int pointer to counter
        (3) int block size of our temporary sieving subinterval
        (4) mpz_t N our key
        (5) mpz_t T from polynomial Q(T) = T^2 - N
        (6) int size of factor base
        (7) prime_element pointer to our factor base array
        (8) int rank of our current node
        (9) polynomial_element pointer to our sieving subinterval array
Return: Nothing
*/
void worker_sieves(int** power_storage, int* counter, int block_size, mpz_t N,
  mpz_t T, int size_FB, prime_element* FB, int* relations_amt, int rank,
  polynomial_element* SI, int max_relations){

  mpz_t p, a, b, r, idx, min1, min2;
  mpz_init(p);
  mpz_init(a);
  mpz_init(b);
  mpz_init(r);
  mpz_init(idx);
  mpz_init(min1);
  mpz_init(min2);
  int step;


  for (int i = 0; i < size_FB; i++){
      //convert p, a, b to mpz types
      mpz_set_ui(p, FB[i].p);
      mpz_set_ui(a, FB[i].a);
      mpz_set_ui(b, FB[i].b);

      unsigned long init1 = prime_find_min(block_size, a, p, min1, T, r, idx,
        rank);
      unsigned long init2 = prime_find_min(block_size, b, p, min2, T, r, idx,
        rank);
      step = mpz_get_ui (p);

      //go ahead and do all of the divisions
      if (init1 < block_size + 1){
          prime_divide(SI, power_storage, block_size, size_FB, init1, step,
            counter, i);
      }
      if (init2 < block_size + 1){
          prime_divide(SI, power_storage, block_size, size_FB, init2, step,
            counter, i);
      }

      if(*counter >= max_relations){
        break;
      }

    }

    // keep counting how many relations are found
    for (int j = 0; j < block_size; j++){
      if (power_storage[size_FB][j] == 1){
        *relations_amt += 1;
      }
    }
    int rels = *relations_amt;

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
polynomial_element * generate_sieving_interval(mpz_t N, int block_size, mpz_t T){

    int size = block_size;
    polynomial_element * SI = new polynomial_element[size];

    //Find smallest value of T such that T^2 - N >= 0
    mpz_t Tsq, res;
    mpz_init(Tsq);
    mpz_init(res);

    //Evaluate for 80,000 values and add to array.
    for (int i = 0; i < size; i++){
        //res = T^2 - N
        mpz_pow_ui(Tsq, T, 2);
        mpz_sub(res, Tsq, N);

        mpz_init(SI[i].poly);
        mpz_set(SI[i].poly, res);

        mpz_add_ui(T, T, 1);
    }
    return SI;
}

/*
Description: Function for worker to pack smooth_nums as strings and then
            send it to the boss node, using MPI_send
Params: (1) int pointer to relations amount
        (2) int size of factor base
        (3) double int pointer to relations
        (4) string pointer (array of strings) to smooth nums
Return: Nothing
*/
void worker_pack_send(int tot_relations, int size_FB, int** all_relations,
  string* all_smooth){

     int* packed_length = new int;
     *packed_length = 0;
     string packed_smooth_nums;
     char* packed;
     int size = tot_relations * size_FB;

     int ret  = MPI_Send(&tot_relations, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
     ret = MPI_Send(&all_relations[0][0],size, MPI_INT, 0, 0, MPI_COMM_WORLD);

     packed_smooth_nums = pack(packed_length, all_smooth, tot_relations);


     packed = new char[*packed_length];
     strcpy(packed, packed_smooth_nums.c_str());
     ret = MPI_Send(packed_length, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
     ret = MPI_Send(&packed[0], *packed_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

     delete packed_length;
     delete[] packed;
}

/*
Description: Helper function to sieving step:
Params: (1) string array all_smooth
        (2) int double pointer to 2D array of all_relations
        (3) int max relations worker is trying to find
        (4) int pointer to total relations we found so far
        (5) int size of factorbase
        (6) int double pointer to 2D array of relations
        (7) string array to smooth_nums
        (8) int relations amount

Return [Explicit]: Nothing
Return [Implicit - Updates]: (1) string array all_smooth
                             (2) int 2D matrix all_relations
                             (3) int pointer total relations

*/
void update_total(string* all_smooth, int** all_relations, int max_relations,
  int* tot_relations, int size_FB, int** relations, string*smooth_nums,
  int relations_amt){

  int leftover;
  int bound;
  leftover = max_relations - *tot_relations;
  bound = min(relations_amt, leftover);

  for (int i = 0; i < bound; i++){
    all_smooth[*tot_relations+ i] = smooth_nums[i];
    for (int j = 0; j < size_FB; j++){
      all_relations[*tot_relations + i][j] = relations[i][j];
    }
  }

  *tot_relations += bound;
}


/*
Description: Helper function used in sieving_step. Here we divide
            elements in the sieving interval by primes when
            applicable
Params: (1) polynomial_element array pointer SI;
        (2) int double pointer to 2D-power_storage array;
        (3) int the block size of our temporary sieving subinterval
        (4) int factorbase size (size_FB)
        (5) int smallest index for which this prime divides
            this sieving element
        (6) int prime (prime of interest)
        (7) int pointer counter (keeps track of number of relations)
        (8) int i (index of the prime of interest)

Return [Explicit]: Nothing
Return [Implicit]: [Note: We write in place using passed pointers]:
                    (1) polynomial_element array pointer SI;
                    (2) int double pointer to 2D-power_storage array;
                    (3) int pointer counter (keeping track of number
                    of relations)
*/
void prime_divide(polynomial_element* SI, int** power_storage,
                  int block_size, int size_FB, int smallest,
                  int prime, int* counter, int i){
    int power = 0;
    int step = prime;

    //iterate through all numbers which will be divisble by prime
    for (int j = smallest; j < block_size; j = j + step){

      //keep dividing until 1
      int q = mpz_cmp_ui(SI[j].poly, 1);
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
        }
      }
    }
}

/*
Description: Helper function: for each prime, figure out the
              smallest polynomial expressed as (a+pk)^2 - N
Params: (1) int block size
        (2) mpz_t a (solution to tonelli shanks (Q(x) is divisble
            by p when evaluated at numbers of form a+kp and b+kp)
        (3) mpz_t prime p
        (4) mpz_t min
        (5) mpz_t T from Q(T)
        (6) mpz_t r
        (7) mpz_t idx
        (8) int rank of current node

Return: unsigned long: smallest polynomial expressed as (a+pk)^2 - N
*/
unsigned long prime_find_min(int block_size, mpz_t a, mpz_t p, mpz_t min,
  mpz_t T, mpz_t r,
   mpz_t idx, int rank){

  unsigned long temp = block_size +1;

  for (unsigned long j = 0; j < block_size; j++){
      //for each prime, figure out the smallest polynomial expressed as
      // (a+pk)^2 - N
      mpz_set_ui(idx, j);
      mpz_add(idx, idx, T);
      mpz_sub(idx, idx, a);

      //check if i congruent to a - T mod p
      mpz_mod(r, idx, p);
      int q = mpz_cmp_ui(r, (unsigned long) 0);
      if (q == 0){
        mpz_set_ui(min, j);
        temp = mpz_get_ui (min);
        break;
      }
    }
    return temp;
}

/*
Description: Helper function: filter out polynomials that we have
            not found to be relations. Also changes the shape
            representation of our relations 2D array
Params: (1) string array of smooth numbers
        (2) int 2D array of relations (represented as prime powers)
        (3) int 2D array of power_storage
        (4) int block size of our temporary sieving sub_interval
        (5) int size of factor base
        (6) polynomial_element array (our sieving sub_interval)
Return [Explicit]: Nothing
Return [Implicit]: (1) smooth_nums array
                   (2) relations 2D matrix
*/
void reduce_and_transpose(string* smooth_nums, int** relations,
                          int** power_storage, int block_size,
                          int size_FB, polynomial_element *SI){
  string temp;

  int s_idx = 0;
  for (int j = 0; j < block_size; j++){
    if (power_storage[size_FB][j] == 1){
      temp = mpz_get_str(NULL, 10, SI[j].poly);
      smooth_nums[s_idx] = temp;
      for (int idx = 0; idx < size_FB; idx++){
        relations[s_idx][idx] = power_storage[idx][j];
      }
      s_idx++;
    }
  }
}


/*
Description: Helper function for worker nodes sieve - packs smooth nums
              into a string to prepare for MPI send
Params: (1) int pointer to string_length
        (2) string array of smooth_nums
        (3) amount of relations found
Return: string - packed smooth nums string
*/
string pack(int* string_length, string* smooth_nums, int relations_amt){

  string packed_smooth_nums = "";
  for (int i = 0; i < relations_amt; i++){
    packed_smooth_nums = packed_smooth_nums + smooth_nums[i] + "|";
    *string_length = *string_length + smooth_nums[i].length() + 1;
  }
  *string_length = *string_length + 1;

  return packed_smooth_nums;

}

/*
Description: Allocates a contiguous 2D matrix given a row and column size
Params: (1) int row size
        (2) int column size
Returns: A int double pointer to a two dimensional contigous array of ints
Taken verbatim from 
https://coderedirect.com/questions/388175/mpi-matrix-multiplication-with-dynamic-allocation-seg-fault

*/
int **alloc_2d_int(int rows, int cols) {
    int *data = (int *)malloc(rows*cols*sizeof(int));
    int **array= (int **)malloc(rows*sizeof(int*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);
    return array;
}
