#include <iostream>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include "windows.h"

typedef long long int numeric_type;

void printOutput(numeric_type* tab, numeric_type tab_len, double* start, double* stop) {
	//printf("\nTime: %fs\n\n", ((double)(*stop - *start) / 1000.0));
	printf("\nTime: %fs\n\n", *stop - *start);
	char bo;
	printf("A total of %lld primes were found.\nPrint output? y/n\n", tab_len);
	std::cin >> bo;
	if (bo == 'y') {
		printf("\n");
		for (numeric_type i = 0; i < tab_len; i++) {
			printf("%lld ", tab[i]);
		}
		printf("\n");
	}
}

int count_primes_in_sieve(char* sieve, numeric_type lower_bound, numeric_type sieve_size) {
	int c = 0;
	if (sieve[0] == 0) {
		c++;
	}
	for (numeric_type i = (lower_bound < 1) ? 1 : lower_bound; i < sieve_size; i++) {
		if (sieve[i] == 0) {
			c++;
		}
	}
	return c;
}

void print_primes_from_sieve(char* sieve, numeric_type lower_bound_index, numeric_type sieve_size) {
	printf("\n");
	if (sieve[0] == 0) {
		printf("2 ");
	}
	for (numeric_type i = (lower_bound_index < 1) ? 1 : lower_bound_index; i < sieve_size; i++) {
		if (sieve[i] == 0) {
			printf("%lld ", i * 2 + 1);
		}
	}
	printf("\n");
}

void sieve_summary(char* sieve, numeric_type lower_bound, numeric_type sieve_size, double* start, double* stop) {
	//printf("\nTime: %fs\n\n", ((double)(*stop - *start) / 1000.0));
	printf("\nTime: %fs\n\n", *stop - *start);
	char bo;
	int c = count_primes_in_sieve(sieve, lower_bound, sieve_size);
	printf("A total of %lu primes were found.\nPrint output? y/n\n", c);
	std::cin >> bo;
	if (bo == 'y') {
		print_primes_from_sieve(sieve, lower_bound, sieve_size);
	}
}

int readInput(numeric_type *a, numeric_type*b, int *mode) {
	printf("Input MODE[0, 1], positive numbers A and B\n");

	std::cin >> *mode;
	if (!std::cin || *mode <= -1 || *mode >= 2) {
		printf("Variable 'MODE' must be 0 or 1\n");
		return 0;
	}
	std::cin >> *a;
	if (!std::cin) {
		printf("'A' must be of type unsigned long int\n");
		return 0;
	}
	std::cin >> *b;
	if (!std::cin) {
		printf("'B' must be of type unsigned long int\n");
		return 0;
	}

	// THERE IS NO INPUT CORRECTION! ( B > A ) ?
	return 1;
}

char *prepare_sieve(numeric_type a, numeric_type b, int mode, 
	numeric_type *sieve_size, numeric_type *sieve_bottom_filter) {
	
	if (mode == 0) { // mode 0
		*sieve_size = (b + 1) / 2;
	}
	else { // mode 1
		*sieve_size = ((unsigned int)sqrt(b) + 1) / 2;
	}

	char* actual_sieve = (char*)calloc((*sieve_size), sizeof(char));
	if (!actual_sieve) {
		printf("Failure calling nTab = calloc(), probably out-of-memory.\n");
		return 0;
	}

	//for no primes b=1;
	if (b < 2) {
		*sieve_size = 0;
		*actual_sieve = 1;
	}

	*sieve_bottom_filter = 0;

	if (mode == 0) { // if mode 0
		*sieve_bottom_filter = a / 2;
		if (a < 3) {
			*sieve_bottom_filter = 0;
		}
		else {
			actual_sieve[0] = 1;
		}
	}
	// if mode 1 - do nothing
	return actual_sieve;
}



void eratostenes(char *actual_sieve, numeric_type *sieve_size) {

	for (numeric_type i = 1; i < *sieve_size; i++) {

		if (actual_sieve[i] == 0) {

			numeric_type num = (i * 2) + 1;

			for (numeric_type j = i + num; j < *sieve_size; j += num) {
				actual_sieve[j] = 1;
			}
		}
	}
}

void eratostenes_parallel_attempt_1(char* actual_sieve, numeric_type* sieve_size) {
	omp_set_num_threads(8);
	for (numeric_type i = 1; i < *sieve_size; i++) {

		if (actual_sieve[i] == 0) {

			numeric_type num = (i * 2) + 1;
			numeric_type ss = *sieve_size;
			numeric_type j;
#pragma omp parallel for schedule(static) firstprivate(j, i, num, ss)
			for (j = i + num; j < ss; j += num) {
				actual_sieve[j]++;
			}
		}
	}
}

void eratostenes_parallel_attempt_2(char* actual_sieve, numeric_type* sieve_size) {

	const int THREAD_COUNT = 8;

	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_COUNT);

	numeric_type bounds[THREAD_COUNT + 1];
	numeric_type chunk_size = ((*sieve_size) / (THREAD_COUNT)) + 1;
	for (int i = 1; i < THREAD_COUNT; i++) {
		bounds[i] = i * chunk_size;
	}
	bounds[0] = 1;
	bounds[THREAD_COUNT] = *sieve_size;

	numeric_type F_CURR_NUM = 0;
	char F_CURR_THREAD = 0;

#pragma omp parallel
	{
		int thread_num = omp_get_thread_num();
		numeric_type array_bottom = bounds[thread_num];
		numeric_type array_top = bounds[thread_num + 1];
		numeric_type local_array_size = array_top - array_bottom;
		char* local_array = (char*)calloc(local_array_size + 1, sizeof(char));
		if (!local_array) {
			perror("ERROR Sieve copy to thread");
		}
		numeric_type num_start = (array_bottom * 2) + 1;
		while (1) {
			if (F_CURR_NUM) {

				numeric_type holder = num_start % F_CURR_NUM;
				holder = (holder == 0) ? 0 : F_CURR_NUM - holder;
				holder = num_start + holder;
				holder = (holder % 2 == 0) ? holder + F_CURR_NUM : holder;
				holder = holder / 2;

				//for (numeric_type j = holder; j < array_top; j += F_CURR_NUM) {
				for (numeric_type j = holder - array_bottom; j < local_array_size; j += F_CURR_NUM) {
					//actual_sieve[j] = 1;
					local_array[j] = 1;
				}
#pragma omp barrier
#pragma omp barrier
			}
			else if (F_CURR_THREAD == thread_num) {
				for (numeric_type i = array_bottom; i < array_top; i++) {
				//for (numeric_type i = 0; i < local_array_size; i++) {

					//if(actual_sieve[i] == 0){
					if (local_array[i] == 0) {
						//numeric_type num = (i * 2) + 1;
						numeric_type num = ((i + array_bottom) * 2) + 1;
						F_CURR_NUM = num;

						for (numeric_type j = i + num; j < array_top; j += num) {
						//for (numeric_type j = i + num; j < local_array_size; j += num) {
							//actual_sieve[j] = 1;
							local_array[j] = 1;
						}
#pragma omp barrier
						F_CURR_NUM = 0;
#pragma omp barrier
					}
				}
				F_CURR_THREAD = thread_num + 1;
			}
			else if (F_CURR_THREAD == THREAD_COUNT) {
				break;
			}
		}
		memcpy(&(actual_sieve[array_bottom]), local_array, local_array_size * sizeof(char));
		if (local_array) {
			free(local_array);
		}
		
	}

	free(bounds);
}



numeric_type* rewrite_primes(char* actual_sieve, numeric_type sieve_size, numeric_type* primes_len) {

	numeric_type* primes = (numeric_type*)calloc(sieve_size, sizeof(numeric_type));
	numeric_type num;
	*primes_len = 0;

	if (!primes) {
		printf("Failure calling primes = calloc(), probably out-of-memory.\n");
		return 0;
	}
	if (actual_sieve[0] == 0) {
		primes[(*primes_len)++] = 2;
	}
	for (numeric_type i = 1; i < sieve_size; i++) {
		if (actual_sieve[i] == 0) {
			num = (i * 2) + 1;
			primes[(*primes_len)++] = num;
		}
	}

	numeric_type min_primes = (*primes_len == 0) ? 1 : *primes_len;	//max(primes_len, 1); - make sure value is > 0
	primes = (numeric_type*)realloc(primes, min_primes * sizeof(numeric_type));
	if (!primes) {
		printf("Failure calling primes = realloc()\n");
		return 0;
	}
	return primes;
}



numeric_type *find_primes_in_range_concept_1(numeric_type a, numeric_type b, numeric_type sieve_size, 
	numeric_type *primes, numeric_type primes_len, 
	numeric_type * outcome_len) {

	numeric_type finalSize = ((b - a) / 2) + 1;
	numeric_type* outcome = (numeric_type*)calloc(finalSize, sizeof(numeric_type));
	*outcome_len = 0;

	if (!outcome) {
		printf("Failure calling outcome = calloc(), probably out-of-memory.\n");
		return 0;
	}

	//copy arimes that are already calculated to ouput if needed;
	// next step evety odd number in between max(a, sieve_size) to b, if checking for modulo against all primes from the sieve;
	// some numbers (copied here) will make modulo itself, resulting in not being added to fianl solution
	for (numeric_type i = 0; i < primes_len; i++) {
		
		if (primes[i] >= a) {
			outcome[(*outcome_len)++] = primes[i];
		}
	}

	//step_by_step calculation to get the beginning point of diviation
	numeric_type max_sieve_value = sieve_size * 2 + 1;	//max value that has been already checked in the sieve
	numeric_type divide_begin = (a > max_sieve_value) ? a : max_sieve_value;	//where to begin; pick the higher from input_a and sieve_max_val
	numeric_type divide_begin_odd = (divide_begin % 2 == 1) ? divide_begin : divide_begin + 1;	//add 1 if needed; we skip even numbers here too

	for (numeric_type x = divide_begin_odd; x <= b; x += 2) { //for every number in set range

		for (numeric_type i = 1; i < primes_len; i++) {
			if (x % primes[i] == 0) { //check for modulo against all primes from sieve
				break;
			}
			else if (i == primes_len - 1) { //if no dividorss - add to outcome
				outcome[(*outcome_len)++] = x;
				break;
			}
		}
	}

	numeric_type min_outcome = (*outcome_len == 0) ? 1 : *outcome_len;	//max(outcome_len, 1); - make sure value is > 0
	outcome = (numeric_type*)realloc(outcome, min_outcome * sizeof(numeric_type));
	if (!outcome) {
		printf("Failure calling outcome = realloc()\n");
		return 0;
	}
	return outcome;
}

numeric_type* find_primes_in_range_concept_1_parallel(numeric_type a, numeric_type b, numeric_type sieve_size,
	numeric_type* primes, numeric_type primes_len,
	numeric_type* outcome_len) {

	numeric_type finalSize = ((b - a) / 2) + 1;
	numeric_type* outcome = (numeric_type*)calloc(finalSize, sizeof(numeric_type));
	*outcome_len = 0;

	if (!outcome) {
		printf("Failure calling outcome = calloc(), probably out-of-memory.\n");
		return 0;
	}

	//copy arimes that are already calculated to ouput if needed;
	// next step evety odd number in between max(a, sieve_size) to b, if checking for modulo against all primes from the sieve;
	// some numbers (copied here) will make modulo itself, resulting in not being added to final solution
	for (numeric_type i = 0; i < primes_len; i++) {

		if (primes[i] >= a) {
			outcome[(*outcome_len)++] = primes[i];
		}
	}

	//step_by_step calculation to get the beginning point of diviation
	numeric_type max_sieve_value = sieve_size * 2 + 1;	//max value that has been already checked in the sieve
	numeric_type divide_begin = (a > max_sieve_value) ? a : max_sieve_value;	//where to begin; pick the higher from input_a and sieve_max_val
	numeric_type divide_begin_odd = (divide_begin % 2 == 1) ? divide_begin : divide_begin + 1;	//add 1 if needed; we skip even numbers here too

	const int THREAD_NUM = 8;
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	int THREAD_FLAG = 0;

#pragma omp parallel
	{
		int LOCAL_IS_PRIME = 0;
		int thread_num = omp_get_thread_num();
		numeric_type local_len = 0;
		numeric_type* local_out = (numeric_type*)malloc(((finalSize/THREAD_NUM)+1) * (sizeof(numeric_type)));
		if (local_out) {
#pragma omp for schedule(static)
			for (numeric_type x = divide_begin_odd; x <= b; x += 2) { //for every *part of* numbers in set range

				for (numeric_type i = 1; i < primes_len; i++) {
					LOCAL_IS_PRIME = 1;
					if (x % primes[i] == 0) { //check for modulo against all primes from sieve
						LOCAL_IS_PRIME = 0;
						break;
					}
				}
				if (LOCAL_IS_PRIME) {
					local_out[local_len++] = x;
				}
				LOCAL_IS_PRIME = 1;
			}
			
			while (THREAD_FLAG < thread_num) {
				Sleep(1);
			}
			memcpy(&outcome[(*outcome_len)], local_out, local_len * (sizeof(numeric_type)));
			free(local_out);
			(*outcome_len) += local_len;
			THREAD_FLAG++;
			
		}
		else {
			printf("malloc error\n");
		}
	}

	numeric_type min_outcome = (*outcome_len == 0) ? 1 : *outcome_len;	//max(outcome_len, 1); - make sure value is > 0
	outcome = (numeric_type*)realloc(outcome, min_outcome * sizeof(numeric_type));
	if (!outcome) {
		printf("Failure calling outcome = realloc()\n");
		return 0;
	}
	return outcome;
}



numeric_type* find_primes_in_range_concept_2(numeric_type a, numeric_type b, numeric_type sieve_size,
	numeric_type* primes, numeric_type primes_len,
	numeric_type* outcome_len) {

	numeric_type finalSize = ((b - a) / 2) + 1;
	numeric_type* outcome = (numeric_type*)calloc(finalSize, sizeof(numeric_type));
	*outcome_len = 0;

	if (!outcome) {
		printf("Failure calling outcome = calloc(), probably out-of-memory.\n");
		return 0;
	}

	//copy arimes that are already calculated to ouput if needed;
	// next step evety odd number in between max(a, sieve_size) to b, if checking for modulo against all primes from the sieve;
	// some numbers (copied here) will make modulo itself, resulting in not being added to final solution
	for (numeric_type i = 0; i < primes_len; i++) {

		if (primes[i] >= a) {
			outcome[(*outcome_len)++] = primes[i];
		}
	}

	//step_by_step calculation to get the beginning point of diviation
	numeric_type max_sieve_value = sieve_size * 2;	//max value that has been already checked in the sieve
	numeric_type divide_begin = (a > max_sieve_value) ? a : max_sieve_value;	//where to begin; pick the higher from input_a and sieve_max_val
	numeric_type divide_begin_odd = (divide_begin % 2 == 1) ? divide_begin : divide_begin + 1;	//add 1 if needed; we skip even numbers here too

	//prepare sieve for odds in range <a,b>
	numeric_type local_tab_size = ((b - divide_begin_odd) / 2) + 1;
	char* final_out = (char*)calloc(local_tab_size, sizeof(char));
	numeric_type local_tab_size_real = (b < divide_begin_odd || b == 1) ? 0 : local_tab_size;

	for (numeric_type i = 1; i < primes_len; i++) { //for every prime

		//find first odd multiply of given prime in range <begin,b>
		numeric_type holder = divide_begin_odd % primes[i];
		holder = (holder == 0) ? 0 : primes[i] - holder;
		holder = divide_begin_odd + holder;
		holder = (holder % 2 == 0) ? holder + primes[i] : holder;

		numeric_type index = (holder - divide_begin_odd) / 2;

		//for (numeric_type x = holder; x <= b; x += primes[i]*2) { // for every number <begin,b> (no evens)
		for (numeric_type x = index; x < local_tab_size_real; x += primes[i]) {
			final_out[x] = -1; //set to non prime;
		}
	}

	for (numeric_type i = 0; i < local_tab_size_real; i++) {
		if (final_out[i] == 0) {
			numeric_type number = (i * 2) + divide_begin_odd;
			outcome[(*outcome_len)++] = number;
		}
	}

	numeric_type min_outcome = (*outcome_len == 0) ? 1 : *outcome_len;	//max(outcome_len, 1); - make sure value is > 0
	outcome = (numeric_type*)realloc(outcome, min_outcome * sizeof(numeric_type));
	if (!outcome) {
		perror("Failure calling outcome = realloc()\n");
		return 0;
	}
	return outcome;
}

numeric_type* find_primes_in_range_concept_2_parallel_domain(numeric_type a, numeric_type b, numeric_type sieve_size,
	numeric_type* primes, numeric_type primes_len,
	numeric_type* outcome_len) {

	numeric_type finalSize = ((b - a) / 2) + 1;
	numeric_type* outcome = (numeric_type*)calloc(finalSize, sizeof(numeric_type));
	*outcome_len = 0;

	if (!outcome) {
		printf("Failure calling outcome = calloc(), probably out-of-memory.\n");
		return 0;
	}

	//copy arimes that are already calculated to ouput if needed;
	// next step evety odd number in between max(a, sieve_size) to b, if checking for modulo against all primes from the sieve;
	// some numbers (copied here) will make modulo itself, resulting in not being added to final solution
	for (numeric_type i = 0; i < primes_len; i++) {

		if (primes[i] >= a) {
			outcome[(*outcome_len)++] = primes[i];
		}
	}

	//step_by_step calculation to get the beginning point of diviation
	numeric_type max_sieve_value = sieve_size * 2;	//max value that has been already checked in the sieve
	numeric_type divide_begin = (a > max_sieve_value) ? a : max_sieve_value;	//where to begin; pick the higher from input_a and sieve_max_val
	numeric_type divide_begin_odd = (divide_begin % 2 == 1) ? divide_begin : divide_begin + 1;	//add 1 if needed; we skip even numbers here too

	char const NUM_THREADS = 8;
	omp_set_dynamic(0);
	omp_set_num_threads(NUM_THREADS);

	numeric_type local_tab_size = ((b - divide_begin_odd) /2)+1;
	char* final_out = (char*)calloc(local_tab_size, sizeof(char));
	numeric_type local_tab_size_real = (b < divide_begin_odd || b == 1) ? 0 : local_tab_size;

	//char* tabs[NUM_THREADS];
	//for (int i = 0; i < NUM_THREADS; i++) {
		//tabs[i] = (char*)calloc(local_tab_size, sizeof(char));
	//}

#pragma omp parallel
	{
		//create local anserw list 
		int thread_num = omp_get_thread_num();
		//char* local_out = tabs[thread_num];

#pragma omp for schedule(static)
		for (numeric_type i = 1; i < primes_len; i++) { //for every *part of* primes

			//find first odd multiply of given prime in range <begin,b>
			numeric_type holder = divide_begin_odd % primes[i];
			holder = (holder == 0) ? 0 : primes[i] - holder;
			holder = divide_begin_odd + holder;
			holder = (holder % 2 == 0) ? holder + primes[i] : holder;

			numeric_type index = (holder - divide_begin_odd) / 2;

			//for (numeric_type x = holder; x <= b; x += primes[i]*2) { 
			for (numeric_type x = index; x < local_tab_size_real; x += primes[i]) { // for every number <begin,b> (no evens)
				final_out[x] = -1; //set to non prime;
			}
		}
		//after all connect outputs
			//for (numeric_type i = 0; i < local_tab_size_real; i++) {
				//if (local_out[i] == -1) {
					//final_out[i] = -1;
				//}
			//}
			//if (local_out) {
				//free(local_out);
			//
	}
	
	
	for (numeric_type i = 0; i < local_tab_size_real; i++) {
		if (final_out[i] == 0) {
			numeric_type number = (i * 2) + divide_begin_odd;
			outcome[(*outcome_len)++] = number;
		}
	}

	numeric_type min_outcome = (*outcome_len == 0) ? 1 : *outcome_len;	//max(outcome_len, 1); - make sure value is > 0
	outcome = (numeric_type*)realloc(outcome, min_outcome * sizeof(numeric_type));
	if (!outcome) {
		printf("Failure calling outcome = realloc()\n");
		return 0;
	}
	return outcome;
}

numeric_type* find_primes_in_range_concept_2_parallel_functional(numeric_type a, numeric_type b, numeric_type sieve_size,
	numeric_type* primes, numeric_type primes_len,
	numeric_type* outcome_len) {

	numeric_type finalSize = ((b - a) / 2) + 1;
	numeric_type* outcome = (numeric_type*)calloc(finalSize, sizeof(numeric_type));
	*outcome_len = 0;

	if (!outcome) {
		printf("Failure calling outcome = calloc(), probably out-of-memory.\n");
		return 0;
	}

	//copy arimes that are already calculated to ouput if needed;
	// next step evety odd number in between max(a, sieve_size) to b, if checking for modulo against all primes from the sieve;
	// some numbers (copied here) will make modulo itself, resulting in not being added to final solution
	for (numeric_type i = 0; i < primes_len; i++) {

		if (primes[i] >= a) {
			outcome[(*outcome_len)++] = primes[i];
		}
	}

	//step_by_step calculation to get the beginning point of diviation
	numeric_type max_sieve_value = sieve_size * 2;	//max value that has been already checked in the sieve
	numeric_type divide_begin = (a > max_sieve_value) ? a : max_sieve_value;	//where to begin; pick the higher from input_a and sieve_max_val
	numeric_type divide_begin_odd = (divide_begin % 2 == 1) ? divide_begin : divide_begin + 1;	//add 1 if needed; we skip even numbers here too

	char const NUM_THREADS = 8;
	omp_set_dynamic(0);
	omp_set_num_threads(NUM_THREADS);

	//prepare sieve size for odds in range <a,b>
	numeric_type local_tab_size = ((b - divide_begin_odd) / 2) + 1;
	//char* final_out = (char*)calloc(local_tab_size, sizeof(char));
	numeric_type local_tab_size_real = (b < divide_begin_odd || b == 1) ? 0 : local_tab_size; //to handle special cases

	//split bounds to threads
	numeric_type bounds[NUM_THREADS + 1];
	numeric_type bounds_numerics[NUM_THREADS + 1];
	numeric_type chunk_size = local_tab_size_real / NUM_THREADS;
	for (int i = 1; i < NUM_THREADS; i++) {
		bounds[i] = i * chunk_size;
		bounds_numerics[i] = divide_begin_odd + (i * chunk_size * 2);
	}
	bounds[0] = 0;
	bounds_numerics[0] = divide_begin_odd;
	bounds[NUM_THREADS] = local_tab_size_real;
	
	//create sieves for each thread
	char* tabs[NUM_THREADS];
	numeric_type tabs_sizes[NUM_THREADS];
	for (int i = 0; i < NUM_THREADS; i++) {
		tabs_sizes[i] = bounds[i + 1] - bounds[i];
		tabs[i] = (char*)calloc(tabs_sizes[i], sizeof(char));
	}

	char THREAD_FLAG = 0;

#pragma omp parallel
	{
		int thread_num = omp_get_thread_num();
		char* local_sieve = tabs[thread_num];
		numeric_type local_num_start = bounds_numerics[thread_num];
		numeric_type local_array_size = tabs_sizes[thread_num];

		for (numeric_type i = 1; i < primes_len; i++) { //for every prime

			//find first odd multiply of given prime in range <begin,b>
			numeric_type holder = local_num_start % primes[i];
			holder = (holder == 0) ? 0 : primes[i] - holder;
			holder = local_num_start + holder;
			holder = (holder % 2 == 0) ? holder + primes[i] : holder;

			numeric_type index = (holder - local_num_start) / 2;

			//#pragma omp parallel for schedule(static)
			for (numeric_type x = index; x < local_array_size; x += primes[i]) { // for every number <begin,b> (no evens)
				local_sieve[x] = -1; //set to non prime;
			}
		}

		while (THREAD_FLAG < thread_num) {
			Sleep(1);
		}
		for (numeric_type i = 0; i < local_array_size; i++) {
			if (local_sieve[i] == 0) {
				numeric_type number = (i * 2) + local_num_start;
				outcome[(*outcome_len)++] = number;
			}
		}
		free(local_sieve);
		THREAD_FLAG++;
	}

	numeric_type min_outcome = (*outcome_len == 0) ? 1 : *outcome_len;	//max(outcome_len, 1); - make sure value is > 0
	if (!realloc(outcome, min_outcome * sizeof(numeric_type))) {
		perror("Failure calling outcome = realloc()\n");
	}
	return outcome;
}



int main(int argc, char* argv[]) {

	/*
	mode 0: determina all primes using one BIG array NTAB
			- always count all prime numbers starting from 1 - filters them afterwards
			- memory consumption dependent from B parameter
				-> good if range(A, B) is wider
				-> linear performance and memory usage, dependent entirely on B value

	mode 1: make prime array NTAB only to sqrt(b), then split every element in range(a, b) by every prime to sqrt(b) 
			- counts primes only to sqrt(B), then divide every range(A, B) by every counted prime
			- memory consumption dependent from sqrt(B) + range(A, B)
				-> good use for tiny range(A, B) at high values
				-> it can target numbers that are out of memory for method 0
				-> the smaller range(A, B) is, the better it is both memory and performance-wise
				-> performance (and memory) dependent mainly on range(A, B), but also on B value
				-> if used for big range(A, B) it will become VERY slow
				-> less memory usage
	*/
	
	//clock_t start, stop;
	double start, stop;

	int mode; 
	numeric_type a, b;

	readInput(&a, &b, &mode);

	//start = clock();
	start = omp_get_wtime();

	numeric_type sieve_bottom_filter;	//index of a in sieve to filter out values below it
	numeric_type sieve_size;	//size of the sieve
	
	char *actual_sieve = prepare_sieve(a, b, mode, &sieve_size, &sieve_bottom_filter);
	if (!actual_sieve) { return 0; }

	eratostenes(actual_sieve, &sieve_size);
	
	if (mode == 0) { // mode 0 ends here
		//stop = clock();
		stop = omp_get_wtime();
		sieve_summary(actual_sieve, sieve_bottom_filter, sieve_size, &start, &stop);
		return 0;
	} // mode 1 continues

	//we will be looping through prime list multiple times, so lets rewrite it and free() the sieve
	numeric_type primes_len;
	numeric_type *primes = rewrite_primes(actual_sieve, sieve_size, &primes_len);
	if (!primes) { return 0; }

	free(actual_sieve);

	start = omp_get_wtime();  /// ------- START TIMER FOR MODE 1 !!!

	//fill and create output
	numeric_type outcome_len;
	numeric_type* outcome = find_primes_in_range_concept_2(a, b, sieve_size, primes, primes_len, &outcome_len);
	if (!outcome) { return -1; }

	//method 1 ends here
	stop = omp_get_wtime();
	printOutput(outcome, outcome_len, &start, &stop);
	return 0;
}