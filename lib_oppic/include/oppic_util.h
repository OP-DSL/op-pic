
/* 
BSD 3-Clause License

Copyright (c) 2022, OP-DSL

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <sys/time.h>

typedef struct oppic_set_core *oppic_set;

#ifndef MIN
#define MIN(a, b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) ((a > b) ? (a) : (b))
#endif

//********************************************************************************
std::string getTimeStr();

//********************************************************************************
inline char *copy_str(char const *src) 
{
    size_t src_len = strlen(src); // Calculate the actual length of src
    size_t dest_len = (src_len > 100) ? 100 : src_len; // Limit the destination length to 100 characters

    char *dest = (char *)malloc((dest_len + 1) * sizeof(char));
    if (dest) {
        memcpy(dest, src, dest_len);
        dest[dest_len] = '\0'; // Ensure the destination string is null-terminated
    }
    return dest;
}

//********************************************************************************
std::vector<size_t> sort_indexes(const int* cell_indices, int size);

//********************************************************************************

int min(int array[], int size);

int binary_search(int a[], int value, int low, int high);

int linear_search(int a[], int value, int low, int high);

void quickSort(int arr[], int left, int right);

void quickSort_2(int arr1[], int arr2[], int left, int right);

void quickSort_dat(int arr[], char dat[], int left, int right, int elem_size);

void quickSort_map(int arr[], int map[], int left, int right, int dim);

int removeDups(int a[], int array_size);

int compare_sets(oppic_set set1, oppic_set set2);

void op_timers(double *cpu, double *et); 

template <typename T>
std::string str(T value, std::string format) {
    char buffer[32];
    //int numChars = snprintf(buffer, sizeof(buffer), format.c_str(), value);
    snprintf(buffer, sizeof(buffer), format.c_str(), value);
    return std::string(buffer);
}

double rnd();

void reset_seed();

double* get_dandom_distriution(int count, int dim);

int file_exist(char const *filename);

bool opp_type_equivalence(const char *a, const char *b);