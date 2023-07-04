/*
 * Copyright [2023] [MALABZ_UESTC Pinglu Zhang]
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Author: Pinglu Zhang
// Contact: zpl010720@gmail.com
// Created: 2023-02-25

#ifndef MEM_FINDER_H
#define MEM_FINDER_H

#include "common.h"
#include "gsacak.h"
#include "utils.h"
#include <cstdint>
#include <cstring>
#include <numeric>
#ifdef __linux__
    #include <omp.h>
#endif
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <sstream>
struct sub_string {
    int_t sequence_index; // the sequence index that substring in
    uint_t position; // the begin position in the seqence
    uint_t* mem_index; // the unique index
};

struct mem {
    int_t mem_length; // substring length
    uint_t* mem_index; // the unique index
    float avg_pos = -1; // average position in sequences, initially set to -1
    std::vector<sub_string> substrings; // the substring set
};

struct IntervalToMemConversionParams {
    const uint_t* SA;
    const int32_t* DA;
    const unsigned char* concat_data;
    std::vector<mem>::iterator result_store;
    int_t min_mem_length;
    std::pair<uint_t, uint_t> interval;
    std::vector<uint_t> joined_sequence_bound;
};

struct FindOptimalChainParams {
    std::vector<std::vector<std::pair<int_t, int_t>>>::iterator chains;
};

/**
* @brief DP Only Once!Filter out overlapping memory regions and generate split points for each sequence.
* Given a vector of memory regions and the number of sequences, this function removes any
* overlapping memory regions and generates split points for each sequence based on the non-overlapping regions.
* @param mems Vector of memory regions.
* @param sequence_num Number of sequences.
* @return Vector of split points for each sequence.
*/
std::vector<std::vector<std::pair<int_t, int_t>>> filter_mem_fast(std::vector<mem>& mems, uint_t sequence_num);

/**
* @brief DP sequence number times!Filter out overlapping memory regions and generate split points for each sequence.
* Given a vector of memory regions and the number of sequences, this function removes any
* overlapping memory regions and generates split points for each sequence based on the non-overlapping regions.
* @param mems Vector of memory regions.
* @param sequence_num Number of sequences.
* @return Vector of split points for each sequence.
*/
std::vector<std::vector<std::pair<int_t, int_t>>> filter_mem_accurate(std::vector<mem>& mems, uint_t sequence_num);

/**
 * @brief Find MEMs in a set of sequences.
 * @param data A vector of strings representing the sequences.
 * @return Vector of split points for each sequence.
 */
std::vector<std::vector<std::pair<int_t, int_t>>> find_mem(std::vector<std::string> data);

/**
 * @brief Concatenates a vector of strings with separator 1 and a terminating 0.
 * @param strings The vector of strings to concatenate.
 * @param n A reference to the total length of the concatenated string.
 * @return A pointer to the concatenated string.
 * @note The returned string must be deleted by the caller.
*/
unsigned char* concat_strings(const std::vector<std::string>& strings, uint_t &n);

/**
 * @brief an LCP (Longest Common Prefix) array and a threshold value,
 * finds all the LCP intervals where each value is greater than or equal to the threshold value,
 * and at least one value in the interval is equal to the threshold value.
 * @param lcp_array The input LCP array
 * @param threshold The threshold value
 * @param min_cross_sequence the min number of crossed sequence
 * @return  The output vector of pairs representing the LCP intervals
*/
std::vector<std::pair<uint_t, uint_t>> get_lcp_intervals(int_t* lcp_array, int_t threshold, int_t min_cross_sequence, uint_t n);

/**
*@brief This function converts an LCP interval to a MEM (Maximal Exact Match).
*@param arg A void pointer to the input parameters.
*@return void* A void pointer to the result, which is stored in the input parameters structure.
*/
void* interval2mem(void* arg);

/**
*Sorts the input vector of MEMs by the average position of each MEM's substrings along the sequences.
*Removes any MEMs that span across multiple sequences.
*Assigns a unique index to each MEM based on its position in the sorted vector.
*@param mems The vector of MEMs to be sorted.
*@param data The vector of sequences used to compute the MEMs.
*/
void sort_mem(std::vector<mem>& mems, std::vector<std::string> data);
#endif