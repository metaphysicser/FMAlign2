/*
 * Copyright [2023] [malab Pinglu Zhang]
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

#include "../include/mem_finder.h"

/**
 * Find MEMs in a set of sequences.
 *
 * @param data A vector of strings representing the sequences.
 * @return A vector of MEMs found in the sequences.
 */
std::vector<mem> find_mem(std::vector<std::string> data){
    uint_t n = 0;
    unsigned char* concat_data = concat_strings(data, n); 
    
    uint_t *SA = NULL;
    SA = (uint_t*) malloc(n*sizeof(uint_t));
    // LCP[0] = 0, LCP[i] = lcp(concat_data[SA[i]], concat_data[SA[i-1]])
    int_t *LCP = NULL;
    LCP = (int_t*) malloc(n*sizeof(int_t));
    int32_t *DA = NULL;
    DA = (int32_t*) malloc(n*sizeof(int32_t));

    gsacak((unsigned char *)concat_data, (uint_t*)SA, LCP, DA, n);
    
    int_t min_mem_length = 1000;
    int_t min_cross_sequence = 336;
    std::vector<uint_t> joined_sequence_bound;
    uint_t total_length = 0;
    for (uint_t i = 0; i < data.size(); i++) {
        joined_sequence_bound.push_back(total_length);
        total_length += data[i].length() + 1;
    }
    // Find all intervals with an LCP >= min_mem_length and <= min_cross_sequence
    std::vector<std::pair<uint_t, uint_t>> intervals = get_lcp_intervals(LCP, min_mem_length, min_cross_sequence, n);
    uint_t interval_size = intervals.size();
    std::vector<mem> mems;
    mems.resize(interval_size);
    uint_t threads = 8;
    // Convert each interval to a MEM in parallel
    IntervalToMemConversionParams* params = new IntervalToMemConversionParams[interval_size];
    threadpool pool;
    threadpool_init(&pool, threads);
    for (uint_t i = 0; i < interval_size; i++) {
        params[i].SA = SA;
        params[i].DA = DA;
        params[i].interval = intervals[i];
        params[i].concat_data = concat_data;
        params[i].result_store = mems.begin() + i;
        params[i].min_mem_length = min_mem_length;
        params[i].joined_sequence_bound = joined_sequence_bound;

        threadpool_add_task(&pool, interval2mem, params+i);
    }
    threadpool_destroy(&pool);
    // Sort the MEMs based on their average positions and assign their indices
    sort_mem(mems, data);

    return mems;
}

/**
 * @brief Concatenates a vector of strings with separator 1 and a terminating 0.
 * @param strings The vector of strings to concatenate.
 * @param n A reference to the total length of the concatenated string.
 * @return A pointer to the concatenated string.
 * @note The returned string must be deleted by the caller.
*/
unsigned char* concat_strings(const std::vector<std::string>& strings, uint_t &n) {
     // Calculate total length of concatenated string
    uint_t total_length = std::accumulate(strings.begin(), strings.end(), 0,
                                          [](uint_t sum, const std::string& s) { return sum + s.length() + 1; });
    total_length++;  // Add 1 for the terminating 0

    // Allocate memory for concatenated string
    unsigned char* concat_data = new unsigned char[total_length];

    // Concatenate all strings with 1 as separator
    uint_t index = 0;
    for (const auto& s : strings) {
        std::copy(s.begin(), s.end(), concat_data + index);
        index += s.length();
        concat_data[index] = 1;
        index++;
    }

    // Set the terminating 0
    concat_data[total_length - 1] = 0;

    n = total_length;

    return concat_data;

}

/**
 * @brief an LCP (Longest Common Prefix) array and a threshold value,
 * finds all the LCP intervals where each value is greater than or equal to the threshold value,
 * and at least one value in the interval is equal to the threshold value.
 * @param lcp_array The input LCP array
 * @param threshold The threshold value
 * @param min_cross_sequence the min number of crossed sequence
 * @return  The output vector of pairs representing the LCP intervals
*/
std::vector<std::pair<uint_t, uint_t>> get_lcp_intervals(int_t* lcp_array, int_t threshold, int_t min_cross_sequence, uint_t n) {

    std::vector<std::pair<uint_t, uint_t>> intervals;

    int_t left = 0, right = 0;
    bool found = false;

    while (right < (int_t)n) {
        if (lcp_array[right] >= threshold) {
            if (lcp_array[right] == threshold) {
                found = true;
            }

            right++;
        } else {
            if (found && right-left+1 >= min_cross_sequence) {
                intervals.emplace_back(left, right);
            }

            left = right = right + 1;
            found = false;
        }
    }

    if (found && right - left + 1 >= min_cross_sequence) {
        intervals.emplace_back(left, right);
    }
    return intervals;
}

void draw_lcp_curve(int_t *LCP, uint_t n){
    std::ofstream outfile("tmp/lcp.bin", std::ios::out | std::ios::binary);
    
    // Write the number of elements in the array
    outfile.write(reinterpret_cast<const char*>(&n), sizeof(n));

    // Write the LCP array
    outfile.write(reinterpret_cast<const char*>(LCP), n * sizeof(int_t));

    outfile.close();

    // corresponding python code
    // import struct
    // import numpy as np
    // import matplotlib.pyplot as plt

    // def read_lcp_array(filename):
    //     with open(filename, "rb") as f:
    //         # Read the number of elements in the array
    //         n_bytes = f.read(4)
    //         n = struct.unpack("I", n_bytes)[0]

    //         # Read the LCP array
    //         lcp_bytes = f.read(n * 4)
    //         lcp_array = np.frombuffer(lcp_bytes, dtype=np.int32)

    //         return lcp_array

    // # Example usage
    // lcp_array = read_lcp_array("tmp/lcp.bin")

    // # Plot the LCP array histogram
    // plt.hist(lcp_array, bins=100)
    // plt.show()

}

/**
*@brief This function converts an LCP interval to a MEM (Maximal Exact Match).
*@param arg A void pointer to the input parameters.
*@return void* A void pointer to the result, which is stored in the input parameters structure.
*/
void* interval2mem(void* arg) {
    // Cast the input parameters to the correct struct type
    IntervalToMemConversionParams* ptr = static_cast<IntervalToMemConversionParams*>(arg);
    // Extract the necessary variables from the struct
    const uint_t* SA = ptr->SA;
    const int32_t* DA = ptr->DA;
    const int_t min_mem_length = ptr->min_mem_length;
    const unsigned char* concat_data = ptr->concat_data;
    const std::vector<uint_t> joined_sequence_bound = ptr->joined_sequence_bound;
    // Initialize the result variables
    std::pair<uint_t, uint_t> interval = ptr->interval;
    mem result;
    uint_t* mem_index = new uint_t;
    *mem_index = 0;
    result.mem_index = mem_index;
    std::vector<sub_string> res_substrings;
    result.mem_length = 0;
    std::vector<uint_t> mem_position;
    // Create the MEM from the input LCP interval
    for (uint_t i = interval.first - 1; i < interval.second; i++) {
        sub_string tmp_substring;
        tmp_substring.sequence_index = DA[i];
        tmp_substring.position = SA[i] - joined_sequence_bound[tmp_substring.sequence_index];
        mem_position.push_back(SA[i]);
        
        tmp_substring.mem_index = mem_index;
        result.substrings.push_back(tmp_substring);
    }
    // Compute the offset of the MEM and adjust the positions of the substrings accordingly
    uint_t offset = 1;
    bool all_char_same = true;
    while (true) {
        char current_char = 0;
        if (mem_position[0] >= offset) {
            current_char = concat_data[mem_position[0] - offset];
        }
        else {
            break;
        }
        all_char_same = true;
        for (uint_t i = 1; i < mem_position.size(); i++) {
            if (mem_position[i] < offset) {
                all_char_same = false;
                break;
            }
            else {
                if (current_char != concat_data[mem_position[i] - offset]) {
                    all_char_same = false;
                    break;
                }
            }

        }
        if (all_char_same == false) {
            break;
        }
        offset++;
    }
    offset -= 1;

    result.mem_length = min_mem_length + offset;
    for (uint_t i = 0; i < result.substrings.size(); i++) {
        result.substrings[i].position -= offset;
    }
    // Store the result in the input parameters structure
    *(ptr->result_store) = result;
  
    return NULL;
}

// function to compute average position of a mem in sequences
void compute_mem_avg_pos(mem& m) {
    float_t sum_pos = 0;
    for (auto& s : m.substrings) {
        sum_pos += s.position;
    }
    m.avg_pos = sum_pos / m.substrings.size();
}

/**
*Sorts the input vector of MEMs by the average position of each MEM's substrings along the sequences.
*Removes any MEMs that span across multiple sequences.
*Assigns a unique index to each MEM based on its position in the sorted vector.
*@param mems The vector of MEMs to be sorted.
*@param data The vector of sequences used to compute the MEMs.
*/
void sort_mem(std::vector<mem> &mems, std::vector<std::string> data) {
    
    auto it = std::remove_if(mems.begin(), mems.end(), [&](const mem& m) {
        if (m.substrings[0].position + m.mem_length < data[m.substrings[0].sequence_index].length()) {
            return false;
        }
        else {
            return true;
        }
       
        });
    mems.erase(it, mems.end());

    // compute average position of each mem
    for (auto& m : mems) {
        compute_mem_avg_pos(m);
    }
    // sort mems by average position
    std::sort(mems.begin(), mems.end(), [](const mem& m1, const mem& m2) {
        return m1.avg_pos < m2.avg_pos;
        });
    // assign mem_index based on position in sorted vector
    for (uint_t i = 0; i < mems.size(); i++) {
        *mems[i].mem_index = i;
    }
    return;
}


