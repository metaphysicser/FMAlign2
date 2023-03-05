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
    
    int_t min_mem_length = 2000;

    
    std::vector<std::pair<uint_t, uint_t>> intervals = get_lcp_intervals(LCP, min_mem_length, n);
    uint_t interval_size = intervals.size();
    std::vector<mem> mems;
    mems.resize(interval_size);
    uint_t threads = 8;

    IntervalToMemConversionParams* params = new IntervalToMemConversionParams[interval_size];
    threadpool pool;
    threadpool_init(&pool, threads);
    for (uint_t i = 0; i < interval_size; i++) {
        params[i].SA = SA;
        params[i].LCP = LCP;
        params[i].DA = DA;
        params[i].interval = intervals[i];
        params[i].concat_data = concat_data;
        params[i].result_store = mems.begin() + i;
        params[i].min_mem_length = min_mem_length;

        threadpool_add_task(&pool, interval2mem, params+i);
    }
    threadpool_destroy(&pool);

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
 * @return  The output vector of pairs representing the LCP intervals
*/
std::vector<std::pair<uint_t, uint_t>> get_lcp_intervals(int_t* lcp_array, int_t threshold, uint_t n) {

    std::vector<std::pair<uint_t, uint_t>> intervals;

    uint_t left = 0, right = 0;
    bool found = false;

    while (right < n) {
        if (lcp_array[right] >= threshold) {
            if (lcp_array[right] == threshold) {
                found = true;
            }

            right++;
        } else {
            if (found) {
                intervals.emplace_back(left, right);
            }

            left = right = right + 1;
            found = false;
        }
    }

    if (found) {
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

void* interval2mem(void* arg) {
    IntervalToMemConversionParams* ptr = static_cast<IntervalToMemConversionParams*>(arg);
    const uint_t* SA = ptr->SA;
    const int_t* LCP = ptr->LCP;
    const int32_t* DA = ptr->DA;
    const int_t min_mem_length = ptr->min_mem_length;
    const unsigned char* concat_data = ptr->concat_data;
    
    std::pair<uint_t, uint_t> interval = ptr->interval;
    mem result;
    std::vector<sub_string> res_substrings;
    result.mem_lengh = 0;
    std::vector<uint_t> mem_position;

    for (uint_t i = interval.first-1; i <= interval.second; i++) {
        sub_string tmp_substring;
        tmp_substring.position = SA[i];
        mem_position.push_back(tmp_substring.position);
        tmp_substring.sequence_index = DA[SA[i]];
        result.substrings.push_back(tmp_substring);
    }

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

    result.mem_lengh = min_mem_length + offset;
    for (uint_t i = 0; i < result.substrings.size(); i++) {
        result.substrings[i].position -= offset;
    }
    *(ptr->result_store) = result;

   
    return NULL;
}


