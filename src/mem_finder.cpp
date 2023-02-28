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
    int_t *LCP = NULL;
    LCP = (int_t*) malloc(n*sizeof(int_t));
    int32_t *DA = NULL;
    DA = (int32_t*) malloc(n*sizeof(int32_t));

    gsacak((unsigned char *)concat_data, (uint_t*)SA, LCP, DA, n);

    int_t min_mem_length = 30;

    std::vector<std::pair<uint_t, uint_t>> intervals = get_lcp_intervals(LCP, min_mem_length, n);
    
    std::vector<mem> mems;
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


