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
// Created: 2023-02-24

//This file defines data types, data structures, and various parameters.
#ifndef COMMON_H
#define COMMON_H

#include <vector>

// Considering that the concatenated input sequence can be too long and to save memory
// different data types are defined through conditional compilation for different sizes of data.
typedef unsigned int joined_data_type;

// This is the default value of mem minimum length and sequence count min proportion
// The actual value used depends on user input
unsigned int mem_minimum_length = 39;
float sequence_count_min_proportion = 0.5;

// To store large string data more efficiently in memory, choose vector<char> instead of string.
std::vector<std::vector<char>> data;

// Intermediate files are stored in the tmp folder
const char* tmp_path = "tmp/";

#endif