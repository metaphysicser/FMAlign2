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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

#ifndef M64
	#define M64 0
#endif
// Considering that the concatenated input sequence can be too long and to save memory.
// different data types are defined through conditional compilation for different sizes of data.
// if data is larger than 2GB, M64 should be selected.
#if M64
	typedef int64_t	int_t;
	typedef uint64_t uint_t;
	#define U_MAX	UINT64_MAX
	#define I_MAX	INT64_MAX
	#define I_MIN	INT64_MIN
#else
	typedef int32_t int_t;
	typedef uint32_t uint_t;
	#define U_MAX	UINT32_MAX
	#define I_MAX	INT32_MAX
	#define I_MIN	INT32_MIN
#endif

#ifndef DEBUG
  #define DEBUG 0
#endif

struct GlobalArgs {
	std::string data_path;
};
extern GlobalArgs global_args;

struct sub_string{
    int_t sequence_index; // the sequence index that substring in
    uint_t position; // the begin position in the seqence
	uint_t* mem_index; // the unique index
};

struct mem{
    int_t mem_length; // substring length
	uint_t* mem_index; // the unique index
	float avg_pos = -1; // average position in sequences, initially set to -1
    std::vector<sub_string> substrings; // the substring set
};
#endif