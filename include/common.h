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
// Created: 2023-02-24

//This file defines data types, data structures, and various parameters.
#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#if (defined(__linux__))
#include "thread_pool.h"
#else
#include <omp.h>
#include <windows.h>
#endif
#include <iomanip>
#include <thread>

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
	int_t thread;
	int_t min_mem_length;
	float min_seq_coverage;
	std::string package;
	std::string output_path;
	int_t degree;
	std::string filter_mode;
	int_t verbose;
	double avg_file_size;
};
extern GlobalArgs global_args;


#endif