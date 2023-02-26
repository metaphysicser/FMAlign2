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

// This header file defines utility functions for reading and outputting data and so on.
// These utility functions are designed to be reusable across different parts of the project and can be easily included in other source files.
#ifndef UTILS_H
#define UTILS_H

#include "common.h"
#include "kseq.h"
#include <fstream>
#include <iomanip> 
#include <chrono> 
#include <math.h>
#if (defined(__linux__) || defined(__APPLE__))
#include <unistd.h>
#else
#include <io.h>
#include <process.h>
#endif

/**
 * @brief A timer class that measures elapsed time. 
 * This class uses C++11 chrono library to measure elapsed time in seconds with double precision. 
 * The timer starts at construction and can be reset to zero by calling reset(). 
 * The elapsed time can be obtained by calling elapsed_time() method. 
 * The timer is based on std::chrono::steady_clock, which is a monotonic clock that is not subject to system clock adjustments.
*/
class Timer{
public:
    Timer();
    void reset();
    double elapsed_time() const;
private:
    std::chrono::time_point<std::chrono::steady_clock> start_time_;
};

/**
 * @brief: read fasta and fastq format data
 * @param data_path   the path to the target data
 * @param data store sequence content
 * @param name store sequence name
 * @return multiple sequence stored in vector 
*/
void read_data(const char* data_path, std::vector<std::string>& data, std::vector<std::string>& name);

/**
 * @brief: Check whether the file exists in the specified path.
 * @param data_path   The file path to check.
 * @return Returns true if the file exists, otherwise false.
*/
bool access_file(const char* data_path);

#endif