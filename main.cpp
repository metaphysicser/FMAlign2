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

// The main function is the entry point of the program. It is where the program starts executing. 
// the program starts executing. 
#include "include/common.h"
#include "include/utils.h"
#include "include/mem_finder.h"
#include "include/mem_filter.h"
#include "include/sequence_split_align.h"
#include "include/thread_pool.h"

#include <chrono>
#include <thread>

const char* data_path = "data/mt1x.fasta";

int main() {
    Timer timer;
    std::vector<std::string> data;
    std::vector<std::string> name;
    read_data(data_path, data, name);
    std::vector<mem> mems = find_mem(data);
    std::vector<mem> filtered_mems = filter_mem(mems);
    split_sequence(data, filtered_mems);
    double total_time = timer.elapsed_time();
    std::cout << "FMAlign2 total time: " << total_time << " seconds." << std::endl;
    return 0;
}