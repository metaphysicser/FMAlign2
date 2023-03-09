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
#include "include/sequence_split_align.h"
#include "include/thread_pool.h"
#include <thread>

GlobalArgs global_args;
int main(int argc, char** argv) {
    Timer timer;
    ArgParser parser;

    parser.add_argument("input", true, "data/mt1x.fasta");
    parser.add_argument_help("input", "The input file name.");
    parser.add_argument("thread", false, "cpu_num");
    parser.add_argument_help("thread", "The maximum number of threads that the program runs, the recommended setting is the number of CPUs");
    parser.add_argument("min_mem_length", false, "30");
    parser.add_argument_help("min_mem_length", "The minimum length of MEM");
    parser.add_argument("min_seq_coverage", false, "0.5");
    parser.add_argument_help("min_seq_coverage", "A floating-point parameter that specifies the minimum coverage across all sequences, with values ranging from 0 to 1.");


    try {
        parser.parse_args(argc, argv);
        global_args.data_path = parser.get("input");
        std::string tmp_thread = parser.get("thread");
        if (tmp_thread == "cpu_num") {
            global_args.thread = std::thread::hardware_concurrency();
        }
        else {
            global_args.thread = std::stoi(tmp_thread);
        }
        std::cout << "The number of threads is set to " << global_args.thread << std::endl;
        global_args.min_mem_length = std::stoi(parser.get("min_mem_length"));
        std::cout << "The minimum length of MEM is set to " << global_args.min_mem_length << std::endl;

        global_args.min_seq_coverage = std::stof(parser.get("min_seq_coverage"));
        if (global_args.min_seq_coverage < 0 || global_args.min_seq_coverage > 1) {
            std::cerr << "min_seq_coverage should be ranged from 0 to 1." << std::endl;
            exit(0);
        }
        else {
            std::cout << "min_seq_coverage is set to " << global_args.min_seq_coverage << std::endl;
        }
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        parser.print_help();
        return 1;
    }

    std::vector<std::string> data;
    std::vector<std::string> name;

    read_data(global_args.data_path.c_str(), data, name);

    std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence = find_mem(data);
    parallel_align(data, split_points_on_sequence);
    double total_time = timer.elapsed_time();
    std::cout << "FMAlign2 total time: " << total_time << " seconds." << std::endl;
    return 0;
}