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
    std::string output = "";

    parser.add_argument("in", true, "data/mt1x.fasta");
    parser.add_argument_help("in", "The input file name.");
    parser.add_argument("t", false, "cpu_num");
    parser.add_argument_help("t", "The maximum number of threads that the program runs, the recommended setting is the number of CPUs");
    parser.add_argument("l", false, "30");
    parser.add_argument_help("l", "The minimum length of MEM");
    parser.add_argument("c", false, "0.5");
    parser.add_argument_help("c", "A floating-point parameter that specifies the minimum coverage across all sequences, with values ranging from 0 to 1.");
    parser.add_argument("p", false, "halign");
    parser.add_argument_help("p", "The MSA method used in parallel align. for example, halign, mafft and so on.");
    parser.add_argument("o", false, "output.aligned.fasta");
    parser.add_argument_help("o", "The output file name with its path");


    try {
        parser.parse_args(argc, argv);
        global_args.data_path = parser.get("in");
        std::string tmp_thread = parser.get("t");
        if (tmp_thread == "cpu_num") {
            global_args.thread = std::thread::hardware_concurrency();
        }
        else {
            global_args.thread = std::stoi(tmp_thread);
        }
      
        global_args.min_mem_length = std::stoi(parser.get("l"));      

        global_args.min_seq_coverage = std::stof(parser.get("c"));
        if (global_args.min_seq_coverage < 0 || global_args.min_seq_coverage > 1) {
            std::string output = "Error: min_seq_coverage should be ranged from 0 to 1!";
            std::cerr << output << std::endl;
            std::cerr << "Program Exit!" << std::endl;
            exit(1);
        }
        global_args.package = parser.get("p");
        if (global_args.package != "halign" && global_args.package != "mafft") {
            std::string output = "Error: " + global_args.package + " is a invalid method!";
            std::cerr << output << std::endl;
            std::cerr << "Program Exit!" << std::endl;
            exit(1);
        }
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        parser.print_help();
        return 1;
    }

    print_algorithm_info();

    std::vector<std::string> data;
    std::vector<std::string> name;

    try {
        read_data(global_args.data_path.c_str(), data, name, true);
        std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence = find_mem(data);
        split_and_parallel_align(data, name, split_points_on_sequence);
    }
    catch (const std::bad_alloc& e) {
        print_table_bound();
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << "Program Exit!" << std::endl;
        exit(1);
    }

    double total_time = timer.elapsed_time();
    std::stringstream s;
    s << std::fixed << std::setprecision(2) << total_time;
    output = "FMAlign2 total time: " + s.str() + " seconds.";
    print_table_line(output);
    print_table_bound();
    return 0;
}