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
#if defined(__linux__)
#include "include/thread_pool.h"
#endif
#include <thread>

GlobalArgs global_args;
int main(int argc, char** argv) {
    // Create a Timer object to record the execution time.
    Timer timer;
    // Create an ArgParser object to parse command line arguments.
    ArgParser parser;
    std::string output = "";
    // Add command line arguments to the ArgParser object.
    parser.add_argument("in", true, "data/mt1x.fasta");
    parser.add_argument_help("in", "The path to the input file.");
    parser.add_argument("t", false, "cpu_num");
    parser.add_argument_help("t", "The maximum number of threads that the program runs, the recommended setting is the number of CPUs.");
    parser.add_argument("l", false, "square root of mean length");
    parser.add_argument_help("l", "The minimum length of MEM.");
    parser.add_argument("c", false, "0.7");
    parser.add_argument_help("c", "A floating-point parameter that specifies the minimum coverage across all sequences, with values ranging from 0 to 1.");
    parser.add_argument("p", false, "halign");
    parser.add_argument_help("p", "The MSA method used in parallel align. for example, halign, mafft and so on.");
    parser.add_argument("o", false, "output.aligned.fasta");
    parser.add_argument_help("o", "The path to the output file.");
    parser.add_argument("d", false, "0");
    parser.add_argument_help("o", "depth of recursion");

    // Add command line arguments to the ArgParser object.
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
        std::string tmp_len = parser.get("l");
        if (tmp_len != "square root of mean length") {
            global_args.min_mem_length = std::stoi(parser.get("l"));
        }
        else {
            global_args.min_mem_length = -1;
        }
        

        global_args.degree = std::stoi(parser.get("d"));
        if (global_args.degree > 1) {
            exit(1);
        }
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

        global_args.output_path = parser.get("o");
    } // Catch any invalid arguments and print the help message.
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        parser.print_help();
        return 1;
    }

    print_algorithm_info();

    std::vector<std::string> data;
    std::vector<std::string> name;

    try {
        // Read data from the input file and store in data and name vectors
        read_data(global_args.data_path.c_str(), data, name, true);
        // Find MEMs in the sequences and split the sequences into fragments for parallel alignment.
        std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence = find_mem(data);
        split_and_parallel_align(data, name, split_points_on_sequence);
    }
    catch (const std::bad_alloc& e) { // Catch any bad allocations and print an error message.
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