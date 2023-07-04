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


#include "../include/utils.h"

KSEQ_INIT(int, read)

/**
 * @brief A timer class that measures elapsed time. 
 * This class uses C++11 chrono library to measure elapsed time in seconds with double precision. 
 * The timer starts at construction and can be reset to zero by calling reset(). 
 * The elapsed time can be obtained by calling elapsed_time() method. 
 * The timer is based on std::chrono::steady_clock, which is a monotonic clock that is not subject to system clock adjustments.
*/
Timer::Timer() {
    start_time_ = std::chrono::steady_clock::now();
}
void Timer::reset() {
    start_time_ = std::chrono::steady_clock::now();
}

double Timer::elapsed_time() const {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time_;
    return elapsed.count();
}

/**
 * @brief: read fasta and fastq format data
 * @param data_path   the path to the target data
 * @param data store sequence content
 * @param name store sequence name
 * @return multiple sequence stored in vector 
*/
void read_data(const char* data_path, std::vector<std::string>& data, std::vector<std::string>& name, bool verbose = true){
    if (verbose && global_args.verbose) {
        std::cout << "#                   Reading Data...                         #" << std::endl;
        print_table_divider();
    }
    std::string output = "";
    std::string str_data_path = data_path;
    // check weather the input path could be accessed 

    if (access_file(data_path)) {
        if (verbose && global_args.verbose) {
            output = str_data_path + " could be accessed";
            print_table_line(output);
        }
    }
    else {
        print_table_bound();
        output = "Error:" + str_data_path + " could not be accessed, Please check if the path of the input data is correct or if the data exists!";
        std::cerr << output << std::endl;
        std::cerr << "Program Exit!" << std::endl;
        exit(1);
    }
    

    FILE* f_pointer = fopen(data_path, "r");
    kseq_t* file_t = kseq_init(fileno(f_pointer));
    
    uint64_t merged_length = 0;
    int64_t tmp_length = 0; 
    // stop loop when tmp_length equals -1
    while ((tmp_length = kseq_read(file_t)) >= 0) // Read one sequence in each iteration of the loop
    {
        std::string tmp_data = clean_sequence(file_t -> seq.s);
        std::string tmp_name = file_t -> name.s;
        if(file_t->comment.s) tmp_name += file_t->comment.s;
        data.push_back(tmp_data);
        name.push_back(tmp_name);
        merged_length += tmp_length;
    }
    kseq_destroy(file_t);
    fclose(f_pointer);

    if(verbose&& global_args.verbose && merged_length + data.size() > UINT32_MAX && M64 == 0){
        print_table_bound();
        std::cerr << "Error: The input data is too large and the 32-bit program may not produce correct results. Please compile a 64-bit program using the M64 parameter." << std::endl;
        std::cerr << "Program Exit!" << std::endl;
        exit(1);
    }
    #if M64
    if (verbose && global_args.verbose) {
        std::stringstream s;
        s << std::fixed << std::setprecision(2) << merged_length / pow(2, 30);
        output = "Data Memory Usage: " + s.str() + " GB";
        print_table_line(output);
    }
    
    #else
    if (verbose && global_args.verbose) {
        std::stringstream s;
        s << std::fixed << std::setprecision(2) << merged_length / pow(2, 20);
        output = "Data Memory Usage: " + s.str() + " MB";
        print_table_line(output);
    }
    #endif
    if (verbose && global_args.verbose) {
        output = "Sequence Number: " + std::to_string(data.size());
        print_table_line(output);
        print_table_divider();
    }
    return;
}

/**
 * @brief: Check whether the file exists in the specified path.
 * @param data_path   The file path to check.
 * @return Returns true if the file exists, otherwise false.
*/
bool access_file(const char* data_path){
    std::ifstream file(data_path);
    if (file.is_open() && file.good()) {
        return true;
    } else {
        return false;
    }
    return true;
}

/**
 * @brief Cleans the input sequence by removing any non-ATCG characters. 
 * This function removes any characters from the input sequence that are not A, T, C, or G (case-insensitive). 
 * The cleaned sequence is then returned as a new string.
 * @param sequence The input DNA sequence to be cleaned.
 * @return The cleaned DNA sequence as a new string.
*/
std::string clean_sequence(std::string sequence){
    std::string result;
    result.reserve(sequence.size());
    for (char& c : sequence) {
        if (c == 'a' || c == 'A') {
            c = 'A';
            result.push_back(c);
        } else if (c == 'c' || c == 'C') {
            c = 'C';
            result.push_back(c);
        } else if (c == 'g' || c == 'G') {
            c = 'G';
            result.push_back(c);
        } else if (c == 't' || c == 'T') {
            c = 'T';
            result.push_back(c);
        }
        else if (c == 'u' || c == 'U') {
            c = 'U';
            result.push_back(c);
        }
        else {
            c = '-';
            result.push_back(c);
        }
    } 
    return result;
}

void ArgParser::add_argument(const std::string& name, bool required = false, const std::string& default_value = "") {
    if (args_.count(name) > 0) {
        throw std::invalid_argument("Duplicate argument name: " + name);
    }
    args_[name] = { required, default_value };
}

void ArgParser::add_argument_help(const std::string& name, const std::string& help_text) {
    if (args_.count(name) == 0) {
        throw std::invalid_argument("Invalid argument name: " + name);
    }
    args_[name].help_text = help_text;
}

void ArgParser::print_help() const {
    std::cout << "Usage: FMAlign2 [OPTIONS]\n\n";
    std::cout << "Options:\n";

    for (const auto& p : args_) {
        std::cout << "  -" << p.first;

        if (!p.second.required) {
            std::cout << " (optional)";
        }

        if (!p.second.default_value.empty()) {
            std::cout << " [default: " << p.second.default_value << "]";
        }

        std::cout << "\n    " << p.second.help_text << "\n\n";
    }
}

void ArgParser::parse_args(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg[0] == '-') {
            std::string name = arg.substr(1); 
            if (name == "help" || name == "h") {
                print_help();
                exit(0);
            }
            if (args_.count(name) == 0) {
                throw std::invalid_argument("Invalid argument: " + arg);
            }
            Arg& a = args_[name];

            if (a.value != "") {
                throw std::invalid_argument("Duplicate argument: " + arg);
            }
            if (i == argc - 1 || argv[i + 1][0] == '-') {
                if (a.required) {
                    throw std::invalid_argument("Missing value for argument: " + arg);
                }
                a.value = a.default_value;
            }
            else {
                a.value = argv[++i];
            }
        }
    }
    for (auto& p : args_) {
        if (p.second.required && p.second.value == "") {
            throw std::invalid_argument("Missing required argument: -" + p.first);  // 修改这里，改为单个破折号
        }
        if (p.second.required == false && p.second.value == "") {
            p.second.value = p.second.default_value;
        }
    }
}


std::string ArgParser::get(const std::string& name) const {
    if (args_.count(name) == 0) {
        throw std::invalid_argument("Invalid argument name: " + name);
    }
    const Arg& a = args_.at(name);
    if (a.value == "") {
        throw std::invalid_argument("Missing value for argument: --" + name);
    }
    return a.value;
}

bool ArgParser::has(const std::string& name) const {
    return args_.count(name) > 0 && args_.at(name).value != "";
}

/**
* @brief Print information about the FMAlign2 algorithm
* This function prints various information about the FMAlign2 algorithm,
* including the mode (32-bit or 64-bit), number of threads, minimum MEM length,
* sequence coverage, and parallel align method.
* @return void
*/
void print_algorithm_info() {
    print_table_bound();
    std::cout << "#               FMAlign2 algorithm info                     #" << std::endl;
    print_table_divider();
#if M64
    std::string output = "Mode: 64 bit";
    print_table_line(output);
#else
    std::string output = "Mode: 32 bit";
    print_table_line(output);
#endif
    std::string thread_output = "Thread: " + std::to_string(global_args.thread);
    print_table_line(thread_output);

    std::string l_output;
    if (global_args.min_mem_length < 0) {
        l_output = "Minimum MEM length: square root of mean length";
    }
    else {
        l_output = "Minimum MEM length: " + std::to_string(global_args.min_mem_length);
    }
    
    print_table_line(l_output);

    std::stringstream s;
    s << std::fixed << std::setprecision(2) << global_args.min_seq_coverage;

    std::string c_output;
    if (global_args.min_seq_coverage < 0) {
        c_output = "Sequence coverage: default";
    }
    else {
        c_output = "Sequence coverage: " + std::to_string(global_args.min_mem_length);
    }

    print_table_line(c_output);

    std::string p_output = "Parallel align method: " + global_args.package;
    print_table_line(p_output);

    print_table_bound();
}

void print_table_line(const std::string &output) {
    std::cout << "# " << std::left << std::setw(TABLE_LEN-2) << std::setfill(' ') << output << "#" << std::endl;
}

void print_table_divider() {
    std::cout << "#" << std::right << std::setw(TABLE_LEN) << std::setfill('-') << "#" << std::endl;
}

void print_table_bound() {
    std::cout << "#" << std::right << std::setw(TABLE_LEN) << std::setfill('#') << "#" << std::endl;
}


