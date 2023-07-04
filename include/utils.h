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
#include <algorithm>
#if (defined(__linux__))
    #include <unistd.h>
#else
    #include <io.h>
    #include <process.h>
#endif
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <iomanip>
#include <sstream>

#define TABLE_LEN 60
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
 * A command-line argument parser.
 * Usage:
 *   1. Create an ArgParser object.
 *   2. Call add_argument() for each command-line argument you want to parse.
 *   3. Call parse_args() to parse the command-line arguments.
 *   4. Use get() or has() to retrieve the values of the parsed arguments.
 *   5. Call print_help() to print a help message with the list of arguments.
 */
class ArgParser {
public:
    /**
     * Adds a command-line argument to the parser.
     *
     * @param name          The name of the argument, without the leading dashes.
     * @param required      Whether the argument is required or optional (default=false).
     * @param default_value The default value of the argument, if it is optional (default="").
     * @throws std::invalid_argument If an argument with the same name has already been added.
     */
    void add_argument(const std::string& name, bool required, const std::string& default_value);

    /**
     * Adds a help message for a command-line argument.
     *
     * @param name      The name of the argument, without the leading dashes.
     * @param help_text The help text for the argument.
     * @throws std::invalid_argument If an argument with the given name has not been added.
     */
    void add_argument_help(const std::string& name, const std::string& help_text);

    /**
     * Parses the command-line arguments.
     *
     * @param argc The number of command-line arguments.
     * @param argv The array of command-line argument strings.
     * @throws std::invalid_argument If an invalid or missing argument is encountered.
     */
    void parse_args(int argc, char** argv);

    /**
     * Returns the value of a parsed argument.
     *
     * @param name The name of the argument, without the leading dashes.
     * @return The value of the argument.
     * @throws std::invalid_argument If the argument is invalid or missing.
     */
    std::string get(const std::string& name) const;

    /**
     * Checks whether a parsed argument has been provided.
     *
     * @param name The name of the argument, without the leading dashes.
     * @return Whether the argument has been provided.
     */
    bool has(const std::string& name) const;

    /**
     * Prints a help message with the list of command-line arguments.
     */
    void print_help() const;

private:
    /**
     * A struct representing a command-line argument.
     */
    struct Arg {
        bool required;              // Whether the argument is required.
        std::string default_value;  // The default value of the argument, if it is optional.
        std::string value;          // The value of the argument, if it has been parsed.
        std::string help_text;      // The help text for the argument.
    };

    std::map<std::string, Arg> args_;  // The map of argument names to argument objects.
};

/**
 * @brief: read fasta and fastq format data
 * @param data_path   the path to the target data
 * @param data store sequence content
 * @param name store sequence name
 * @return multiple sequence stored in vector 
*/
void read_data(const char* data_path, std::vector<std::string>& data, std::vector<std::string>& name, bool verbose);

/**
 * @brief: Check whether the file exists in the specified path.
 * @param data_path   The file path to check.
 * @return Returns true if the file exists, otherwise false.
*/
bool access_file(const char* data_path);

/**
 * @brief Cleans the input sequence by removing any non-ATCG characters. 
 * This function removes any characters from the input sequence that are not A, T, C, or G (case-insensitive). 
 * The cleaned sequence is then returned as a new string.
 * @param sequence The input DNA sequence to be cleaned.
 * @return The cleaned DNA sequence as a new string.
*/
std::string clean_sequence(std::string sequence);

/**
* @brief Print information about the FMAlign2 algorithm
* This function prints various information about the FMAlign2 algorithm,
* including the mode (32-bit or 64-bit), number of threads, minimum MEM length,
* sequence coverage, and parallel align method.
* @return void
*/
void print_algorithm_info();

void print_table_line(const std::string& output);

void print_table_divider();

void print_table_bound();

#endif