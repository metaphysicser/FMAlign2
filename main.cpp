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
#include <filesystem>
namespace fs = std::filesystem;

#include <fstream>
#include <sstream>
#include <string>

std::string readFile(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) {
        // 打不开文件时，直接返回原字符串（说明它可能本身就是模板）
        return path;
    }
    std::ostringstream buf;
    buf << ifs.rdbuf();
    return buf.str();
}

#include <string>
#include <stdexcept>



int test_cmd(const std::string& cmdTemplate, int thread = 1) {
    const fs::path tempDir = "./temp";
    const fs::path inPath = tempDir / "tiny.fasta";
    const fs::path outPath = tempDir / "aligned.fasta";

    int rc = 0;  // 最终返回值（命令退出码或我们的错误码）

    // 用 do{...}while(false) 做“单出口”
    do {
        // 1) 创建 ./temp
        std::error_code ec;
        fs::create_directories(tempDir, ec);
        if (ec) {
            std::cerr << "Failed to create " << tempDir << ": " << ec.message() << "\n";
            rc = 1; break;
        }

        // 2) 写一个很小的 FASTA
        {
            std::ofstream ofs(inPath);
            if (!ofs) {
                std::cerr << "Cannot open " << inPath << " for writing.\n";
                rc = 1; break;
            }
            ofs <<
                R"(>seq1
ACGTACGTGA
>seq2
ACGTTGCA
>seq3
ACGTACGA
)";
        }

        // 3) 生成命令
        std::string cmd = buildCommand(cmdTemplate, inPath.string(), outPath.string(), thread);
        std::cout << "Running: " << cmd << "\n";

        // 4) 执行（注意 system 返回的是 wait 状态码，简单起见直接用）
        rc = std::system(cmd.c_str());

        // 5) 打印结果状态 & 预览输出
        if (rc == 0) {
            std::cout << "cmd finished successfully.\n";
        }
        else {
            std::cout << "cmd failed with code: " << rc << "\n";
        }

        if (fs::exists(outPath)) {
            std::ifstream ifs(outPath);
            std::string line; int cnt = 0;
          
        }
        else {
            std::cout << "Output file not found: " << outPath << "\n";
        }

    } while (false);

    // 统一清理 ./temp（无论上面成功或失败）
    std::error_code delEc;
    fs::remove_all("./temp", delEc);
    if (delEc) {
        std::cerr << "Warning: failed to remove ./temp: " << delEc.message() << "\n";
    }

    return rc;
}

GlobalArgs global_args;
int main(int argc, char** argv) {
    // Create a Timer object to record the execution time.
    Timer timer;
    // Create an ArgParser object to parse command line arguments.
    ArgParser parser;
    std::string output = "";
    // Add command line arguments to the ArgParser object.
    parser.add_argument("i", true, "/path/to/input.fasta");
    parser.add_argument_help("i", "The path to the input fasta file.");
    parser.add_argument("o", true, "/path/to/output.fasta");
    parser.add_argument_help("o", "The path to the output fasta file.");
    parser.add_argument("p", false, "mafft");
    parser.add_argument_help("p", "MSA method (mafft, halign3, halign4) or Path to MSA command file.");

    parser.add_argument("t", false, "max_cpu_num");
    parser.add_argument_help("t", "The maximum number of threads that the program runs, the recommended setting is the number of CPUs.");
    parser.add_argument("l", false, "30");
    parser.add_argument_help("l", "The minimum length of MEM, the default value is 30.");
   

    parser.add_argument("f", false, "accurate");
    parser.add_argument_help("f", "The filter MEMs mode. The default is accurate mode.");
    parser.add_argument("v", false, "1");
    parser.add_argument_help("v", "Verbose option, 0 or 1. You could ignore it.");
    parser.add_argument("h", false, "help");
    parser.add_argument_help("h", "print help information");

    // Add command line arguments to the ArgParser object.
    try {
        parser.parse_args(argc, argv);
        global_args.data_path = parser.get("i");
        std::string tmp_thread = parser.get("t");
        if (tmp_thread == "max_cpu_num") {
            global_args.thread = std::thread::hardware_concurrency();
        }
        else {
            global_args.thread = std::stoi(tmp_thread);
        }
        std::string tmp_len = parser.get("l");
        if (tmp_len != "default") {
            global_args.min_mem_length = std::stoi(parser.get("l"));
        }
        else {
            global_args.min_mem_length = -1;
        }

        std::string tmp_filter_mode = parser.get("f");
        if (tmp_filter_mode == "default") {
            global_args.filter_mode = tmp_filter_mode;
        }
        else if (tmp_filter_mode == "fast" || tmp_filter_mode == "accurate") {
            global_args.filter_mode = tmp_filter_mode;
        }
        else {
            throw "filer mode --f parameter should be accurate or fast!";
        }

        global_args.verbose = std::stoi(parser.get("v"));
        if (global_args.verbose != 0 && global_args.verbose != 1) {
            throw "verbose should be 1 or 0";
        }
        std::string cmd_template;
        std::string cmd_path = parser.get("p");
        if (cmd_path == "mafft") {
			cmd_template = "mafft --thread {thread} {input} > {output}";
        }
        else if (cmd_path == "halign3") {
            cmd_template = "halign -t {thread} -o {output} {input}";
		} else if(cmd_path == "halign4"){
            cmd_template = "halign4 {input} {output} -t {thread}";
            #ifdef _WIN32
            std::cout << "halign4 only supports linux" << std::endl;
            exit(1);
            #endif

        }
        else {
            cmd_template = readFile(cmd_path);
        }

        int return_code = test_cmd(cmd_template);
        if (return_code != 0) {
            throw "The command template in " + cmd_path + " is invalid or the MSA software cannot be executed!";
            exit(1);
        }

        global_args.package = cmd_template;
        //if (global_args.package != "halign2" && global_args.package != "halign3" && global_args.package != "mafft") {
        //    throw ("Error: " + global_args.package + " is a invalid method!");
        //}



        global_args.output_path = parser.get("o");
    } // Catch any invalid arguments and print the help message.
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Program Exit!" << std::endl;
        parser.print_help();
        return 1;
    }
    if (global_args.verbose) {
        print_algorithm_info();
    }

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
    if (global_args.verbose) {
        output = "FMAlign2 total time: " + s.str() + " seconds.";
        print_table_line(output);
        print_table_bound();
    }

    return 0;
}