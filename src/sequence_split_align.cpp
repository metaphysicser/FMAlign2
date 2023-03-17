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

#include "../include/sequence_split_align.h"

void split_and_parallel_align(std::vector<std::string> data, std::vector<std::string> name, std::vector<std::vector<std::pair<int_t, int_t>>> chain){
    std::cout << "#                Parallel Aligning...                       #" << std::endl;
    print_table_divider();
    std::string output = "";
    Timer timer;
    uint_t chain_num = chain[0].size();
    uint_t seq_num = data.size();
    std::vector<std::vector<std::string>> chain_string(chain_num); // chain_num * seq_num

    std::vector<ExpandChainParams> params(chain_num);

    for (uint_t i = 0; i < chain_num; i++) {
        params[i].data = &data;
        params[i].chain = &chain;
        params[i].chain_index = i;
        params[i].result_store = chain_string.begin() + i;
    }

    for (uint_t i = 0; i < chain_num; i++) {
        expand_chain(&params[i]);
    }
    params.clear();

#if DEBUG   
    for (uint_t j = 0; j < chain_num; ++j) {
        uint_t len = chain_string[j][0].length();
        for (uint_t i = 1; i < seq_num; ++i) {
            if (chain_string[j][i].length() != len) {
                std::cout <<"chain_index: " << j << " seq_index: " << i << " " << chain_string[j][i].length() << "!=" << len << std::endl;
            }
        }
    }
#endif    

    double SW_time = timer.elapsed_time();
    std::stringstream s;
    s << std::fixed << std::setprecision(2) << SW_time;
    output = "SW expand time: " + s.str() + " seconds.";
    print_table_line(output);
    timer.reset();

    if (0 != access(TMP_FOLDER, 0))
    {
        if (0 != mkdir(TMP_FOLDER)) {
            std::cerr << "Fail to create file folder " << TMP_FOLDER << std::endl;
        }
    }


    std::vector<std::vector<std::pair<int_t, int_t>>> parallel_align_range = get_parallel_align_range(data, chain);
    uint_t parallel_num = parallel_align_range.size();
    std::vector<std::vector<std::string>> parallel_string(parallel_num, std::vector<std::string>(seq_num));

    std::vector<ParallelAlignParams> parallel_params(parallel_num);
#if (defined(__linux__))
    threadpool pool;
    threadpool_init(&pool, global_args.thread);
    for (uint_t i = 0; i < parallel_num; i++) {
        parallel_params[i].data = &data;
        parallel_params[i].parallel_range = parallel_align_range.begin()+i;
        parallel_params[i].task_index = i;
        parallel_params[i].result_store = parallel_string.begin() + i;
        threadpool_add_task(&pool, parallel_align, &parallel_params[i]);
    }
    threadpool_destroy(&pool);
#else
#pragma omp parallel for num_threads(global_args.thread)
    for (uint_t i = 0; i < parallel_num; i++) {
        parallel_params[i].data = &data;
        parallel_params[i].parallel_range = parallel_align_range.begin() + i;
        parallel_params[i].task_index = i;
        parallel_params[i].result_store = parallel_string.begin() + i;
        parallel_align(&parallel_params[i]);
    }
#endif
    delete_tmp_folder(parallel_num);

    print_table_divider();
    return;
}

/**
* @brief Get a vector of integers that are not in the selected_cols vector and have a maximum value of n.
* @param n The maximum value of the integers in the resulting vector.
* @param selected_cols A vector of integers that are already selected.
* @return A vector of integers that are not in selected_cols and have a maximum value of n.
*/
std::vector<int_t> get_remaining_cols(int_t n, const std::vector<int_t> selected_cols) {
    // Initialize a boolean vector to indicate whether an integer is selected.
    std::vector<bool> is_selected(n, false);
    // Mark the integers in the selected_cols vector as selected.
    for (int_t i : selected_cols) {
        if (i < n) {
            is_selected[i] = true;
        }
    }
    // Get the integers that are not selected and store them in a new vector.
    std::vector<int_t> remaining;
    for (int i = 0; i < n; i++) {
        if (!is_selected[i]) {
            remaining.push_back(i);
        }
    }
    return remaining;
}
/**
* @brief Selects columns from a sequence of split points to enable multi thread.
* @param split_points_on_sequence A vector of vectors of pairs, where each pair represents the start and mem length
* @return A vector of indices of the selected columns.
*/
std::vector<int_t> select_columns(std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence) {
    // Get the number of columns and rows in the split points sequence.
    uint_t col_num = split_points_on_sequence[0].size();
    uint_t row_num = split_points_on_sequence.size();
    // Create a vector to keep track of whether each column needs to be changed.
    std::vector<bool> col_need_change(col_num, false);
    // Create a vector to store the indices of the selected columns.
    std::vector<int_t> selected_cols;
    int_t count = col_num;
    // Create a vector to store the number of effect columns for each column.
    std::vector<std::pair<uint_t, uint_t>> effect_col_num(col_num);

    for (uint_t i = 0; i < col_num; i++) {
        effect_col_num[i].first = 0;
        effect_col_num[i].second = i;
        for (uint_t j = 0; j < row_num; j++) {
            if (split_points_on_sequence[j][i].first == -1) {
                col_need_change[i] = true;
                count--;
                break;
            }
        }
    }

    // Compute the number of effect columns for each column.
    for (uint_t i = 0; i < col_num; i++) {
        if (col_need_change[i]) {
            continue;
        }
        bool has_left = (i == 0 || col_need_change[i - 1] == false);
        bool has_right = (i == col_num - 1 || col_need_change[i + 1] == false);

        if (has_left && has_right) {
            col_need_change[i] = true;
            count--;
            continue;
        }
        effect_col_num[i].first += 1;

        if (has_left && i < col_num-1) {
            effect_col_num[i + 1].first += 1;
        }

        if (has_right && i > 0) {
            effect_col_num[i - 1].first += 1;
        }
    }
    // Sort the columns based on the number of effect columns.
    std::sort(effect_col_num.begin(), effect_col_num.end(), [](const std::pair<uint_t, uint_t>& a, const std::pair<uint_t, uint_t>& b) { return a.first > b.first; });
    
    for (uint_t i = 0; i < col_num; i++) {
        if (count <= 0) {
            break;
        }
        uint_t effect_num = effect_col_num[i].first;
        uint_t col_index = effect_col_num[i].second;
        if (effect_num <= 0) {
            std::cerr << "some bugs occur in select column." << std::endl;
            exit(-1);
        }
        if (col_need_change[col_index] == false) {
            col_need_change[col_index] = true;
            selected_cols.push_back(col_index);
            count--;
        }
        if (effect_num > 1) {
            if ((col_index == 1 || (col_index >= 2 && col_need_change[col_index - 2]) == true) && (col_index >= 1 && col_need_change[col_index - 1] == false)) {
            col_need_change[col_index - 1] = true;
            count--;
            }

            if ((col_index == col_num - 2 || (col_index + 2 < col_num && col_need_change[col_index + 2]) == true) && (col_index + 1 < col_num && col_need_change[col_index + 1] == false)) {
                col_need_change[col_index + 1] = true;
                count--;
            }
        }
        
    }
    return selected_cols;
}

/**
@brief Expands the chain at the given index for all sequences in the input data.
This function takes a void pointer to input arguments and casts it to the correct struct type.
It then retrieves the required variables, which include the data and chain input parameters, and the chain index.
The query sequence is obtained from the chain for each sequence in the data input, and a reference sequence is obtained
using the neighboring chains. If there is no neighboring chain, the reference sequence is obtained from the start of
the sequence to the end of the previous neighboring chain or to the end of the sequence if there is no previous neighboring chain.
The aligner is then used to align the query sequence and the reference sequence, and the result is stored in an alignment object.
If the chain at the given index for the current sequence is empty, the result is stored in the chain. Otherwise, the query sequence
is already aligned, and the aligned fragment is stored in the aligned_fragment vector.
Finally, the function stores the aligned fragments in the result_store vector.
@param arg A void pointer to input arguments.
@return NULL
*/
void* expand_chain(void* arg) {
    // Cast the input parameters to the correct struct type
    ExpandChainParams* ptr = static_cast<ExpandChainParams*>(arg);
    // Get data, chain, and chain_index from the input parameters
    const std::vector<std::string> data = *(ptr->data);
    std::vector<std::vector<std::pair<int_t, int_t>>> chain = *(ptr->chain);
    const uint_t chain_index = ptr->chain_index;
    // std::cout << "in" << chain_index << '\n';
    // Get the number of sequences in the data vector and the number of chains in the current chain
    uint_t seq_num = data.size();
    uint_t chain_num = chain[0].size();
    
    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;

    uint_t query_length = 0;
    std::string query = "";
    std::vector<std::string> aligned_fragment(seq_num);
    // Find the query sequence and its length in the current chain
    for (uint_t i = 0; i < seq_num; i++) {
        if (chain[i][chain_index].first != -1) {
            query_length = chain[i][chain_index].second;
            query = data[i].substr(chain[i][chain_index].first, query_length);
            break;
        }
    }
    

    for (uint_t i = 0; i < seq_num; i++) {
        int_t begin_pos = chain[i][chain_index].first;
        // If the begin position is -1, the current subsequence is unaligned
        if (begin_pos == -1) {
            uint_t tmp_index = chain_index;
            // Find the beginning and end positions of the unaligned subsequence
            uint_t ref_begin_pos = 0;
            uint_t ref_end_pos = 0;
            int_t maskLen = query_length / 2;
            maskLen = maskLen < 15 ? 15 : maskLen;

            for (; tmp_index > 0 && chain[i][tmp_index - 1].first == -1; --tmp_index);

            ref_begin_pos = tmp_index <= 0 ? 0 : chain[i][tmp_index-1].first + chain[i][tmp_index - 1].second;

            tmp_index = chain_index;
            for (; tmp_index < chain_num - 1 && chain[i][tmp_index + 1].first == -1; ++tmp_index);

            ref_end_pos = tmp_index >= chain_num - 1 ? data[i].length() - 1 : chain[i][tmp_index + 1].first;

            std::string ref = data[i].substr(ref_begin_pos, ref_end_pos - ref_begin_pos);

            // Get the reference subsequence and align it with the query subsequence
            aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);

            std::pair<int_t, int_t> p = store_sw_alignment(alignment, ref, query, aligned_fragment, i);
       
            if (p.first != -1) {
                p.first += ref_begin_pos;
                (*(ptr->chain))[i][chain_index] = p;
            }
           
        }
        else {
            aligned_fragment[i] = query;
        }
    }
    *(ptr->result_store) = aligned_fragment;
#if DEBUG
    std::cout << chain_index << '\n';
#endif    
    return NULL;
}

/**
* @brief Store the Smith-Waterman alignment results in a vector of aligned sequences.
* This function takes the alignment results generated by the StripedSmithWaterman algorithm and
* stores the aligned sequences in a vector of strings. The function also returns the alignment
* start and length as a pair of integers. If the alignment failed, the function returns (-1,-1).
* @param alignment The alignment results generated by StripedSmithWaterman algorithm.
* @param ref The reference sequence.
* @param query The query sequence.
* @param res_store A vector of aligned sequences.
* @param seq_index The index of the query sequence in the vector of aligned sequences.
* @return A pair of integers representing the alignment start and length.
*/
std::pair<int_t, int_t> store_sw_alignment(StripedSmithWaterman::Alignment alignment, std::string& ref, std::string& query,
    std::vector<std::string>& res_store, uint_t seq_index)
{
    // Extract cigar string from the alignment
    std::vector<unsigned int> cigar = alignment.cigar;
    // Extract the start and end positions of the alignment on the reference sequence
    int_t ref_begin = alignment.ref_begin;
    int_t ref_end = alignment.ref_end;

    uint_t query_begin = 0;

    std::string aligned_result = "";
    // If the alignment failed, return (-1,-1)
    if (ref_begin <= -1) {
        res_store[seq_index] = aligned_result;
        return std::make_pair(-1, -1);
    }

    int_t S_count = 0;
    int_t total_length = alignment.query_end;
    if (cigar_int_to_op(cigar[0]) == 'S') {
        S_count += cigar_int_to_len(cigar[0]);
    }
    if (cigar_int_to_op(cigar[cigar.size() - 1]) == 'S') {
        S_count += cigar_int_to_len(cigar[cigar.size() - 1]);
        total_length += cigar_int_to_len(cigar[cigar.size() - 1]);
    }

    if (S_count > ceil(0.8 * total_length)) {
            res_store[seq_index] = aligned_result;
            return std::make_pair(-1, -1);
    }
    uint_t p_ref = 0;
    uint_t p_query = 0;

    for (uint_t i = 0; i < cigar.size(); i++) {
        char op = cigar_int_to_op(cigar[i]);
        int_t len = cigar_int_to_len(cigar[i]);

        switch (op)
        {
        case 'S': {
            // Handle soft clipping at the beginning and end of the alignment
            if (i == 0) {
                int_t tmp_len = len;
                if (ref_begin <= len) {
                    while (tmp_len > ref_begin) {
                        aligned_result += "-";
                        tmp_len--;
                    }                 
                }
                for (int_t j = ref_begin - tmp_len; j < ref_begin; j++) {
                    aligned_result += ref[j];
                }        
            }
            else {
                int_t tmp_len = len;
                for (uint_t j = ref_end; j < ref.length() && tmp_len > 0; j++) {
                    aligned_result += ref[j];
                    tmp_len--;
                }
                while (tmp_len > 0) {
                    aligned_result += "-";
                    tmp_len--;
                }
            }
            p_query += len;
            break;
        }
        case 'M':
        case 'X':
        case '=': {
            // Handle match, mismatch, and substitution operations
            for (int_t j = 0; j < len; j++) {
                aligned_result += ref[ref_begin + p_ref + j];
            }
            p_ref += len;
            p_query += len;
            break;
        }
        case 'I': {
            // Handle insertion operations
            for (int_t j = 0; j < len; j++) {
                aligned_result += '-';
            }
            p_query += len;
            break;
        }
        case 'D': {
            // Handle deletion operations
            std::string gaps = "";
            for (int_t j = 0; j < len; j++) {
                gaps += '-';
            }
            for (int_t j = 0; j < len; j++) {
                aligned_result += ref[ref_begin + p_ref + j];
            }
            p_ref += len;

            query.insert(query_begin + p_query, gaps);
            for (uint_t j = 0; j < seq_index; j++) {
                if (res_store[j].length() != 0) {
                    res_store[j].insert(query_begin + p_query, gaps);
                }
                
            }
            p_query += len;
            break;
        }

        default:
            break;
        }
    }
    res_store[seq_index] = aligned_result;

    std::pair<int, int> p(ref_begin, ref_end - ref_begin + 1);
    return p;
}

/**
 * @brief Get the range of each sequence in parallel alignment
 * @param data The vector of sequences to be aligned
 * @param chain The vector of chains representing the alignment
 * @return The vector of ranges for each sequence in the alignment
 */
std::vector<std::vector<std::pair<int_t, int_t>>> get_parallel_align_range(std::vector<std::string> data, std::vector<std::vector<std::pair<int_t, int_t>>> chain) {
    // Get the number of sequences and chains
    uint_t seq_num = data.size();
    uint_t chain_num = chain[0].size();
    // Initialize the vector to store the ranges
    std::vector<std::vector<std::pair<int_t, int_t>>> parallel_align_range(seq_num);

    // Iterate through each sequence
    for (uint_t i = 0; i < seq_num; i++) {
        // Initialize the last position as 0 and a temporary vector to store the ranges
        int_t last_pos = 0;
        std::vector<std::pair<int_t, int_t>> tmp_range;
        // Iterate through each chain for the current sequence
        for (uint_t j = 0; j < chain_num; j++) {
            // Get the begin position for the current chain
            int_t begin_pos = chain[i][j].first;
            // If the chain cannot be aligned, add (-1,-1) to the range and set the last position as -1
            if (begin_pos == -1) {
                tmp_range.push_back(std::make_pair(-1, -1));
                last_pos = -1;
            }
            else {
                // If the last chain could not be aligned, add (-1,-1) to the range
                if (last_pos == -1) {
                    tmp_range.push_back(std::make_pair(-1, -1));
                }
                else {
                    // Add the range between the last position and the current chain's begin position
                    tmp_range.push_back(std::make_pair(last_pos, begin_pos - last_pos));
                }
                // Update the last position to the end of the current chain
                last_pos = begin_pos + chain[i][j].second;
            }
        }
        // If the last chain could not be aligned, add (-1,-1) to the range
        if (last_pos == -1) {
            tmp_range.push_back(std::make_pair(-1, -1));
        }
        else {
            // Add the range between the last position and the end of the sequence
            tmp_range.push_back(std::make_pair(last_pos, data[i].length() - last_pos));
        }
        // Add the ranges for the current sequence to the vector of ranges
        parallel_align_range[i] = tmp_range;
    }
    std::vector<std::vector<std::pair<int_t, int_t>>> transpose_res(parallel_align_range[0].size(), std::vector<std::pair<int_t, int_t>>(seq_num));
    for (uint_t i = 0; i < seq_num; i++) {
        for (uint_t j = 0; j < parallel_align_range[0].size(); j++) {
            transpose_res[j][i] = parallel_align_range[i][j];
        }
    }
    parallel_align_range.clear();
    // Return the vector of ranges
    return transpose_res;
}

void* parallel_align(void* arg) {
    // Cast the input parameters to the correct struct type
    ParallelAlignParams* ptr = static_cast<ParallelAlignParams*>(arg);
    // Get data, chain, and chain_index from the input parameters
    const std::vector<std::string> data = *(ptr->data);
    std::vector<std::pair<int_t, int_t>> parallel_range = *(ptr->parallel_range);
    const uint_t task_index = ptr->task_index;
    // Get the number of sequences in the data vector and the number of chains in the current chain
    uint_t seq_num = data.size();
    std::string file_name = "./tmp/task-" + std::to_string(task_index) + ".fasta";
    std::ofstream file;
    file.open(file_name);

    if (!file.is_open()) {
        std::cerr << file_name << " fail to open!" << std::endl;
        exit(1);
    }
    for (uint_t i = 0; i < seq_num; i++) {  
        if (parallel_range[i].first >= 0) {
            std::string seq_content = data[i].substr(parallel_range[i].first, parallel_range[i].second);
            std::stringstream sstreasm;
            sstreasm << ">SEQENCE" << i << "\n" << seq_content << "\n";
            file << sstreasm.str();
        }       
    }
    file.close();
    std::string res_file_name = align_fasta(file_name);
    

   
    return NULL;
}

std::string align_fasta(std::string file_name) {
    std::string cmnd = "";
    std::string res_file_name = file_name.substr(0, file_name.find(".fasta")) + ".aligned.fasta";
    if (global_args.package == "halign") {
         cmnd.append("java -jar ./ext/halign3/share/halign-stmsa.jar ")
         .append("-t 1").append(" -o ").append(res_file_name).append(" ").append(file_name);
    }
    else if (global_args.package == "mafft") {

    }

#if (defined(__linux__))
    cmnd.append(" > /dev/null");
#else
    cmnd.append(" > NUL");
#endif

    system(cmnd.c_str());
    return res_file_name;
}

/**
* @brief Deletes temporary files generated during sequence alignment tasks.
* @param task_count The number of tasks for which temporary files were created.
*/
void delete_tmp_folder(uint_t task_count) {
    for (uint_t i = 0; i < task_count; i++) {
        std::string file_name = "./tmp/task-" + std::to_string(i) + ".fasta";
        std::string res_file_name = "./tmp/task-" + std::to_string(i) + ".aligned.fasta";
        if (remove(file_name.c_str()) != 0) {
            std::cerr << "Error deleting file " << file_name << std::endl;
        }
        if (remove(res_file_name.c_str()) != 0) {
            std::cerr << "Error deleting file " << res_file_name << std::endl;
        }
    }
}