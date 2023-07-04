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

#ifndef SEQUENCE_SPLIT_ALIGN_H
#define SEQUENCE_SPLIT_ALIGN_H

#include "common.h"
#include "../ext/SW/ssw_cpp.h"
#include "../ext/SW/ssw.h"
#include "utils.h"
#include <algorithm>
#include <sstream>
#ifdef __linux__
#include <sys/stat.h>
#else
#include <direct.h>
#endif
#include <random>
#include <climits>

const std::string TMP_FOLDER = "./temp/";

struct ExpandChainParams {
	std::vector<std::string>* data;
	std::vector<std::vector<std::pair<int_t, int_t>>>* chain;
	uint_t chain_index;
	std::vector<std::vector<std::string>>::iterator result_store;
};

struct ParallelAlignParams {
	std::vector<std::string>* data;
	std::vector<std::vector<std::pair<int_t, int_t>>>::iterator parallel_range;
	uint_t task_index;
	std::vector<std::vector<std::string>>::iterator result_store;
};

/**
* @brief Generates a random string of the specified length.
* This function generates a random string of the specified length. The generated string
* consists of lowercase English letters ('a' to 'z') for Linux platforms, and random bytes
* for Windows platforms.
* @param length The length of the generated string.
* @return A random string of the specified length, or an empty string if an error occurs.
*/
std::string generateRandomString(int length);

/**
* @brief Split and parallel align multiple sequences using a vector of chain pairs.
* This function takes in three parameters: a vector of input sequences (data), a vector of sequence names (name),
* and a vector of chain pairs (chain) that represent initial pairwise alignments between sequences.
* It then splits the chain pairs into smaller regions and performs parallel sequence alignment on these regions.
* Finally, it concatenates the aligned regions and performs sequence-to-profile alignment to generate a final alignment.
* @param data A vector of input sequences to be aligned
* @param name A vector of sequence names corresponding to the input sequences
* @param chain A vector of chain pairs representing initial pairwise alignments between sequences
* @return void
*/
void split_and_parallel_align(std::vector<std::string> data, std::vector<std::string> name, std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence);
/**
* @brief Selects columns from a sequence of split points to enable multi thread.
* @param split_points_on_sequence A vector of vectors of pairs, where each pair represents the start and mem length
* @return A vector of indices of the selected columns.
*/
std::vector<int_t> select_columns(std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence);

/**
* @brief Get a vector of integers that are not in the selected_cols vector and have a maximum value of n.
* @param n The maximum value of the integers in the resulting vector.
* @param selected_cols A vector of integers that are already selected.
* @return A vector of integers that are not in selected_cols and have a maximum value of n.
*/
std::vector<int_t> get_remaining_cols(int_t n, const std::vector<int_t> selected_cols);

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
void* expand_chain(void* arg);

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
	std::vector<std::string>& res_store, uint_t seq_index);

/**
 * @brief Get the range of each sequence in parallel alignment
 * @param data The vector of sequences to be aligned
 * @param chain The vector of chains representing the alignment
 * @return The vector of ranges for each sequence in the alignment
 */
std::vector<std::vector<std::pair<int_t, int_t>>> get_parallel_align_range(std::vector<std::string> data, std::vector<std::vector<std::pair<int_t, int_t>>> chain);

/**
* @brief Function for parallel alignment of sequences.
* This function aligns a subset of input sequences in parallel using multiple threads.
* @param arg Pointer to a ParallelAlignParams struct which contains the input data, the range of sequences to align,
* the index of the current task, and a pointer to the storage for the aligned sequences.
* @return NULL
*/
void* parallel_align(void* arg);

/**
* @brief Align sequences in a FASTA file using either halign or mafft package.
* @param file_name The name of the FASTA file to align.
* @return The name of the resulting aligned FASTA file.
*/
std::string align_fasta(std::string file_name);

/**
* @brief Deletes temporary files generated during sequence alignment tasks.
* @param task_count The number of tasks for which temporary files were created.
*/
void delete_tmp_folder(uint_t task_count);

/**
* @brief Concatenate multiple sequence alignments into a single alignment and write the result to an output file.
* @param concat_string A 2D vector of strings containing the aligned sequences to concatenate.
* @param name A vector of strings containing the names of the sequences.
*/
void concat_alignment(std::vector<std::vector<std::string>>&concat_string, std::vector<std::string> &name);

/**
* @brief Convert sequence fragments into profile by aligning missing fragments with existing ones.
* @param concat_string A reference to a vector of vectors of strings representing concatenated sequence fragments.
* @param data A reference to a vector of strings representing the sequence names.
* @param concat_range A reference to a vector of vectors of pairs of integers representing the start and end positions of the sequence fragments.
* @param fragment_len A reference to a vector of unsigned integers representing the lengths of the sequence fragments.
* @return None.
*/
void seq2profile(std::vector<std::vector<std::string>>& concat_string, std::vector<std::string>& data, std::vector<std::vector<std::pair<int_t, int_t>>>& concat_range, std::vector<uint_t>& fragment_len);

/**
* @brief: Aligns a sequence and a profile using a third-party tool called profile_two_align and returns the iterator pointing to the next position in the 2D vector of strings.
* The function first checks for gaps between the fragments and adds them as necessary.
* Then it writes the sequence to align and a profile file, passes them to profile_two_align, and reads the results.
* Finally, it updates the concatenated string, range and fragment length information accordingly.
* @param seq_index: Index of the sequence to align.
* @param left_index: Index of the left-most fragment.
* @param right_index: Index of the right-most fragment.
* @param concat_string: 2D vector of strings containing the concatenated fragments.
* @param data: Vector of strings containing the sequences.
* @param concat_range: 2D vector of pairs of integers representing the range of each fragment in each sequence.
* @param fragment_len: Vector of unsigned integers representing the length of each fragment.
* @return std::vector<std::vectorstd::string>::iterator: Iterator pointing to the next position in the 2D vector of strings.
*/
std::vector<std::vector<std::string>>::iterator seq2profile_align(uint_t seq_index, uint_t left_index, uint_t right_index, std::vector<std::vector<std::string>>& concat_string, std::vector<std::string>& data, std::vector<std::vector<std::pair<int_t, int_t>>>& concat_range, std::vector<uint_t>& fragment_len);
/**
* @brief Concatenate two sets of sequence data (chain and parallel) into a single set of concatenated data.
* @param chain_string A vector of vectors containing the chain sequence data.
* @param parallel_string A vector of vectors containing the parallel sequence data.
* @return std::vector<std::vectorstd::string> A vector of vectors containing the concatenated sequence data.
*/
std::vector<std::vector<std::string>> concat_chain_and_parallel(std::vector<std::vector<std::string>>& chain_string, std::vector<std::vector<std::string>>& parallel_string);

/**
* @brief Get the length of the first non-empty string in each row of a 2D vector of strings
* @param concat_string The input 2D vector of strings
* @return A vector of unsigned integers representing the length of the first non-empty string in each row
*/
std::vector<uint_t> get_first_nonzero_lengths(const std::vector<std::vector<std::string>>& concat_string);

/**
* @brief Concatenates two sets of ranges in a chain and parallel manner
* Given two vectors of vectors, chain and parallel, this function concatenates the ranges
* in a chain and parallel manner. Specifically, it takes the i-th range from each vector in parallel
* and concatenates them into a single vector. Then, it takes the i-th range from the chain vector
* and concatenates it with the previous vector to form a new concatenated vector. The resulting
* concatenated vector is stored in a new vector of vectors and returned.
* @param chain A vector of vectors representing the chains of ranges to concatenate
* @param parallel A vector of vectors representing the parallel ranges to concatenate
* @return A new vector of vectors representing the concatenated ranges
*/
std::vector<std::vector<std::pair<int_t, int_t>>> concat_chain_and_parallel_range(std::vector<std::vector<std::pair<int_t, int_t>>>& chain, std::vector<std::vector<std::pair<int_t, int_t>>>& parallel);

#endif