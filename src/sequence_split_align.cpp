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

void parallel_align(std::vector<std::string> data, std::vector<std::string> name, std::vector<std::vector<std::pair<int_t, int_t>>> split_points_on_sequence){
    select_columns(split_points_on_sequence);
    return;
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
