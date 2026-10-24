#include "mem_finder.h"
#include <stdexcept>

GlobalArgs global_args;

static void check(const std::vector<std::pair<int_t, int_t>>& candidates,
                  const std::vector<std::pair<int_t, int_t>>& expected) {
    for (bool fast : {false, true}) {
        std::vector<mem> mems;
        for (auto candidate : candidates) {
            mem anchor{};
            anchor.mem_length = candidate.second;
            anchor.avg_pos = candidate.first;
            anchor.substrings = {{0, static_cast<uint_t>(candidate.first), nullptr},
                                 {1, static_cast<uint_t>(candidate.first), nullptr}};
            mems.push_back(anchor);
        }
        auto result = fast ? filter_mem_fast(mems, 2) : filter_mem_accurate(mems, 2);
        if (result.size() != 2 || result[0] != expected || result[1] != expected) {
            throw std::runtime_error(std::string(fast ? "fast" : "accurate") +
                                     " MEM backtracking differs from expected chain");
        }
    }
}

int main() {
    try {
        global_args.thread = 2;
        global_args.min_seq_coverage = 1;
        check({}, {});
        check({{10, 0}}, {});
        check({{10, 30}}, {{10, 30}});
        check({{10, 30}, {60, 20}}, {{10, 30}, {60, 20}});
        check({{10, 30}, {20, 10}}, {{10, 30}});
        check({{10, 10}, {15, 30}}, {{15, 30}});
        // An exclusive end equal to the sequence length is still in bounds.
        uint_t index = 0;
        for (uint_t position : {2u, 3u}) {
            mem anchor{};
            anchor.mem_length = 4;
            anchor.mem_index = &index;
            anchor.substrings = {{0, position, &index}, {1, position, &index}};
            std::vector<mem> mems{anchor};
            sort_mem(mems, {"ACGTAC", "ACGTAC"});
            if (mems.size() != (position == 2 ? 1u : 0u)) {
                throw std::runtime_error("Incorrect MEM end-boundary filtering");
            }
        }
        std::cout << "Empty, single, index-zero, disjoint and overlapping chains passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
