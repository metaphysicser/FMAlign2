// Exercise the production split/align/output path with a supplied MEM result.
// Keeps fragment/anchor handling independent of MEM discovery; see the LCP test.
#include "sequence_split_align.h"
#include <fstream>
#include <iterator>

GlobalArgs global_args;

int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::runtime_error("Expected input, output, command file, and anchor mode");
        std::ifstream command(argv[3]);
        global_args.package = std::string(std::istreambuf_iterator<char>(command), {});
        global_args.output_path = argv[2];
        global_args.temp_folder = "work";
        global_args.thread = 2;
        global_args.avg_file_size = 1;
        global_args.min_seq_coverage = 1;
        std::vector<std::string> original, names;
        read_data(argv[1], original, names, false, false);
        std::vector<std::string> normalized;
        for (const auto& sequence : original) normalized.push_back(clean_sequence(sequence));
        std::vector<std::vector<std::pair<int_t, int_t>>> chains(original.size());
        if (std::string(argv[4]) == "anchor") {
            for (auto& chain : chains) chain.push_back({120, 30});
        }
        split_and_parallel_align(normalized, names, chains, original);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
