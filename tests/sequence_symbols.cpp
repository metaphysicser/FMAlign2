#include "utils.h"
#include <functional>

GlobalArgs global_args;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static void rejects(const std::function<void()>& action) {
    try { action(); } catch (const std::runtime_error&) { return; }
    throw std::runtime_error("Invalid alignment was accepted");
}

int main(int argc, char** argv) {
    try {
        require(argc == 2, "Expected the reader fixture path");
        global_args.verbose = 0;
        std::vector<std::string> data, names;
        read_data(argv[1], data, names, false, true);
        require(data == std::vector<std::string>({"A-RYSWKMBDHVNU-", "ACGT----------"}), "Reader lost symbols or gaps");
        require(names[0] == "s1 descriptive title", "Reader lost title separator");

        const std::vector<std::string> original = {"RACGY", "NTUHV"};
        std::vector<std::string> aligned = {"-nacg-n", "n-ttnn-"};
        restore_sequence_symbols(original, aligned);
        require(aligned == std::vector<std::string>({"-RACG-Y", "N-TUHV-"}), "Symbol restoration failed");
        for (char symbol : std::string("RYSWKMBDHVN")) {
            std::vector<std::string> source = {std::string(1, symbol) + "AC" + symbol + symbol};
            std::vector<std::string> row = {"-N-A-CNN-"};
            restore_sequence_symbols(source, row);
            require(row[0] == std::string("-") + symbol + "-A-C" + symbol + symbol + "-", "IUPAC restoration failed");
        }
        for (const std::vector<std::string>& bad : std::vector<std::vector<std::string>>{
                {"NACG"}, {"NACGNN"}, {"NATGN"}, {"NACGX"}, {"-----"}}) {
            auto rows = bad;
            rejects([&] { restore_sequence_symbols({"RACGY"}, rows); });
        }
        std::vector<std::string> bad_shape = {"AC", "AC-"};
        rejects([&] { restore_sequence_symbols({"AC", "AC"}, bad_shape); });
        std::vector<std::string> bad_count = {"AC"};
        rejects([&] { restore_sequence_symbols({"AC", "AC"}, bad_count); });
        std::vector<std::string> reordered = {"TN-", "-AC"};
        order_alignment({"SEQENCE0", "SEQENCE2"}, {"SEQENCE2 description", "SEQENCE0"}, reordered);
        require(reordered == std::vector<std::string>({"-AC", "TN-"}), "Row order not restored");
        rejects([&] { order_alignment({"a", "b"}, {"a", "a"}, reordered); });
        rejects([&] { order_alignment({"a", "b"}, {"a", "c"}, reordered); });
        rejects([&] { order_alignment({"a"}, {"a", "b"}, reordered); });
        std::cout << "Sequence symbol and alignment integrity checks passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
