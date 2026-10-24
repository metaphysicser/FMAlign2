// Exercise the real caller with nonzero malloc storage, independent of heap state.
// Load all headers before the macros so only mem_finder.cpp is instrumented.
#include "mem_finder.h"
#include <random>
#include <stdexcept>

GlobalArgs global_args;

static unsigned checked_calls = 0;
static unsigned poisoned_allocations = 0;

static void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

static void* poisoned_malloc(size_t size) {
    void* memory = std::malloc(size);
    if (memory) {
        std::memset(memory, 0x41, size);
        ++poisoned_allocations;
    }
    return memory;
}

static int checked_gsacak(unsigned char* text, uint_t* sa, int_t* lcp,
                          int_da* da, uint_t n) {
    require(sa && lcp && da, "Suffix-array allocation failed");
    for (uint_t i = 0; i < n; ++i) {
        require(lcp[i] == 0, "find_mem passed a nonzero initial LCP entry at index " +
                            std::to_string(i));
    }

    // Run the unmodified algorithm, then check its outputs against naive oracles.
    int depth = gsacak(text, sa, lcp, da, n);
    require(depth >= 0, "gsacak failed");
    std::vector<bool> seen(n, false);
    std::vector<int_da> document(n);
    int_da current_document = 0;
    for (uint_t i = 0; i < n; ++i) {
        document[i] = current_document;
        if (text[i] == 1) ++current_document;
    }
    std::string previous_suffix;
    for (uint_t i = 0; i < n; ++i) {
        require(sa[i] < n, "SA entry out of range");
        require(!seen[sa[i]], "SA is not a permutation");
        seen[sa[i]] = true;
        uint_t end = sa[i];
        while (end + 1 < n && text[end] > 1) ++end;
        // Generalized suffixes end at their document separator (or sentinel).
        std::string suffix(reinterpret_cast<const char*>(text + sa[i]),
                           end - sa[i] + 1);
        require(i == 0 || previous_suffix <= suffix, "SA is not sorted");
        previous_suffix = suffix;
        uint_t expected_lcp = 0;
        if (i > 0) {
            uint_t left = sa[i - 1], right = sa[i];
            while (left + expected_lcp < n && right + expected_lcp < n &&
                   text[left + expected_lcp] > 1 &&
                   text[left + expected_lcp] == text[right + expected_lcp]) {
                ++expected_lcp;
            }
        }
        require(lcp[i] == static_cast<int_t>(expected_lcp),
                "LCP differs from direct prefix comparison at index " + std::to_string(i));
        require(da[i] == document[sa[i]], "DA differs from the source document");
    }
    ++checked_calls;
    return depth;
}

// The shim affects this translation unit only. calloc and gsacak's own
// allocations stay unchanged; no test hooks are added to the production binary.
#define malloc poisoned_malloc
#define gsacak checked_gsacak
#include "../src/mem_finder.cpp"
#undef gsacak
#undef malloc

static void check_case(const char* name, const std::vector<std::string>& sequences) {
    global_args = GlobalArgs{};
    global_args.thread = 2;
    global_args.filter_mode = "accurate";
    for (const auto& sequence : sequences) {
        global_args.min_mem_length = std::max(global_args.min_mem_length,
                                              static_cast<int_t>(sequence.size() + 1));
    }
    unsigned calls_before = checked_calls;
    unsigned allocations_before = poisoned_allocations;
    auto chains = find_mem(sequences);
    require(checked_calls == calls_before + 1, "find_mem did not call gsacak");
    require(poisoned_allocations >= allocations_before + 2, "Allocation shim was bypassed");
    require(chains.size() == sequences.size(), "Unexpected MEM result row count");
    for (const auto& chain : chains) {
        require(chain.empty(), "MEM exceeds the longest input sequence");
    }
    std::cout << name << ": SA/LCP/DA and caller checks passed\n";
}

int main() {
    try {
        // G/CAC is the four-base reproducer that used garbage as a prefix length.
        check_case("minimal_G_CAC", {"G", "CAC"});
        check_case("single_base_documents", {"A", "C"});
        std::string repeated;
        for (int i = 0; i < 25; ++i) repeated += "ACGT";
        check_case("identical_repeated_documents", {repeated, repeated});

        std::mt19937 random(20261024);
        for (const auto& shape : {std::make_pair(2, 100), std::make_pair(2, 1000),
                                  std::make_pair(20, 100)}) {
            std::vector<std::string> sequences(shape.first, std::string(shape.second, 'A'));
            for (auto& sequence : sequences) {
                for (auto& base : sequence) base = "ACGT"[random() % 4];
            }
            std::string name = "random_" + std::to_string(shape.first) + "x" +
                               std::to_string(shape.second);
            check_case(name.c_str(), sequences);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
