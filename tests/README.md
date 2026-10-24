# Regression tests

Build with CMake and run the regression suite:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_M64=ON
cmake --build build -j 2
ctest --test-dir build --output-on-failure
```

On Windows use a MinGW generator and append `.exe` to the executable path.
Python uses its standard library only. Its temporary directory must be writable.
When CMake finds a Python 3 interpreter, CTest also runs both Python-driven tests
below. Otherwise, only the C++ tests are registered. Use
`-DPython3_EXECUTABLE=/path/to/python3` to select an interpreter explicitly.

## LCP initialization

`lcp_initialization_test` exercises the production `find_mem` caller with
test-only `malloc` storage filled with `0x41`. It checks the zero-initialization
contract at the `gsacak` boundary, then runs the unmodified algorithm and verifies
SA permutation/order, every LCP entry by direct prefix comparison, and DA against
the original document boundaries. Generalized suffixes and LCP comparisons stop
at the document separator. Equal suffixes may appear in either order.

The cases include the four-base reproducer (`G` and `CAC`), one-base documents,
identical repeated sequences, and fixed-seed random inputs of 2 x 100, 2 x 1000,
and 20 x 100 bases. A threshold above the longest input also checks that MEM
processing returns the correct number of empty rows. The test fails on the old
`malloc` allocation regardless of whether the system allocator returns zeroed
pages. Instrumentation is confined to the test translation unit; it neither
replaces gSACA-K nor affects the production binary.

Run separate builds with `-DUSE_M64=OFF` and `-DUSE_M64=ON` to check both index
widths. The fix zero-initializes LCP in `find_mem`; no gSACA-K algorithm changes
or small-input bypass are required.

## Sequence preservation

The C++ test checks all eleven ambiguity symbols, lowercase output, U/T restoration,
gap retention, title separation, row order, and rejection of invalid content,
lengths, and IDs. Restoration validates all rows before changing any symbols.

The pipeline test calls the production `split_and_parallel_align` with either no
anchors or a supplied exact anchor. It exercises fragment files, parallel backend
calls, backend row reordering, gap-containing output, concatenation, and restoration
of the original records, including repeated user IDs with distinct descriptions.
Four negative cases check dropped rows, duplicate internal IDs, changed bases, and
truncated rows. All must fail before writing the final output.

`mock_aligner.py` deliberately pads and reorders rows. It is an I/O fixture, not a
biological aligner or an alignment-quality test. The `symbols_pipeline` CTest
uses supplied anchors and does not call MEM/LCP construction. The `symbols_cli`
CTest runs the actual executable, including MEM/LCP construction, with a threshold
above the input lengths (no discovered anchors). Both check two successful cases
and four invalid-backend-output cases. They can also be run manually:

```sh
python tests/test_symbols_cli.py build/symbols_pipeline_test --pipeline
python tests/test_symbols_cli.py build/fmalign2
```

These fixtures validate record, gap, content, and error handling. They do not
validate MAFFT / HAlign integration or biological alignment quality.

## Workspace, input, argument and MEM regression checks

`mem_chains` checks both accurate and fast filtering with empty chains, zero-length
candidates, a single anchor at index zero, disjoint anchors, and overlapping anchors
whose best chain starts at either index zero or a later index. It also checks that
anchors ending exactly at the sequence boundary survive, while overrun anchors do not.

When Python is available, `cli_safety` tests the actual executable with synthetic
FASTA/FASTQ inputs and the mock backend. It verifies that files already present in
`-w` (including nested files and probe-like filenames) survive success, backend
failure, missing output, and malformed command templates. Backend availability
checks use a uniquely created child directory and clean only that child.

Malformed FASTQ quality strings must fail with exit status 1 before alignment output
is written; valid FASTQ must retain all records. Invalid modes, numeric overflow,
nonpositive thread/MEM-length values, trailing numeric junk and missing values must
produce handled errors. Use `-l default` for automatic MEM length. The suite also
checks real MEM discovery and final output with identical sequences in both filter
modes, so losing the sole anchor cannot pass merely because output remains intact.

Custom command templates should quote `"{input}"` and `"{output}"` when paths can
contain spaces. These tests use quoted placeholders and a workspace containing spaces.
