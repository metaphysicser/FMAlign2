# Sequence preservation regression tests

Build with CMake and run the reader/restoration/ID-mapping checks:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_M64=ON
cmake --build build -j 2
ctest --test-dir build --output-on-failure
python tests/test_symbols_cli.py build/symbols_pipeline_test --pipeline
```

On Windows use a MinGW generator and append `.exe` to the executable path.
Python uses its standard library only. Its temporary directory must be writable.

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
biological aligner or an alignment-quality test. The pipeline test does not call
MEM/LCP construction. To exercise that additional path, run:

```sh
python tests/test_symbols_cli.py build/fmalign2
```

The existing uninitialized-LCP defect can crash this last test before fragment
alignment begins. It is not fixed by the symbol-preservation change. Do not treat
the supplied-anchor pipeline test as a successful full-program run or as MAFFT /
HAlign integration validation.
