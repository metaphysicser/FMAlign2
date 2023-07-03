# FMAlign2

## Environment

The program is supported both on Linux and Windows. The difference between them is the use of multithreading, on Linux we use pthread, and on Windows we use openmp. So for a better performance experience, we strongly recommend that you use a Linux system to run this project.

## Usage

1. DownLoad

   ```shell
   git clone https://gitee.com/zpl010720/fmalign2
   ```

2. Build

   ```shell
   cd FMAlign2 && make [M64=1]
   ```

   We prepare 32bit  and 64bit mode. Under normal circumstances, 32bit mode is enough to handle most data. But if the length of all sequences spliced exceeds the range of uint32_t (4294967295), you should add the M64 parameter when compiling the program to generate 64bit program.

3. Usage

   ```shell
   ./FMAlign2 --in /path/to/data [options]
   ```

   Parameters Details:

   - -in [file path] [required] The path to the input file.
   - -o [output path] [default: ouput.aligned.fasta] The path to the output file.
   - -t [int] [default: cpu number]  The maximum number of threads that the program runs, the recommended setting is the number of CPUs.
   - -l [int] [default: square root of mean length] The minimum length of MEM, the default value is square root of mean length.
   - -c [float] [default: 0.7] A floating-point parameter that specifies the minimum coverage across all sequences, with values ranging from 0 to 1.
   - -f [mode] [default: fast or accurate] The filter MEM mode. The default setting is that if sequence number less 100, accurate mode otherwise fast mode.
   - -d [int] [default:0] Depth of recursion, you could ignore it.
