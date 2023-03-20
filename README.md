# FMAlign2

## Usage

1. DownLoad

   ```shell
   git clone https://gitee.com/zpl010720/fmalign2
   ```

2. Build

   ```
   cd FMAlign2 && make [M64=1]
   ```

3. Usage

   ```shell
   ./FMAlign2 --in /path/to/data [options]
   ```

   Parameters Details:

   - --in [file path] [required] The path to the input file.
   - --o [output path] [default: ouput.aligned.fasta] The path to the output file.
   - --t [int] [default: cpu number]  The maximum number of threads that the program runs, the recommended setting is the number of CPUs
   - --l [int] [default: 30] The minimum length of MEM
   - --c [float] [default: 0.7] A floating-point parameter that specifies the minimum coverage across all sequences, with values ranging from 0 to 1.
