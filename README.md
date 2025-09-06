# FMAlign2: A novel fast multiple nucleotide sequence alignment method for ultra-long datasets

[FMAlign2](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btae014/7515251?searchresult=1) is a novel multiple sequence alignment algorithm. It is designed to efficiently align ultra-long nucleotide sequences fast and accurately.

## Installation

### Installation via Conda (recommended)

We recommend installing **FMAlign2** with **conda**.
Since FMAlign2 can be combined with any MSA method and customized parameters, conda makes it convenient to manage dependencies and install different MSA tools.

```bash
conda create -n fmalign2_env
conda activate fmalign2_env
conda install -c malab fmalign2
```

### Installation from Source

Alternatively, you can build FMAlign2 from source with a few simple steps:

```bash
git clone https://github.com/metaphysicser/FMAlign2
cd FMAlign2
make -j
```



## Usage

> Reminder: Please ensure that all external files (such as MAFFT, HALIGN, etc.) are properly copied to their corresponding directories. Pay close attention to the relative paths between FMAlign2 and the ext folder to avoid issues during execution,

if you are Linux user:

   ```shell
   ./FMAlign2 -i /path/to/input.fasta -o /path/to/output.fasta [other options]
   ```

if you are Windows user:

   ```shell
   ./FMAlign2.exe -i /path/to/input.fasta -o /path/to/output.fasta [other options]
   ```
### Parameter detail

 ~~~
 ./FMAlign2 -h
 ~~~

  Parameters Details:

   - -i [file path] **[required]** The path to the input file.
   - -o [output path] [default: ouput.fmaligned2.fasta] The path to the output file.
   - -p [file path] [default: mafft] MSA method (mafft, halign3) or Path to MSA commad file.
   - -t [int] [default: cpu number]  The maximum number of threads that the program runs, the recommended setting is the number of CPUs.
   - -l [int] [default: square root of mean length] The minimum length of MEMs, the default value is square root of mean length.
   - -f [mode] [default: automatic] The filter MEMs mode. If sequence number < 100 → local mode. Otherwise → global modember
   - -v [int] [default:1] Verbose option, 0 or 1. You could ignore it.
   - -h [help] print help information

-----


If you want to evaluate the generated alignment results, you can run the `sp.py` script (requires a Python environment) with the following parameters:

```shell
python sp.py --input output.fmalign2.fasta --match 0 --mismatch 1 --gap1 2 --gap2 2
```

This command will calculate and print the SP (Sum-of-Pairs) score for the multiple sequence alignment results. The `--input` parameter specifies the input alignment file (`output.fmalign2.fasta` in this case), and the `--match`, `--mismatch`, `--gap1`, and `--gap2` parameters define the scoring scheme for matches, mismatches, and gap penalties.

By running this command, you will obtain the SP score, which provides an evaluation of the alignment quality.


## Issue

FMAlign2 is supported by [ZOU's Lab](https://github.com/malabz). If you have any suggestions or feedback, we encourage you to provide them through the issue page on the project's repository. You can also reach out via email to zpl010720@gmail.com.

We value your input and appreciate your contribution to improving the project. Thank you for taking the time to provide feedback, and we will address your concerns as soon as possible.


## Citation

Pinglu Zhang, Huan Liu, Yanming Wei, Yixiao Zhai, Qinzhong Tian, Quan Zou, FMAlign2: a novel fast multiple nucleotide sequence alignment method for ultralong datasets, *Bioinformatics*, 2024;, btae014, https://doi.org/10.1093/bioinformatics/btae014

## License

[Apache 2.0](https://github.com/metaphysicser/FMAlign2/blob/master/LICENSE) © [[MALABZ_UESTC](https://github.com/malabz) [Pinglu Zhang](https://github.com/metaphysicser)]
