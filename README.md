# FMAlign2: A novel fast multiple nucleotide sequence alignment method for ultra-long datasets

[FMAlign2](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btae014/7515251?searchresult=1) is a novel multiple sequence alignment algorithm. It is designed to efficiently align ultra-long nucleotide sequences fast and accurately.


## Installation via Conda (recommended)

We recommend installing **FMAlign2** with **conda**.
By default, the Conda environment installs three MSA methods: [MAFFT](https://mafft.cbrc.jp/alignment/software/), HAlign3, and [HAlign4](https://github.com/metaphysicser/HAlign-4). Since FMAlign2 can integrate with any MSA method and supports custom parameters, Conda provides a convenient way to manage dependencies and install alternative MSA tools.

```bash
conda create -n fmalign2_env
conda activate fmalign2_env
conda install -c malab fmalign2
```

## Usage

if you are Linux user:

   ```shell
   FMAlign2 -i /path/to/input.fasta -o /path/to/output.fasta [other options]
   ```

if you are Windows user:

   ```shell
   FMAlign2.exe -i /path/to/input.fasta -o /path/to/output.fasta [other options]
   ```

### Parameter Details

```
./FMAlign2 -h
```

**Parameters**

* `-i <file>` **Required.** Path to the input FASTA file.
* `-o <file>` **Required.** Path to the output FASTA file.
* `-p <method|file>` (default: `mafft`). MSA backend — `mafft`, `halign3`, or `halign4` — or a path to a custom MSA command file.
* `-t <int>` (default: number of available CPU cores). Maximum number of threads to use.
* `-l <int>` (default: 30). Minimum MEM length.
* `-f <mode>` (default: `accurate`). MEM filtering mode; use `fast` to speed up at the cost of sensitivity.
* `-v <0|1>` (default: 1). Verbosity flag.
* `-h` Show help information and exit.

**Notes**

* If too few MEMs are found, try lowering `-l` (e.g., below 30).
* If MEM detection takes too long, consider `-f fast`; it runs faster but may find fewer MEMs.

If you want to evaluate the generated alignment results, you can run the `sp.py` script (requires a Python environment) with the following parameters:

```shell
python sp.py --input output.fmalign2.fasta --match 0 --mismatch 1 --gap1 2 --gap2 2
```

This command will calculate and print the SP (Sum-of-Pairs) score for the multiple sequence alignment results. The `--input` parameter specifies the input alignment file (`output.fmalign2.fasta` in this case), and the `--match`, `--mismatch`, `--gap1`, and `--gap2` parameters define the scoring scheme for matches, mismatches, and gap penalties.

By running this command, you will obtain the SP score, which provides an evaluation of the alignment quality.

### Installation from Source

Alternatively, you can build FMAlign2 from source with a few simple steps:

```bash
git clone https://github.com/metaphysicser/FMAlign2
cd FMAlign2
make -j
```


## Issue

FMAlign2 is supported by [ZOU's Lab](https://github.com/malabz). If you have any suggestions or feedback, we encourage you to provide them through the issue page on the project's repository. You can also reach out via email to zpl010720@gmail.com.

We value your input and appreciate your contribution to improving the project. Thank you for taking the time to provide feedback, and we will address your concerns as soon as possible.


## Citation

Pinglu Zhang, Huan Liu, Yanming Wei, Yixiao Zhai, Qinzhong Tian, Quan Zou, FMAlign2: a novel fast multiple nucleotide sequence alignment method for ultralong datasets, *Bioinformatics*, 2024;, btae014, https://doi.org/10.1093/bioinformatics/btae014

## License

[Apache 2.0](https://github.com/metaphysicser/FMAlign2/blob/master/LICENSE) © [[MALABZ_UESTC](https://github.com/malabz) [Pinglu Zhang](https://github.com/metaphysicser)]
