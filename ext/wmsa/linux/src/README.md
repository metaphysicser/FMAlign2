# WMSA

A method of Multiple Sequence Alignment, which the writer is Wym6912.

![version](https://anaconda.org/wym6912/wmsa/badges/version.svg) ![last updated](https://anaconda.org/wym6912/wmsa/badges/latest_release_date.svg) ![supported platforms](https://anaconda.org/wym6912/wmsa/badges/platforms.svg) ![conda](https://anaconda.org/wym6912/wmsa/badges/installer/conda.svg)

## Recommended environment

Liunx-based systems, like `Ubuntu`, `CentOS` .

If you are a `Windows 10` user, you can use it by `WSL` ([Manually download Windows Subsystem for Linux (WSL) Distros | Microsoft Docs](https://docs.microsoft.com/en-us/windows/wsl/install-manual)).

## How to use this program

### Directly clone and complie the project (Recommended for macOS users)

```bash
git clone github.com/wym6912/WMSA --recursive 
# add --recursive arugment to get cd-hit and mafft subprograms
make all THREADS=16 # not -j16 because the THREADS=16 can be used by the subprograms
make install # program is installed on /usr/bin by default
wmsa -H # help message on this program
```

You also can download source code on the [release](https://github.com/malabz/WMSA/releases).

If you are an MACOS user, please use the following `make` command:

```bash
make all THREADS=16 NOOPENMP=yes
```

in order to avoid the bug on `openmp` library.

### Download from `conda` (Recommended for Linux users)

We recommended use `conda` to configuration for Linux users. See [Install Miniconda](https://docs.conda.io/en/latest/miniconda.html) if not installed `conda` or `miniconda` before. The command of configuration like this:

```bash
conda create -n wmsa # make an new environment for running wmsa
conda activate wmsa
conda install wmsa -c wym6912
wmsa -H
```

## How to run the method on a test set

The test data is published [here](https://github.com/malabz/WMSA-dataset).

You can clone and use the following project by using this command:

```bash
git clone https://github.com/malabz/WMSA-dataset
# test mt1x dataset on WMSA
cd WMSA-dataset/mt/
unzip mt1x.zip
/usr/bin/wmsa -i mt1x.fasta -o mt1x.wmsa.fasta -T 16 -c 0.9
```

The arugments in last line means the input file ( `-i` ) is `mt1x.fasta`, the output file ( `-o` ) is `mt1x.wmsa.fasta`, use `16` threads ( `-T` ) and the similarity of cd-hit ( `-c` ) in  `wmsa `is `0.9`.

For dataset from www.drive5.com/bench, see [here](https://github.com/malabz/WMSA-dataset/blob/main/benchmark/README.md) for testing in `wmsa` .

## How to interpret the results

For `mt` and `SARS-COV-2` test case, we use `SP Score` to measure the result. The `SP Score` script can be found [here](https://github.com/malabz/MSATOOLS/tree/main/SPscore).

Use the script test the result by changing the arguments with modifing the match ( `--match` ) score equals to 1 and all mismatch conditions (like mismatch `--mismatch` , gap with character `--gap1` and gap with gap `--gap2` ) score equals to 0:

```bash
wget https://github.com/malabz/MSATOOLS/raw/main/SPscore/SPscore.py
python3 SPscore.py --input mt1x.wmsa.fasta --match 1 --mismatch 0 --gap1 0 --gap2 0
```

## How to upgrade this program

```bash
git pull
git submodule foreach 'git pull'
```

...or you can download source code on the release.

## How to remove this program

```bash
make uninstall
```
