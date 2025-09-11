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


---



### Integrating Custom MSA Methods with FMAlign2

FMAlign2 can be combined with any MSA method as long as the input and output are in FASTA format. For example, if you want to run MAFFT with additional parameters, you can create a text file (e.g., `mafft-cmd.txt`) containing a command in the following format:

```
mafft --retree 1 --thread {thread} {input} > {output}
```

Here, `{input}` and `{output}` are placeholders for the input and output FASTA files, respectively. The `{thread}` placeholder is optional — if omitted, the MSA tool will fall back to its own default number of threads. Compared to the default MAFFT command, this example adds `--retree 1` as a custom parameter to better suit certain use cases.

Importantly, the same mechanism works with any MSA tool, not just MAFFT.

To use a custom command file, simply pass the path to it with the `-p` parameter. For example:

```
FMAlign2 -i /path/to/input.fasta -o /path/to/output.fasta -p /path/to/mafft-cmd.txt
```

### Evaluation
If you want to evaluate the generated alignment results, you can run the `sp.py` script (requires a Python environment) with the following parameters:

```shell
python sp.py --input output.fasta --match 0 --mismatch 1 --gap1 2 --gap2 2
```

This command will calculate and print the SP (Sum-of-Pairs) score for the multiple sequence alignment results. The `--input` parameter specifies the input alignment file (`output.fasta` in this case), and the `--match`, `--mismatch`, `--gap1`, and `--gap2` parameters define the scoring scheme for matches, mismatches, and gap penalties.

By running this command, you will obtain the SP score, which provides an evaluation of the alignment quality.

### Installation from Source (detailed)

You can build FMAlign2 from source on Linux and Windows (MSYS2/MinGW). Below are step-by-step instructions, optional flags, and install targets.

#### 0) Prerequisites

* **Build tools**: `g++` (GCC ≥ 9), `make`
* **Recommended runtime tools** (pick any you plan to use):

  * **MAFFT** (for `-p mafft`), writes alignment to **stdout**
  * **HAlign3 / HAlign4** (for `-p halign3` / `-p halign4`)
  * **OpenJDK 11** (only needed if you use the HAlign JAR fallback)

Install examples:


# Installation Guide

## Ubuntu/Debian

```bash
sudo apt update
# Optional runtime dependency
sudo apt install -y mafft
# Or via conda
# conda install -c conda-forge -c bioconda mafft halign openjdk=11
```

---

## Install HAlign3 (JAR)

### System-wide installation (with sudo)

```bash
# Download and move to /usr/local/bin
wget https://github.com/malabz/HAlign-3/releases/download/v3.0.0-rc1/HAlign-3.0.0_rc1.jar
sudo mv HAlign-3.0.0_rc1.jar /usr/local/bin/

# Create a wrapper script
sudo tee /usr/local/bin/halign3 >/dev/null <<'EOF'
#!/usr/bin/env bash
exec java -jar /usr/local/bin/HAlign-3.0.0_rc1.jar "$@"
EOF
sudo chmod +x /usr/local/bin/halign3

# Test
halign3 -h
```

### **User installation (no sudo)**

```bash
# Install under $HOME/bin
mkdir -p $HOME/bin
wget -O $HOME/bin/HAlign-3.0.0_rc1.jar \
  https://github.com/malabz/HAlign-3/releases/download/v3.0.0-rc1/HAlign-3.0.0_rc1.jar

# Create a wrapper script
cat > $HOME/bin/halign3 <<'EOF'
#!/usr/bin/env bash
exec java -jar $HOME/bin/HAlign-3.0.0_rc1.jar "$@"
EOF
chmod +x $HOME/bin/halign3

# Add $HOME/bin to PATH if not already
export PATH=$HOME/bin:$PATH

# Test
halign3 -h
```

---

## Install HAlign4 (C++)

### System-wide installation (with sudo)

```bash
git clone https://github.com/metaphysicser/HAlign-4.git
cd HAlign-4
make -j
sudo install -m 0755 halign4 /usr/local/bin/

# Test
halign4 -h
```

### **User installation (no sudo)**

```bash
git clone https://github.com/metaphysicser/HAlign-4.git
cd HAlign-4
make -j

# Move to user bin directory
mkdir -p $HOME/bin
cp halign4 $HOME/bin/

# Add $HOME/bin to PATH if not already
export PATH=$HOME/bin:$PATH

# Test
halign4 -h
```

---

## Summary

* **System-wide installation**: requires `sudo`, binaries are placed under `/usr/local/bin`.
* **User installation**: no `sudo` required, binaries go under `$HOME/bin`, make sure `$HOME/bin` is added to your `PATH`.

After installation, both `halign3` and `halign4` will be available in your system `PATH` and can be directly used with FMAlign2.

---


#### 1) Clone and build

```bash
git clone https://github.com/metaphysicser/FMAlign2
cd FMAlign2

# Build (default: optimized; Linux defaults to static linking if available)
make -j
```

Optional flags:

* `DEBUG=1` → add `-O0 -g -DDEBUG`
* `M64=1` → define `-DM64` (and on x86\_64 adds `-m64`)
* `STATIC_LINK=0` → dynamic linking (recommended for most users)

Examples:

```bash
make STATIC_LINK=0 M64=1
make DEBUG=1
```

---

#### 2) Install (optional)

The Makefile provides `install`/`uninstall` targets:

```bash
# Install to /usr/local/bin (default PREFIX)
sudo make install

# Custom prefix (e.g., /opt/fmalign2)
make install PREFIX=/opt/fmalign2

# Uninstall (use the same PREFIX/DESTDIR you installed with)
sudo make uninstall
```

This installs the `fmalign2` binary; you can then run `fmalign2 -h` anywhere on your system.

---

#### 3) Verify

```bash
./fmalign2 -h
# or after install:
fmalign2 -h
```
---

#### 4) Choosing an MSA backend

FMAlign2 can use:

* **MAFFT** (`-p mafft`) — **writes to stdout** → FMAlign2 redirects to your `-o` file
* **HAlign3/HAlign4** (`-p halign3`, `-p halign4`) — support `-o <file>` natively
* **Custom command file** (`-p /path/to/cmd.txt`) — must include `{input}`, `{output}`, and optional `{thread}` placeholders




## Issue

FMAlign2 is supported by [ZOU's Lab](https://github.com/malabz). If you have any suggestions or feedback, we encourage you to provide them through the issue page on the project's repository. You can also reach out via email to zpl010720@gmail.com.

We value your input and appreciate your contribution to improving the project. Thank you for taking the time to provide feedback, and we will address your concerns as soon as possible.


## Citation

Pinglu Zhang, Huan Liu, Yanming Wei, Yixiao Zhai, Qinzhong Tian, Quan Zou, FMAlign2: a novel fast multiple nucleotide sequence alignment method for ultralong datasets, *Bioinformatics*, 2024;, btae014, https://doi.org/10.1093/bioinformatics/btae014

## License

[Apache 2.0](https://github.com/metaphysicser/FMAlign2/blob/master/LICENSE) © [[MALABZ_UESTC](https://github.com/malabz) [Pinglu Zhang](https://github.com/metaphysicser)]
