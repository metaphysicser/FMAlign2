# FMAlign2: A novel fast multiple nucleotide sequence alignment method for ultra-long datasets

[FMAlign2](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btae014/7515251?searchresult=1) is a novel multiple sequence alignment algorithm based on [FMAlign](https://github.com/iliuh/FMAlign). It is designed to efficiently align ultra-long nucleotide sequences fast and accurately.

## Table of Contents
- [Installation](#Installation)
- [Usage](#Usage)
- [Data](#Data)
- [Issue](#Issue)
- [Related](#Related)
- [Citation](#Citation)
- [License](#License)

## Installation

The program is supported both on Linux and Windows(**Linux is strongly recommended for its convenience and better performance**). Please make sure your computer meets the following requirements:

- To compile the executable program for the entire project, **please ensure that you have the `make` command available.**

  - Verify `make` availability: Open your command-line interface and type ```make --version``` to check if the `make` command is installed on your system. If it is installed, you will see the version information. If not, you need to install it before proceeding.
  - Install `make` on Windows: If you are using Windows, you may need to install the appropriate tool to enable `make` functionality. One popular option is GNU Make for Windows , which provides a Windows-compatible version of `make`. You can download it from the [official website](https://www.gnu.org/software/make/) and follow the installation instructions.

- To compile the project, **you need to have the `g++` compiler available on your system.** Here are the steps to ensure `g++` support:

  - Check `g++` availability: Open your command-line interface and type `g++ --version` to check if the `g++` command is installed. If it is installed, you will see the version information. If not, you need to install it before proceeding.

  - Install `g++` on Windows: If you are using Windows, you can install `g++` by using a compiler suite such as MinGW or Cygwin. These packages provide a Windows-compatible version of `g++` along with other essential tools. You can download MinGW from the official website (https://mingw-w64.org/) or Cygwin from their official website (https://www.cygwin.com/). Follow the installation instructions provided by the respective package to set up `g++` on your system.

  - Install `g++` on Linux: On most Linux distributions, the `g++` compiler is included as part of the GNU Compiler Collection (GCC). To install `g++`, open your terminal and run the following command:

    ```
    sudo apt-get install g++
    ```

    This will install `g++` and its dependencies on your system.

  - Verify `g++` installation: After installation, run `g++ --version` again to verify that `g++` is installed correctly and accessible from the command line. Please note that if you are a Windows user, make sure that the installed version(>4.2) of `g++` supports OpenMP. On Windows systems, we utilize OpenMP for parallel computing.

---

If you have ensured that your system meets the requirements mentioned above, you can proceed with the following steps to compile the executable file. **However, you also have the option to directly use the pre-compiled executable file available in the [Release](https://github.com/metaphysicser/FMAlign2/releases/tag/v1.0).**

1. **DownLoad**

   ```shell
   git clone https://github.com/metaphysicser/FMAlign2.git
   cd FMAlign2
   # for Linux
   chmod 777 ./ext/mafft/linux/usr/libexec/mafft/disttbfast
   ```

2. **Build**

   ```shell
   cd FMAlign2 && make [M64=1]
   ```

   Switch to the FMAlign2 directory in your terminal and execute the above command to build the project. We provide two compilation modes: 32-bit and 64-bit. In most cases, the 32-bit mode is sufficient to handle most data. However, if the concatenated length of all sequences exceeds the range of uint32_t (4294967295), you should add the M64 parameter when compiling the program to generate a 64-bit executable.

   - If you don't need the 64-bit mode, simply execute the `make` command.
   - If you need the 64-bit mode, execute the `make M64=1` command.

   During the compilation process, please be patient as the time required depends on the size and complexity of the project.

   Once the compilation is complete, you will find the generated executable file in the specified output directory.

   Note: If you want to remove all the generated `.o` files, you can execute the following command:
   ```shell
   make clean
   ```
   This command will clean up the intermediate object files and leave only the source code and executable file in the project directory. Use this command when you want to start a fresh build or clean up unnecessary files to save disk space.

---

**Please note that if you choose halign2 and halign3 as your multiple sequence alignment methods, make sure you have Java environment installed.** To check the version of Java installed on your system, you can open a command prompt or terminal and execute the following command:

```bash
java -version
```

This will display the installed Java version information.

If you don't have Java installed or if the installed version is not compatible, you can follow these steps to install Java:

To install Java on Windows:

1. Visit the official Java website at [java.com](https://www.java.com/) or the OpenJDK website at [openjdk.java.net](https://openjdk.java.net/).
2. Download the appropriate Java Development Kit (JDK) for Windows.
3. Run the downloaded installer and follow the on-screen instructions to complete the installation.
4. After the installation is complete, open a new Command Prompt and run `java -version` to verify that Java is installed and the correct version is displayed.

To install Java on Linux:

1. Update Package Lists: Run the command `sudo apt update` to update the package lists on your system.
2. Install OpenJDK: Run the command `sudo apt install default-jdk` to install the default version of OpenJDK.
3. Verify Installation: After the installation is complete, run `java -version` to verify that Java is installed and the correct version is displayed.

Once you have Java installed and verified the version, you should be able to use halign2 and halign3 for multiple sequence alignment.


## Usage

> Reminder: Please ensure that all external files (such as MAFFT, HALIGN, etc.) are properly copied to their corresponding directories. Pay close attention to the relative paths between FMAlign2 and the ext folder to avoid issues during execution,

if you are Linux user:

   ```shell
   ./FMAlign2 -i /path/to/data [other options]
   ```

if you are Windows user:

   ```shell
   ./FMAlign2.exe -i /path/to/data [other options]
   ```

if you want to show the parameters details:

 ~~~
 ./FMAlign2 -h
 ~~~

  Parameters Details:

   - -i [file path] **[required]** The path to the input file.
   - -o [output path] [default: ouput.fmaligned2.fasta] The path to the output file.
   - -p [package] [default: mafft] The MSA method used in parallel align. for example, [**halign3**](https://github.com/malabz/HAlign-3), [**halign2**](https://github.com/ShixiangWan/HAlign2.0) and [**mafft**](https://mafft.cbrc.jp/alignment/software/).
   - -t [int] [default: cpu number]  The maximum number of threads that the program runs, the recommended setting is the number of CPUs.
   - -l [int] [default: square root of mean length] The minimum length of MEMs, the default value is square root of mean length.
   - -c [float] [default: 1] A floating-point parameter that specifies the minimum coverage across all sequences, with values ranging from 0 to 1.
   - -f [mode] [default: global or local] The filter MEMs mode. The default setting is that if sequence number less 100, **accurate** mode otherwise **global** mode.
   - -d [int] [default:0] Depth of recursion, you could ignore it.
   - -v [int] [default:1] Verbose option, 0 or 1. You could ignore it.
   - -h [help] print help information

-----

We will demonstrate with the example data `mt1x.fasta`, assuming you are running on a Linux system.

```shell
./FMAlign2 -i ./data/mt1x.fasta -l 20 -c 1 -p mafft -f gloabl -o output.fmaligned2.fasta
```

This command specifies the following options:
- Input data: `mt1x.fasta` located in the `data` folder.
- Minimum length of MEMs: 20.
- Sequence coverage of MEMs: 1.
- Parallel alignment method: mafft
- Alignment mode: global mode.
- Output file: `output.fmaligned2.fasta` will be generated in the FMAlign2 directory.

After running this command, you will obtain the aligned output in the `output.fmaligned2.fasta` file.

-----

If you want to evaluate the generated alignment results, you can run the `sp.py` script (requires a Python environment) with the following parameters:

```shell
python sp.py --input output.fmalign2.fasta --match 0 --mismatch 1 --gap1 2 --gap2 2
```

This command will calculate and print the SP (Sum-of-Pairs) score for the multiple sequence alignment results. The `--input` parameter specifies the input alignment file (`output.fmalign2.fasta` in this case), and the `--match`, `--mismatch`, `--gap1`, and `--gap2` parameters define the scoring scheme for matches, mismatches, and gap penalties.

By running this command, you will obtain the SP score, which provides an evaluation of the alignment quality.

## Data

Data can be assessed in [data](https://github.com/metaphysicser/FMAlign2/tree/master/data) fold. All the data is compressed using xz compression. Before using it, please decompress the files.

Here are the methods to decompress the files on different operating systems:

**Decompressing on Linux:**

1. Open the terminal.

2. Navigate to the directory where the compressed file is located.

3. Run the following command to decompress the file:

   ```
   xz -d filename.xz
   ```

   Replace `filename.xz` with the name of the file you want to decompress.

**Decompressing on Windows:**

1. Download and install an xz compression tool for Windows, such as [7-Zip](https://www.7-zip.org/) or [WinRAR](https://www.win-rar.com/).
2. Right-click on the compressed file.
3. Select "Extract to" or a similar option to decompress the file.

Please note that the decompressed files will occupy more disk space. Make sure you have enough disk space to store the uncompressed files.

If you need more data, you can visit http://lab.malab.cn/~cjt/MSA/datasets.html for more datasets.

## Issue

FMAlign2 is supported by [ZOU's Lab](https://github.com/malabz). If you have any suggestions or feedback, we encourage you to provide them through the issue page on the project's repository. You can also reach out via email to zpl010720@gmail.com.

We value your input and appreciate your contribution to improving the project. Thank you for taking the time to provide feedback, and we will address your concerns as soon as possible.

## Related

- [FMAlign](https://github.com/iliuh/FMAlign): a fast multiple nucleotide sequence alignment method based on FM-index
- [HAlign3](https://github.com/malabz/HAlign-3) and [HAlign2](https://github.com/ShixiangWan/HAlign2.0)
- [WMSA](https://github.com/malabz/WMSA) and [WMSA2](https://github.com/malabz/WMSA2)
- [TPRA](https://github.com/malabz/TPRA): A refinement tool for ensembling different multiple sequence alignment results
- [MSATOOLS](https://github.com/malabz/MSATOOLS): Some tools for MSA, like SP score, Q score and so on.

## Citation

Pinglu Zhang, Huan Liu, Yanming Wei, Yixiao Zhai, Qinzhong Tian, Quan Zou, FMAlign2: a novel fast multiple nucleotide sequence alignment method for ultralong datasets, *Bioinformatics*, 2024;, btae014, https://doi.org/10.1093/bioinformatics/btae014

## License

[Apache 2.0](https://github.com/metaphysicser/FMAlign2/blob/master/LICENSE) © [[MALABZ_UESTC](https://github.com/malabz) [Pinglu Zhang](https://github.com/metaphysicser)]
