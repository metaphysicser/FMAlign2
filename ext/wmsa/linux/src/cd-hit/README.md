# For cd-hit

## Module

This is the part of `cd-hit` modified by [gitee@wym6912](https://gitee.com/wym6912) / [github@wym6912](https://github.com/wym6912), which is a module of the other program. This module only have these parts in below:

```bash
cd-hit
cd-hit-est
cd-hit-454
```

## Requirements
Since 4.8.1, cd-hit supports `gz` format input file. This requires `zlib` library. `zlib` should
be install in most Linux systems, so cd-hit should be compiled without issue. If your system
don't have `zlib`, please install it first. 

- On Ubuntu, to install `zlib`:
```bash
sudo apt install zlib1g-dev
```
- On CentOS, to install `zlib`:
```bash
sudo yum install zlib-devel 
```

## How to compile
1. Compile with multi-threading support (default):  `make`
2. Compile without multi-threading support (if you are on very old systems): `make openmp=no`
3. Compile without `zlib` (if you can not install `zlib`): `make zlib=no`

Having problems to compile, please contact the author.


## Compile cd-hit on MacOS
To install CD-HIT on MacOS, first install `gcc` on your system.
To use [Homebrew](https://brew.sh/), see [https://formulae.brew.sh/formula/gcc@6](https://formulae.brew.sh/formula/gcc@6). Then locate the path to your `g++` executable, (e.g. `/usr/local/Cellar/gcc/6.3.0_1/bin/g++-6`, note: yours `g++` path is likely to be different), then use command like this:

```bash
make CC=/usr/local/Cellar/gcc/6.3.0_1/bin/g++-6
```
For more information, please visit http://cd-hit.org.

Most up-to-date documents are available at https://github.com/weizhongli/cdhit/wiki.

cd-hit is also available as web server, visit http://cd-hit.org for web server address.
