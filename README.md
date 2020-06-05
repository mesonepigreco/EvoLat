# EvoLat
Evolution on a Lattice simulation program

# Requirements

The requirements to compile the code are a C++ compiler (tested with GNU compilers, version 5) and the library libconfig++ with headers.
On Ubuntu linux, the following command should be able to install all the requirements.
```bash
>>> sudo apt-get install build-essential g++ libconfig++-dev
```


# Installation

The code can be compiled using the Makefile, with the standard procedure:
```bash
>>> make
```

You can edit the Makefile to include the libconfig++ headers installation path, and eventually change the path of the compiler.

# Running 

The code will be installed in the bin directory. It needs an inputfile in json format, as readed by the libconfig library.
An example input file can be found into the Test directory

To execute the code with the input.in input file, run
```bash
bin/wator.exe -i input.in 
```

