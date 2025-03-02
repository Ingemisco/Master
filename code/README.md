# Polyline Simplification Algorithms

## About

 TODO

## Prerequisites

The following dependencies are required to build and run this project

- Boost

```bash
To install on Arch based distributions 
> sudo pacman -S boost

To install on Debian based distributions
> sudo apt install boost
```

## Build and Run

To build the project CMake and a C++ Compiler that supports C++23
(GCC 13+, Clang 16+, ...) are required as well as all additional
libraries mentioned in [Prerequisites](#prerequisites).

Clone this project and navigate to the `code` directory.
Create a `build` directory and navigate into it:

```bash
> mkdir build 
> cd build 
```

After that you can build the project in debug mode:

```bash
> cmake .. -DCMAKE_BUILD_TYPE=Debug
> cmake --build .
```

Or in Release mode:

```bash
> cmake .. -DCMAKE_BUILD_TYPE=Release 
> cmake --build .
```

Additionally, there is a test build (does not do anything as of yet):

```bash
> cmake .. -DCMAKE_BUILD_TYPE=Test 
> cmake --build .
```

## Usage

The debug and release build create two executables: polyline and datagen.
datagen allows creation of polyline data and store them in files.
It has various options which can be seen using the -h flag

```bash
> ./datagen -h
```

polyline can read polyline data files and perform algorithms on them.
For more information see the -h flag

```bash
> ./polyline -h
```
