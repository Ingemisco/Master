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

To build the project CMake and a C++ Compiler that supports C++23 (GCC >= 13, Clang >= 16, ...) are required
as well as all additional libraries mentioned in [Prerequisites](#prerequisites).

Clone this project and navigate to the `code` directory.
Create a `build` directory and navigate into it:

```bash
> mkdir build 
> cd build 
```

After that you can build the project:

```bash
> cmake ..
> cmake --build .
```

By default, the debug build will be created. To create the
release build (i.e., without debug flags and with optimizations)
do the following:

```bash
> cmake .. -DCMAKE_BUILD_TYPE=Release 
> cmake --build .
```

## Usage

TODO
