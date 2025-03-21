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

## Visualization

If the project has been built in debug, applying any algorithm will
create a visualization file in the directory visualizations.
This can be loaded in the visualizer. To start it use the following command.

```bash
> python vis.py
```

TODO: add required python libraries

After that open a browser at ```127.0.0.1:8050```.
This will (if the python program is running)
open the visualization with a standard polyline.
To load the visualization click on the respective button
and navigate to the visualizations directory and open the one you want.
On the top you will see the data flow tree
that shows the dependence of the entries of the dynamic program.
You can click on a node in the tree the show the visualization for
that respective entry. This is shown on the bottom of the page and
shows the polyline with the respective simplification up to this point
and the part of the polyline that is traversed.
The other button returns you to the original view that shows the polyline
and the simplification. By clicking on a point on the polyline
a circle/square/diamond is drawn around it with the specified radius of
epsilon to easily find intersections with the nearby line segments.
