
# Summary
--------------------------------------------------------------------------------
Bö is a collection of basic and advanced methods and data structures for 2D image processing and segmentation, 3D object analysis, surface reconstruction. Written entirely in C++, the library provides loosely-coupled though completely compatible modules for various problems in computer vision. Bö is [BSD licensed](http://opensource.org/licenses/bsd-license.php). Please be advised that Bö is mainly a research project and might contain bugs and rough edges.

## Main Features
 * Highly customizable [Markov random field (MRF)](http://en.wikipedia.org/wiki/Markov_random_field) representation for regular 2D lattices, as well as several predefined prior and likelihood energy functionals ready to use in MRF models.
 * Optimization algorithms for Markov random fields: [Metropolis–Hastings algorithm](http://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm), [Modified Metropolis Dynamics (MMD)](http://www.sciencedirect.com/science/article/pii/0262885695010726), [Iterated conditional modes (ICM)](http://en.wikipedia.org/wiki/Iterated_conditional_modes).
 * [Modified Butterfly](http://mrl.nyu.edu/~dzorin/papers/zorin1996ism.pdf) subdivision surface method.
 * Robust [mesh-growing surface reconstruction](http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1353243) method.
 * [Iterative closest point (ICP)](http://en.wikipedia.org/wiki/Iterative_closest_point) registration algorithm.
 * Triangular mesh with support of export/import to [PLY file format](http://en.wikipedia.org/wiki/PLY_%28file_format%29).
 * Various point-to-point and point-to-plane metrics in 3D.
 * Transformations in 3D space based on quaternions or transformation matrices.
 * Generic representation of multidimensional vectors and points.
 * Generic representation of 2D images.
 * Linear filtering for 2D images.
 * Matrix operations (e.g. matrix inversion, eigenvectors and eigenvalues of real matrices).

## C++11 
No C++11 features are currently used in order to support older compilers.

# Project Dependencies
--------------------------------------------------------------------------------
 * [Boost](http://www.boost.org/). Minimal required version of boost is 1.43 because of `boost::array::fill()`. This function is used only once and can be rewritten, that lowers the version needed to 1.38 because of the ScopeExit library. If you also want to build tests, boost version 1.44 is required because of the Filesystem version 3 library.
 * [Google Test](http://code.google.com/p/googletest/) is used for testing. If you are not going to build tests, you may ignore this dependency. 
 * There is some header-only utility code for [OpenCV library](http://sourceforge.net/projects/opencvlibrary/), however, it is disabled by default and if you don't plan to use it in your project, you may ignore this dependency.

# Getting Started
--------------------------------------------------------------------------------
## Getting Bö
No pre-built binaries are currently available. Use the following command to get the source code:

    $ hg pull https://bitbucket.org/rukletsov/b

## Building Bö
CMake is currently used to build Bö and its tests. If you use CMake from command line, you may find the following options useful:

    -DCMAKE_BUILD_TYPE=Debug -DBuildTests=ON -DUsePch=OFF -DBUILD_SHARED_LIBS=ON
    -DCMAKE_BUILD_TYPE=Release -DBuildTests=ON -DUsePch=OFF -DBUILD_SHARED_LIBS=OFF

The code has been successfully built with, tested and run on:

 * Ubuntu 10.04 32bit with gcc 4.4.3.
 * Fedora 17 32bit with gcc 4.7.0.
 * openSUSE 12.2 64bit with gcc 4.6.2 and clang 3.1.
 * Mac OS X Lion (10.7.2) 64bit with gcc 4.2.1 and clang 3.0.
 * Windows XP 32bit with msvc2005 and msvc2008.
 * Windows 7 64bit with msvc2010.

## Using Bö
Code samples and overview of some features are available on [project wiki](https://bitbucket.org/rukletsov/b/wiki/Home).
