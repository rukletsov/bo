
Summary
--------
Bö is a collection of basic and advanced methods and data structures for image processing and 3D reconstruction. Written entirely in C++, the library provides loosely-coupled though completely compatible modules for various problems in computer vision. Bö is [BSD licensed] (http://opensource.org/licenses/bsd-license.php).

The code has been compiled with, tested and run on 
- Ubuntu Linux with gcc, Fedora Linux with gcc;
- Mac OS X Lion (10.7) with gcc and clang;
- Windows XP with msvc2005, Windows XP with msvc2008, Windows 7 with msvc2010.

Project Dependencies
---------------------
- [Boost] (http://www.boost.org/). Minimal required version of boost is 1.43 because of `boost::array::fill()`. This function is used only once and can be rewritten, which lower the version needed to 1.38 because of the ScopeExit library. If you also want to build tests, boost version 1.44 is required because of the Filesystem version 3 library.
- [Google Test] (http://code.google.com/p/googletest/) is used for testing. If you are not going to build tests, you may ignore this dependency. 
- There is some header-only utility code for [OpenCV library] (http://sourceforge.net/projects/opencvlibrary/), however, it is disabled by default and if you don't plan to use it your project, you may ignore this dependecy.

Getting Started
----------------
### Getting Bö
No pre-built binaries are currently available. Use the following command to get the source code:
    $ hg pull https://bitbucket.org/rukletsov/b

### Building Bö
CMake is currently used to build Bö and its tests. If you use cmake from command line, you may find the following options useful:
    -DCMAKE_BUILD_TYPE=Debug -DBuildTests=ON -DUsePch=OFF -DBUILD_SHARED_LIBS=ON
    -DCMAKE_BUILD_TYPE=Release -DBuildTests=ON -DUsePch=OFF -DBUILD_SHARED_LIBS=OFF
