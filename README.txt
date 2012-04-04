
Minimal required version of boost is 1.43 because of boost::array::fill(). 
By rewriting this function you can lower the version needed, but not below 
1.38 because of the ScopeExit library. If you also want to build tests,
boost version 1.44 is needed because of the Filesystem version 3 library. 
You may try to lower the version by changing the code to use version 2 of
the Filesystem library, however that has not been tested.

BitBucket cmds:
hg pull https://bitbucket.org/rukletsov/b
hg push https://bitbucket.org/rukletsov/b

CMake command-line options:
-DCMAKE_BUILD_TYPE=Debug -DBuildTests=ON -DUsePch=OFF -DBUILD_SHARED_LIBS=ON
-DCMAKE_BUILD_TYPE=Release -DBuildTests=ON -DUsePch=OFF -DBUILD_SHARED_LIBS=OFF
