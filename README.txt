
Boost minimal version is 1.43 because of boost::array::fill(). By rewriting 
this function you can lower the version needed, but not below 1.38 because
of the ScopeExit library.

BitBucket cmds:
hg pull https://bitbucket.org/rukletsov/b
hg push https://bitbucket.org/rukletsov/b

CMake command-line options:
-DCMAKE_BUILD_TYPE=Debug -DBuildTests=ON -DUsePch=OFF -DBUILD_SHARED_LIBS=ON
-DCMAKE_BUILD_TYPE=Release -DBuildTests=ON -DUsePch=OFF -DBUILD_SHARED_LIBS=OFF
