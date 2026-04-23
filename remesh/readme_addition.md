
## Dependencies, Installation Instructions, and Mex
The only dependency for this code is the library OpenMesh, specifically version 5.0, which can be found at https://www.graphics.rwth-aachen.de/software/openmesh/download/. Here we will provide step-by-step instructions for proper installation of OpenMeshand usage with MEX and the included Isotropic Remesher library (https://github.com/christopherhelf/isotropicremeshing).


1. Download OpenMesh 5.0. Extract it to anywhere you want.
2. Navigate in terminal to the location of the extracted folder (OpenMesh-5.0)
3. Use the following commands to cmake and install the library to the correct locations on your computer:
      1. mkdir build && cd build
      2. cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-std=c++11 -fPIC" -DBUILD_APPS=OFF -DOM_FORCE_STATIC_CAST=1
      3. sudo make install
4. Now that OpenMesh is properly installed, you can access it from matlab’s mex framework. Simply navigate to the isomesh folder, and run
      1. mex -I./ -I/usr/local/include -L/usr/local/lib -lOpenMeshTools -lOpenMeshCore LDFLAGS='\$LDFLAGS -Wl,-rpath,/usr/local/lib' remeshing.c ./src/BSP.cpp ./src/BSPTraits.cpp ./src/IsotropicRemesher.cpp
5. Now you can use the remesher per the instructions laid out in the library’s github, i.e.
      1. [facetsOut, pointsOut] = remeshing(facets, points, features, targetedgelength, iterations); 
      2. note: features can simply be int32([]). You also will need to int32([facets]) and int32([iterations])
6. On some systems, a version error may arise after using mex, saying 
      1. Invalid mex file, libstdc++.so.6: version GLIBCXX_3.4.29' not found
      2. In this case, the fix is to boot matlab from the terminal with the command
      3. export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
      4. matlab
