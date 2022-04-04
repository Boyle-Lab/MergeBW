# Installation
Clone a copy of the repository and submodules:

```
git clone --recurse-submodules https://github.com/Boyle-Lab/MergeBW.git
```

Build bamtools API (please see bamtools documentation for more information)
Note: bamtools requires zlib to be installed
```
cd MergeBW/bamtools/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=$(cd ..; pwd)/install ..
make
make install
cd ../..
```

Build 
```
make
```

