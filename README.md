# Installation
Clone a copy of the repository and submodules:

```
git clone --recurse-submodules https://github.com/Boyle-Lab/MergeBW.git
```

Build libBigWig
```
cd MergeBW/libBigWig/
make
cd ..
```

Build 
```
make
```

# Use
Input path to the input folder with bigwig files and the text file containing the name and size of chromosomes in command line: 
(note: the text file should be in the same folder as the bigwig files)

```

./MergeBW /path/inputfolder/ filename
```
