Compression
===========
This is the MATLAB code accompanying the M.Eng. thesis "A Probabilistic Graphical Model Based Data Compression Architecture for Gaussian Sources" by Wai Lok Lai, Sep 2016, the pdf of which is also included in this repository.

This software is released under the MIT License (included with the software). Note, however, that using this code (and/or the results of running it) to support any form of publication (e.g., a book, a journal paper, a conference paper, a patent application, etc.) requires you to cite the aforementioned thesis:

```
@mastersthesis{lai2016compression,
  title={A Probabilistic Graphical Model Based Compression Architecture for Gaussian Sources},
  author={Lai, Wai Lok},
  year={2016},
  school={Massachusetts Institute of Technology}
}
```

Author:
-------

Wai Lok Lai (email: anguslai@alum.mit.edu)


Requirements:
-------------
- MATLAB 2011a or later. Not tested with MATLAB 2010 or before. 
- A MATLAB recognized C compiler. 

Getting Started:
----------------

- Compile the C-MEX files with MATLAB: from the root of the git repository, run the following commands in the MATLAB console:

```
cd lib/ldpc
mex -largeArrayDims generate_ldpc.c rand.c alloc.c intio.c open.c mod2sparse.c mod2dense.c mod2convert.c rcode.c distrib.c
cd ../..
```

- Running the Code: The entry point of the code is `./main.m`.

Folder structure: 
-----------------
Files:
`./main.m`: main entry point of the code
`./randfile`: used for generating randomness given a seed
`./startup.m`: initializes the environment for running our code in MATLAB
`./thesis.pdf`: the pdf of the aforementioned thesis

Folders:
`./decode`: contains functions for decoding (decompression)
`./encode`: contains functions for encoding (compression)
`./lib`: contains libraries that the code uses
`./source`: contains the definitions and BP implementations of different source models
