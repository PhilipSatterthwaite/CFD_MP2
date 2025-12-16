Hi Izzy!

Compiling:
    mp2.f90 can be compiled with gfortran using the following command:
    gfortran precision.f90 compressible.f90 mp2.f90 -llapack -lopenblas -O3 -o mp2

Running:
    mp2.exe can be run with inputs formatted as follows:
    ./mp1 [x] [t] [end][output]
    [x] = number of grid points (in both directions)
    [t] = time step size
    [end] = end time t
    [output] = output directory and file name where the .txt file of data will be saved
    example: ./mp2 100 0.001 1.0 ./test.txt

Code structure:
    mp1.f90 reads in the input variables and calls scheme1.f90 for scheme 1 
    and scheme2.f90 for scheme 2.
    Note: 
    The solution matrix is allocated before the run (not written to file during the run),
    so extremely small spatial and temporal resolution or runs to extremely long end times
    will cause a memory error.

Assignment organization:
    Running mp1.py will sequentially run 4 python functions that call mp1.exe as necessary 
    and produce the plots used in the the final submission.
    Each function can also be run individually by commenting out necessary lines.

