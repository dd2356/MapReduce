# MapReduce
The final project of DD2356, implementing a word counter using MapReduce

# Building 
Building can be done using `make`. This will generate a `mapreduce.out` binary under `bin/`. 
## Building in beskow
When building in beskow, the `CC = mpic++` line of the `Makefile` needs to be changed to `CC = CC`. Additionally, the following modules need to be loaded or swapped: 

* PrgEnv-gnu
* craype/2.5.14
* cce/8.6.5
* cray-mpich/7.7.0


# Running

The executable takes two parameters, the file to process and the location to place the result files in (plus any additional prefix). Running the following on beskow:
`aprun -n 8 ./bin/mapreduce.out file.in filepath/out`
will for instance yeild two output files, `filepath/out.8.bench` containing performance data and `filepath/out.8.words` containing the results. 
