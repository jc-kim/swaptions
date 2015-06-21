# swaptions
PARSEC Benchmark Swaption's OpenCL, MPI, and SnuCL Implementation

## Make
```make version=<option>```
### Options list
* pthreads: make file with Pthreads.
* cpu: make file with opencl implementation, using CPU.
* gpu: make file with opencl implementation, using GPU.
* mpi: make file with opencl and MPI implementation, using GPU.
* snucl: make file with SnuCL implementation, using GPU.

## use
When you execute makefile, you can see differnt file by compile option.

```./swaptions_<compile_option> -ns <swaptions_number> -sm <simulation_number> -nt <thread_number>```

### Options list
* swaptions_number: number of swaptions. default value is 1.
* simulation_number: number of simulation for each swaption. default value is 102400.
* thread_number. number of threads. it's only avaiable with pthread compile.
