DEF =
INCLUDE =

EXEC = swaptions 
ALL_EXEC = swaptions run_cpu run_gpu run_mpi run_snucl

OBJS= CumNormalInv.o MaxFunction.o RanUnif.o nr_routines.o \
	HJM_SimPath_Forward_Blocking.o HJM.o HJM_Swaption_Blocking.o  \
	HJM_Securities.o

ifdef version
  ifeq "$(version)" "pthreads" 
    DEF := $(DEF) -DENABLE_THREADS
    CXXFLAGS := $(CXXFLAGS) -pthread
	endif
	ifeq "$(version)" "cpu"
		DEF := $(DEF) -DENABLE_OPENCL -DDEVICE_CPU
		CXXFLAGS := $(CXXFLAGS) -lOpenCL
		EXEC := run_cpu
	endif
	ifeq "$(version)" "gpu"
		DEF := $(DEF) -DENABLE_OPENCL -DDEVICE_GPU
		CXXFLAGS := $(CXXFLAGS) -lOpenCL
		EXEC := run_gpu
	endif
	ifeq "$(version)" "mpi"
		CXX := mpig++
		DEF := $(DEF) -DENABLE_OPENCL -DENABLE_MPI
		EXEC := run_mpi
	endif
	ifeq "$(version)" "snucl"
		DEF := $(DEF) -DENABLE_OPENCL -DENABLE_SNUCL
		CXXFLAGS := $(CXXFLAGS) -lOpenCL
		EXEC := run_snucl
	endif
endif

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(DEF) $(OBJS) $(INCLUDE) $(LIBS) -o $(EXEC)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEF) -c $*.cpp -o $*.o

.c.o:
	$(CXX) $(CXXFLAGS) $(DEF) -c $*.c -o $*.o

clean:
	rm -f $(OBJS) $(ALL_EXEC)

