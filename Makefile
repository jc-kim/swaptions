DEF =
INCLUDE =
CXXFLAGS := $(CXXFLAGS) -g -Wall -O4

EXEC = swaptions 
ALL_EXEC = swaptions swaptions_cpu swaptions_gpu swaptions_mpi swaptions_snucl

OBJS= CumNormalInv.o MaxFunction.o RanUnif.o nr_routines.o \
	HJM_SimPath_Forward_Blocking.o HJM.o HJM_Swaption_Blocking.o  \
	HJM_Securities.o timers.o

ALL_OBJS= $(OBJS) HJM_Securities_cpu.o HJM_Securities_gpu.o \
					HJM_Securities_mpi.o HJM_Securities_snucl.o

ifdef version
  ifeq "$(version)" "pthreads" 
    DEF := $(DEF) -DENABLE_THREADS
    CXXFLAGS := $(CXXFLAGS) -pthread
	else
		ifeq "$(version)" "cpu"
			LDFLAGS := $(LDFLAGS) -lOpenCL
			DEF := $(DEF) -DENABLE_OPENCL -DDEVICE_CPU
			OBJS := CumNormalInv.o MaxFunction.o RanUnif.o nr_routines.o \
				HJM_SimPath_Forward_Blocking.o HJM.o HJM_Swaption_Blocking.o  \
				HJM_Securities_cpu.o timers.o
			EXEC := swaptions_cpu
		endif
		ifeq "$(version)" "gpu"
			LDFLAGS := $(LDFLAGS) -lOpenCL
			DEF := $(DEF) -DENABLE_OPENCL -DDEVICE_GPU
			OBJS := CumNormalInv.o MaxFunction.o RanUnif.o nr_routines.o \
				HJM_SimPath_Forward_Blocking.o HJM.o HJM_Swaption_Blocking.o  \
				HJM_Securities_gpu.o timers.o
			EXEC := swaptions_gpu
		endif
		ifeq "$(version)" "mpi"
			LDFLAGS := $(LDFLAGS) -lOpenCL
			CXX := mpic++
			DEF := $(DEF) -DENABLE_OPENCL -DENABLE_MPI
			OBJS := CumNormalInv.o MaxFunction.o RanUnif.o nr_routines.o \
				HJM_SimPath_Forward_Blocking.o HJM.o HJM_Swaption_Blocking.o  \
				HJM_Securities_mpi.o timers.o
			EXEC := swaptions_mpi
		endif
		ifeq "$(version)" "snucl"
			CXX := mpic++
			DEF := $(DEF) -DENABLE_OPENCL -DENABLE_SNUCL
			OBJS := CumNormalInv.o MaxFunction.o RanUnif.o nr_routines.o \
				HJM_SimPath_Forward_Blocking.o HJM.o HJM_Swaption_Blocking.o  \
				HJM_Securities_snucl.o timers.o
			EXEC := swaptions_snucl
			INCLUDE := $(INCLUDE) -I$(SNUCLROOT)/inc
			LDLIBS := $(LDLIBS) -L$(SNUCLROOT)/lib 
			LDFLAGS := $(LDFLAGS) -lsnucl_cluster
		endif
		ifeq "$(version)" "snucl_no_reduc"
			CXX := mpic++
			DEF := $(DEF) -DENABLE_OPENCL -DENABLE_SNUCL -DNON_REDUCTION
			OBJS := CumNormalInv.o MaxFunction.o RanUnif.o nr_routines.o \
				HJM_SimPath_Forward_Blocking.o HJM.o HJM_Swaption_Blocking.o  \
				HJM_Securities_snucl_without_reduc.o timers.o
			EXEC := swaptions_snucl_no_reduc
			INCLUDE := $(INCLUDE) -I$(SNUCLROOT)/inc
			LDLIBS := $(LDLIBS) -L$(SNUCLROOT)/lib 
			LDFLAGS := $(LDFLAGS) -lsnucl_cluster
		endif
	endif
endif

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(DEF) $(OBJS) $(INCLUDE) $(LIBS) $(LDLIBS) -o $(EXEC)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEF) -c $*.cpp -o $*.o

.c.o:
	$(CXX) $(CXXFLAGS) $(DEF) -c $*.c -o $*.o

HJM_Securities_cpu.o: HJM_Securities.cpp
	$(CXX) $(CXXFLAGS) $(DEF) -c $^ -o $@

HJM_Securities_gpu.o: HJM_Securities.cpp
	$(CXX) $(CXXFLAGS) $(DEF) -c $^ -o $@

HJM_Securities_mpi.o: HJM_Securities.cpp
	$(CXX) $(CXXFLAGS) $(DEF) -c $^ -o $@

HJM_Securities_snucl.o: HJM_Securities.cpp
	$(CXX) $(CXXFLAGS) $(DEF) -c $^ -o $@

HJM_Securities_snucl_without_reduc.o: HJM_Securities.cpp
	$(CXX) $(CXXFLAGS) $(DEF) -c $^ -o $@

clean:
	rm -f $(ALL_OBJS) $(ALL_EXEC)
