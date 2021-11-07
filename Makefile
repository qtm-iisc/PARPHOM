SRC_DIR = ./SRC/
EXAMPLE_DIR = ./Examples/

MODS=$(wildcard $(SRC_DIR)mod*.f90)
MODOBJS=$(patsubst %.f90, %.o,$(MODS))

SRCS=$(filter-out $(MODS), $(wildcard $(SRC_DIR)*.f90))
OBJS=$(patsubst %.f90, %.o,$(SRCS))

EXE = PhonFreq

FC = h5pfc
FLFLAGS = -L${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a  \
		  -L${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
		  -L${MKLROOT}/lib/intel64 -L${MKLROOT}/lib/intel64/libmkl_avx512.so \
		  -L${MKLROOT}/lib/intel64/libmkl_def.so \
		  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
		  -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl
FCFLAGS = -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
		  

all: $(MODOBJS) $(OBJS)
	$(FC) $(FCFLAGS) -o $(EXE) $(OBJS) $(MODOBJS) $(FLFLAGS)

$(MODOBJS): %.o: %.f90
	$(FC) $(FCFLAGS) -module $(SRC_DIR) -c $^ -o $@

$(OBJS): %.o: %.f90
	$(FC) $(FCFLAGS) -c $^ -o $@ 

debug:
	@echo "SRCS=$(SRCS)"
	@echo "OBJS=$(OBJS)"
	@echo "MODS=$(MODS)"
	@echo "MODOBJS=$(MODOBJS)"
	@echo "EXE=$(EXE)"
	@echo "FC =$(FC)"
	@echo "FLFLAGS=$(FLFLAGS)"
	@echo "FCFLAGS=$(FCFLAGS)"

.PHONY: clean
clean:
	rm -f $(SRC_DIR)*.o $(SRC_DIR)*.mod $(EXE) $(EXAMPLE_DIR)P_vec_*
