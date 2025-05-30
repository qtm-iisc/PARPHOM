# Try to detect Intel Fortran (ifort/ifx), else fall back to gfortran
FC := $(shell (command -v ifort >/dev/null 2>&1 && echo ifort) || (command -v ifx >/dev/null 2>&1 && echo ifx) || echo gfortran)
F2PY := FC=$(FC) python -m numpy.f2py
PYTHON_EXT_SUFFIX := $(shell python3-config --extension-suffix)

all: pyphutil/bz_integration$(PYTHON_EXT_SUFFIX) pyphutil/neighbor$(PYTHON_EXT_SUFFIX)

pyphutil/bz_integration$(PYTHON_EXT_SUFFIX): pyphutil/fortran_routines/bz_integration.f90
	$(F2PY) -c $< -m bz_integration 

pyphutil/neighbor$(PYTHON_EXT_SUFFIX): pyphutil/fortran_routines/neighbor_list.f90
	$(F2PY) -c $< -m neighbor 

install: all
	pip install .

clean:
	rm -f pyphutil/bz_integration*.so pyphutil/neighbor*.so

.PHONY: all install clean
