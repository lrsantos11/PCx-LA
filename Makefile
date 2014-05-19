# Makefile.  Generated from Makefile.in by configure.

export PCx_ARCH = darwin
export WSSMP_LIB = 
export NG_LIB = ../Ng-Peyton/cholesky.a
export TARGETDIR = .

export CC ?= icc
export CFLAGS ?= -g -O2
export FC ?= @FC@
export FCFLAGS ?= @FCFLAGS@

# default
PCx:: NgPeyton
	@echo "***************************************************";
	@echo "Building PCx with the Ng Peyton Cholesky solver....";
	@echo "***************************************************";
	$(MAKE) $(MFLAGS) -C SRC PCx_NgPeyton

NgPeyton:: F2C
	@echo "***************************************************";
	@echo "Building the Ng Peyton Cholesky solver....";
	@echo "***************************************************";
	$(MAKE) $(MFLAGS) -C Ng-Peyton

F2C::
	@echo "***************************************************";
	@echo "Building the F2C library ....";
	@echo "***************************************************";
	$(MAKE) $(MFLAGS) -C F2C

PCx_wssmp::
	@echo "***************************************************";
	@echo "Building PCx with the IBM WSSMP Cholesky solver....";
	@echo "***************************************************";
	$(MAKE) $(MFLAGS) -C SRC PCx_wssmp

PCx_mysolver::
	@echo "***********************************************";
	@echo "Building PCx with the user-supplied solver     ";
	@echo "Requires file ./SRC/mysolver.c and             ";
	@echo "    library ./mysolver/libmysolver.a           ";
	@echo "***********************************************";
	$(MAKE) $(MFLAGS) -C SRC PCx_mysolver

mex:: NgPeyton
	@echo "***********************************************";
	@echo "The Matlab Interface is only available with    ";
	@echo "    the Ng-Peyton solver.                      ";
	@echo "***********************************************";
	$(MAKE) $(MFLAGS) -C mex

doc::
	-mkdir DOC/_doxygen
	cd SRC; doxygen
	$(MAKE) $(MFLAGS) -C DOC html

clean::
	rm -f PCx
	rm -f *.log
	rm -f mps/*.log
	rm -f mps/*.out	
	rm -f Ng-Peyton/*.[ao]
	rm -f SRC/*.[ao]
	rm -f mex/*.mex*
	rm -f mex/*.[ao]
	rm -f F2C/*.[ao]
