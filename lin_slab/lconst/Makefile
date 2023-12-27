FFC = mpinfort

FFLAGS  = -report-all -O3 -fpp $(EXTRA) #-mparallel
FFLAGS  += -fdiag-vector=2 -fdiag-inline=2 -fdiag-parallel=2

LIB = -L/system/apps/rhel7/sx/fftw/3.3.8/lib -lfftw3 -lm
INC = -I/system/apps/rhel7/sx/fftw/3.3.8/include
SRC = ../src/
MYL = ../src/
CWD = ./
OUT = rmhdfs.exe

MATH = rmhdfs_r0.1_math
# MATRIX = LB_m2.0_matrix
MATRIX = LB_m2.1_matrix

rmhdfs:	$(SRC)rmhdfs_r0.2_header.f90\
	$(SRC)rmhdfs_r0.1_mpienv.f90\
	$(MYL)$(MATH).f90\
	$(MYL)$(MATRIX).f90\
	$(SRC)rmhdfs_r0.1_clock.f90\
	$(SRC)rmhdfs_r0.1_intgrl.f90\
	$(SRC)rmhdfs_r0.1_bndry.f90\
	$(SRC)rmhdfs_r0.2_pssn_fftw.f90\
	$(SRC)rmhdfs_r0.2_fld.f90\
	$(SRC)rmhdfs_r0.1_feedback.f90\
	$(SRC)rmhdfs_r0.2_set.f90\
	$(SRC)rmhdfs_r0.2_diag.f90\
	$(SRC)rmhdfs_r0.2_advnc.f90\
	$(SRC)rmhdfs_r0.1_main.f90
	
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.2_header.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.1_mpienv.f90
	$(FFC) $(FFLAGS) -c $(MYL)$(MATH).f90 $(INC)
	$(FFC) $(FFLAGS) -c $(MYL)$(MATRIX).f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.1_clock.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.1_intgrl.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.1_bndry.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.2_pssn_fftw.f90 $(INC)
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.2_fld.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.1_feedback.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.2_set.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.2_diag.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.2_advnc.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r0.1_main.f90
	
	$(FFC) $(FFLAGS) rmhdfs_r0.1_main.o\
			rmhdfs_r0.2_header.o\
			rmhdfs_r0.1_mpienv.o\
			$(MATH).o\
			$(MATRIX).o\
			rmhdfs_r0.1_clock.o\
			rmhdfs_r0.1_intgrl.o\
			rmhdfs_r0.1_bndry.o\
			rmhdfs_r0.2_fld.o\
			rmhdfs_r0.1_feedback.o\
			rmhdfs_r0.2_set.o\
			rmhdfs_r0.2_diag.o\
			rmhdfs_r0.2_pssn_fftw.o\
			rmhdfs_r0.2_advnc.o\
			-o $(OUT) $(LIB)

clean:
	rm -f ./*.o ./*.mod ./*.lst ./*.exe ./*.err ./*.out ./sub.q.00* ./rmhdfs_r0.1.nml.00*

clear:
	rm -f ./*.o ./*.mod ./*.lst ./*.err ./*.out
