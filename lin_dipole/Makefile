FFC = mpinfort

FFLAGS = -warn all -check bounds -check uninit -qopenmp -traceback -check all
FFLAGS  = -report-all -O3 -fpp $(EXTRA) #-mparallel
FFLAGS  += -fdiag-vector=2 -fdiag-inline=2 -fdiag-parallel=2

# FFLAGS = -xHOST -O3 -ipo -no-prec-div -fp-model fast=2 # Optimization
# FFLAGS += -qopenmp # OpenMP

# FFLAGS = -fw -Kfast,parallel -Kocl -Qt -Cpp -Haefosux

LIB = -L/system/apps/rhel7/sx/fftw/3.3.8/lib -lfftw3 -lm
#LIB = -L/usr/local/lib -lfftw3
# LIB = -lfftw3 -lm
# LIB = -SSL2 -lfftw3 -lm -SSL2BLAMP

#INC = -I/usr/local/Cellar/fftw/3.3.8_1/include
INC = -I/system/apps/rhel7/sx/fftw/3.3.8/include

SRC = ../src/
MYL = ../src/
CWD = ./
OUT = rmhdfs.exe

MATH = rmhdfs_r1.0_math
# MATRIX = LB_m2.0_matrix
MATRIX = LB_m2.1_matrix

rmhdfs:	$(SRC)rmhdfs_r1.0_header.f90\
	$(SRC)rmhdfs_r1.0_mpienv.f90\
	$(MYL)$(MATH).f90\
	$(MYL)$(MATRIX).f90\
	$(SRC)rmhdfs_r1.0_clock.f90\
	$(SRC)rmhdfs_r1.0_intgrl.f90\
	$(SRC)rmhdfs_r1.0_bndry.f90\
	$(SRC)rmhdfs_r1.0_pssn_fftw.f90\
	$(SRC)rmhdfs_r1.0_fld.f90\
	$(SRC)rmhdfs_r1.0_feedback.f90\
	$(SRC)rmhdfs_r1.0_diag.f90\
	$(SRC)rmhdfs_r1.0_set.f90\
	$(SRC)rmhdfs_r1.0_advnc.f90\
	$(SRC)rmhdfs_r1.0_main.f90
	
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_header.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_mpienv.f90
	$(FFC) $(FFLAGS) -c $(MYL)$(MATH).f90 $(INC)
	$(FFC) $(FFLAGS) -c $(MYL)$(MATRIX).f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_clock.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_intgrl.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_bndry.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_pssn_fftw.f90 $(INC)
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_fld.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_feedback.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_diag.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_set.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_advnc.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhdfs_r1.0_main.f90
	
	$(FFC) $(FFLAGS) rmhdfs_r1.0_main.o\
			rmhdfs_r1.0_header.o\
			rmhdfs_r1.0_mpienv.o\
			$(MATH).o\
			$(MATRIX).o\
			rmhdfs_r1.0_clock.o\
			rmhdfs_r1.0_intgrl.o\
			rmhdfs_r1.0_bndry.o\
			rmhdfs_r1.0_fld.o\
			rmhdfs_r1.0_feedback.o\
			rmhdfs_r1.0_diag.o\
			rmhdfs_r1.0_set.o\
			rmhdfs_r1.0_pssn_fftw.o\
			rmhdfs_r1.0_advnc.o\
			-o $(OUT) $(LIB)

clean:
	rm -f ./*.o ./*.mod ./*.lst ./*.exe ./*.err ./*.out ./*.L ./sub.q.00* ./rmhdfs_r1.0.nml.*

clear:
	rm -f ./*.o ./*.mod ./*.lst ./*.err ./*.out
