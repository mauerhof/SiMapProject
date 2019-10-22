#############################################################################
# Fortran compiler options and directives

#SYSTYPE="CCF"
#SYSTYPE="lyocral"
#SYSTYPE="lyocralP"
SYSTYPE="macP"
#SYSTYPE="mac"
#SYSTYPE="ccage"

OPTIONS = -cpp


ifeq ($(SYSTYPE),"lyocral")
F90         = ifort
#FFLAGS      = -g -traceback -fpp -check all -debug -warn all -ftrapuv 
##FFLAGS      = -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS      = -O3 -fpp -ftz -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
##FFLAGS      = -p -g
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
endif

ifeq ($(SYSTYPE),"lyocralP")
F90         = mpif90 -ffree-line-length-none -ffree-form
FFLAGS      = #-g -fcheck=all -ftrapv -fbacktrace -Wall
FFLAGS      = -O2
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
endif


ifeq ($(SYSTYPE),"CCF")
F90         = ifort
#FFLAGS      = -g -fpp -check all -debug -warn all -ftrapuv -traceback
##FFLAGS      = -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS      = -O3 -fpp -ftz  -xavx -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
##FFLAGS      = -p -g
FFLAGS += $(OPTIONS)
LDFLAGS     =  -lmpi
endif

ifeq ($(SYSTYPE),"mac")
#F90         = gfortran-mp-4.7 -ffree-line-length-none -ffree-form  ###-fopenmp
F90         = gfortran -ffree-line-length-none -ffree-form -fopenmp
FFLAGS      = -g -fcheck=all -ftrapv -fbacktrace -fbounds-check
#FFLAGS      = -O2
FFLAGS += $(OPTIONS)
LDFLAGS     = 
endif


ifeq ($(SYSTYPE),"macP")
#F90         = mpif90-openmpi-mp
F90	= /usr/bin/mpif90 -ffree-line-length-none -frecord-marker=4 -O3 -g -fbacktrace -fopenmp -fbounds-check
#FFLAGS      = -g -fcheck=all -ftrapv -fbounds-check -fbacktrace 
FFLAGS      = -O3
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
endif


ifeq ($(SYSTYPE),"ccage")
F90         = mpif90
#FFLAGS      = -g -fcheck=all -ffree-line-length-none -ffree-form -ftrapv -fbounds-check -fbacktrace
FFLAGS      = -O3 -ffree-line-length-none -ffree-form
FFLAGS += $(OPTIONS)
LDFLAGS     = 
endif


#############################################################################
VPATH = :/home/cral/mauerhofer/rascas/f90/:/home/cral/mauerhofer/krome/KromeFiles/All_paper/:./

OBJS = module_parallel_mpi_mine.o opkda2.o opkda1.o opkdmain.o krome_user_commons.o krome_all.o module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_spectra.o module_ramses.o module_krome.o module_cell.o module_file.o module_worker_mine.o module_master_mine.o
#OBJS = module_parallel_mpi_mine.o opkda2.o opkda1.o opkdmain.o krome_user_commons.o krome_commons.o krome_constants.o krome_user_commons.o krome_fit.o krome_getphys.o krome_gadiab.o krome_grfuncs.o krome_phfuncs.o krome_subs.o krome_stars.o krome_dust.o krome_photo.o krome_tabs.o krome_coolingGH.o krome_cooling.o krome_heating.o krome_ode.o krome_user.o krome_reduction.o krome.o module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_spectra.o module_ramses.o module_krome.o module_cell.o module_file.o module_worker_mine.o module_master_mine.o
#OBJS = module_parallel_mpi_mine.o opkda2.o opkda1.o opkdmain.o krome_user_commons.o krome_commons.o krome_constants.o krome_user_commons.o krome_fit.o krome_getphys.o krome_gadiab.o krome_grfuncs.o krome_phfuncs.o krome_dust.o krome_photo.o krome_subs.o krome_stars.o krome_tabs.o krome_coolingGH.o krome_cooling.o krome_heating.o krome_ode.o krome_user.o krome_reduction.o krome.o module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_spectra.o module_ramses.o module_krome.o module_cell.o module_file.o module_worker_mine.o module_master_mine.o

#############################################################################
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
#rule for f
%.o:%.f
	$(F90) -ffree-line-length-none -w -c $^ -o $@

#############################################################################
ion_fractions: $(OBJS) ion_fractions.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) ion_fractions.o

clean:
	rm -f $(OBJS) *.o *.mod
#############################################################################
