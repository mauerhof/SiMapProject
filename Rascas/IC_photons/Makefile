OPTIONS = -cpp

F90         = gfortran -ffree-line-length-none -ffree-form
FFLAGS      = -g -fcheck=all -ftrapv -fbacktrace -fbounds-check
FFLAGS += $(OPTIONS)



VPATH = :/home/cral/mauerhofer/rascas/f90/:./

OBJS = module_constants.o module_random.o module_utils.o


%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@


My_point_source: $(OBJS) My_point_source.o
	$(F90) $(OBJS) My_point_source.o -o My_point_source

clean:
	rm -f $(OBJS) *.o *.mod
