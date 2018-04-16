#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@cea.fr
#############################################################################
# Compilation time parameters
NDIM = 3
NVAR = 22 #d, u, v, w, p, (z,) xHii, xHeii, xHeiii + 4*npac = 5 (+ 1) + nIons + 4*nPacs
NPACS= 3
NPRE = 8
SOLVER = hydro #hydro or mhd
EXEC = ramses
UNITS= units.o

#--NORMAL RT
PATCH = ../jpatch/rt ../jpatch/proper_time

#--NORMAL RT with simple model for stellar UV feedback:
#PATCH = ../jpatch/rt_star_simple ../jpatch/rt ../jpatch/proper_time 

# Radiative Transfer objects------------------------------------------------
RTPARAMS = rt_parameters.o coolrates_module.o rt_spectra.o \
	   rt_cont_sources.o rt_read_params.o
RTOBJ = rt_flux_module.o rt_godunov_fine.o rt_init.o      \
        rt_init_xion.o rt_step.o rt_interpol_hydro.o

#############################################################################
# --- MPI, ifort syntax ------------------------------
F90 = mpif90
FFLAGS = -cpp -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DNVAR=$(NVAR) -DNPACS=$(NPACS)

# --- MPI, ifort syntax, additional checks -----------
#F90 = mpif90
#FFLAGS = -CB -cpp -fpe0 -ftrapuv -g -O3          \
#         -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DNOSYSTEM -DNVAR=$(NVAR) -DNPACS=$(NPACS)

#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI = -lmpi
#LIBMPI = -lfmpi -lmpi -lelan
LIBS = $(LIBMPI)
#############################################################################
# Sources directories are searched in this exact order
VPATH = $(PATCH):../$(SOLVER):../hydro:../pm:../poisson:../amr
#############################################################################
# All objects
MODOBJ = amr_parameters.o amr_commons.o random.o pm_parameters.o pm_commons.o poisson_parameters.o poisson_commons.o hydro_parameters.o hydro_commons.o  $(RTPARAMS) cooling_module.o bisection.o
AMROBJ = read_params.o init_amr.o init_time.o init_refine.o adaptive_loop.o amr_step.o update_time.o output_amr.o flag_utils.o physical_boundaries.o virtual_boundaries.o refine_utils.o nbors_utils.o hilbert.o load_balance.o title.o sort.o cooling_fine.o $(UNITS) light_cone.o
# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o move_fine.o newdt_fine.o particle_tree.o add_list.o remove_list.o star_formation.o sink_particle.o feedback.o
# Poisson solver objects
POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o multigrid_coarse.o multigrid_fine_commons.o multigrid_fine_fine.o multigrid_fine_coarse.o gravana.o boundary_potential.o rho_ana.o output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o courant_fine.o godunov_fine.o uplmde.o umuscl.o interpol_hydro.o godunov_utils.o condinit.o hydro_flag.o hydro_boundary.o boundana.o read_hydro_params.o synchro_hydro_fine.o 
# All objects
AMRLIB = $(MODOBJ) $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ) $(RTOBJ)
#############################################################################
ramses:	$(AMRLIB) ramses.o
	$(F90) $(FFLAGS) $(AMRLIB) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
#############################################################################
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
#############################################################################
clean :
	rm *.o *.$(MOD)
#############################################################################
