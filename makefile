# YouJunHu's makefile
#include ${PETSC_DIR}/lib/petsc/conf/variables
program_name=TEK
#main_src=test.f90
#main_src=t.f90
#main_src=main_adiabatic_electrons.f90
main_src=main.f90
#main_src=test_ion_trajectory.f90
#main_src=test_electron_trajectory.f90
#main_src=test_electron_trajectory2.f90
#main_src=test_electron_trajectory3.f90
#main_src=test_gc_trajectory.f90
#main_src=test_electron_deposition.f90
#main_src=test_sine_filter.f90
#field_solver_src=guess_em_field.f90
#field_solver_src= field_solver_sp_fd.f90
#field_solver_src= petsc_field_solver.F90
#field_solver_src=field_solver_flux_tube.f90
#source_term_src=prepare_source_terms_for_simple_iteration_scheme.f90
#source_term_src=prepare_source_terms_matrix_inversion.f90

ifeq  (${which_machine},edison)
  COMPILER=	ftn
   #OPTION = -FR -r8 -heap-arrays -O2 -g -traceback -check bounds
#OPTION = -FR -r8 -heap-arrays -O2 -g -traceback 
fftw_include=/global/homes/y/yjhu/installed/fftw-3.3.7/include
fftw_lib=-L /global/homes/y/yjhu/installed/fftw-3.3.7/lib/
else ifeq (${which_machine},cori)
  COMPILER=	ftn
  #OPTION = -FR -r8 -heap-arrays -O2 -g -traceback -check bounds -check all
#OPTION = -FR -r8 -heap-arrays -O2 -g -traceback -qopenmp
fftw_include=/global/homes/y/yjhu/installed/fftw-3.3.7/include
fftw_lib=-L /global/homes/y/yjhu/installed/fftw-3.3.7/lib/
#petsc_included:=-I/global/u2/y/yjhu/installed/petsc-3.9.3_cori/include -I/global/u2/y/yjhu/installed/petsc-3.9.3_cori/cori_yj/include
else ifeq (${which_machine},tianhe)
   COMPILER=	mpiifort
   #OPTION= -g -O0  -traceback -CB
   OPTION= -O3
   fftw_include=/THL8/software/fftw/3.3.4-icc16-IMPI5.1/include/
   fftw_lib=/THL8/software/fftw/3.3.4-icc16-IMPI5.1/lib/libfftw3.a
  lapack_location:=/THL6/software/lapack/3.8.0-icc16/lib64/liblapack.a
  blas_location:=  /THL6/software/lapack/3.8.0-icc16/lib64/libblas.a

else ifeq (${which_machine},3m)
     COMPILER=	mpif90
    fftw_include=/thfs3/software/fftw/3.3.8-gcc9.3.0-mpi-x/include/
    fftw_lib=/thfs3/software/fftw/3.3.8-gcc9.3.0-mpi-x/lib/libfftw3.a
  lapack_location:= /thfs3/software/lapack/3.8.0-gcc4.9.3/lib/liblapack.a
  blas_location:=  /thfs3/software/lapack/3.8.0-gcc4.9.3/lib/libblas.a
     OPTION=      -fmax-errors=1 -O3 -fimplicit-none -fbounds-check -fbacktrace  -Wmaybe-uninitialized -fno-realloc-lhs  -pedantic

else ifeq (${which_machine},hfc)
   COMPILER=	mpiifort
   OPTION= -O3 -lfftw3 
fftw_include=/usr/include
lapack_location:= /public/software/mathlib/lapack/intel/3.10.0/lib/liblapack.a
  blas_location:=  /public/software/mathlib/lapack/intel/3.10.0/lib/libblas.a

else
  COMPILER= mpif90
fftw_include=/usr/include
#  lapack_location:=/usr/lib/lapack/liblapack.so.3
#  blas_location:=/usr/lib/libblas/libblas.so.3
  lapack_location:=/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3
  blas_location:=  /usr/lib/x86_64-linux-gnu/blas/libblas.so.3

    fftw_lib:=/usr/lib/x86_64-linux-gnu/libfftw3.so.3
OPTION=    -g  -fmax-errors=1 -O0 -fimplicit-none  -Wsurprising -fbounds-check -fbacktrace  -Wno-unused-variable -Wmaybe-uninitialized -fno-realloc-lhs  -pedantic -Wconversion
#COMPILER=  ${FLINKER}
  #OPTION= -g -fbounds-check -fopenmp -fimplicit-none  -Wsurprising -Wall  -Wno-unused-variable
  #OPTION= -g -Wall -Wextra -Wconversion -pedantic  -ffpe-trap=zero,overflow,underflow -fbounds-check 
  #OPTION=-Wall -Wextra -Wimplicit-interface  -fmax-errors=1 -g -fcheck=all -fbacktrace -ffree-line-length-none 
  #OPTION=   -Og  -fmax-errors=1 -g -fimplicit-none -fbounds-check -fbacktrace -Wall -Wextra  -pedantic -Wno-unused -Wno-unused-dummy-argument -Warray-temporaries
#OPTION=   -g  -fmax-errors=1 -g -fimplicit-none -fbounds-check -fbacktrace  -pedantic -Wno-unused -Wno-unused-dummy-argument -Warray-temporaries -fcheck=all
#OPTION=     -fmax-errors=1 -g -fimplicit-none -fbounds-check -fbacktrace  -pedantic -Wno-unused -Wno-unused-dummy-argument  -fcheck=all
#OPTION=      -fmax-errors=1 -Og -fimplicit-none -fbounds-check -fbacktrace  -Wmaybe-uninitialized  -Wimplicit-interface


#-Wunused-parameter or -Wno-unused-parameter, -Wall
#-finit-integer=214483647 -finit-real=nan  -ffpe-trap=zero,overflow,underflow
#Here the option "-finit-integer=214483647 -finit-real=nan" is to initialize stack local variables to an unusual value to aid error detection.
#the option "-ffpe-trap=zero,overflow,underflow" tells Fortran to trap the listed floating-point error (fpe). Having zero on the list means that if you divid by zero the code will die rather than setting the result to +Infinity and continuing. Similarly, if overflow is on the list it will halt if you try to store a number larger than can be stored for the type of real number you are using because the exponent is too large.

#OPTION= -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
#OPTION=  -O3 -ftree-vectorizer-verbose=1 -fopenmp -fimplicit-none

#petsc_included:=-I/home/yj/install_source/petsc-3.9.3/include -I/home/yj/install_source/petsc-3.9.3/arch-linux2-c-debug/include -I/home/yj/installed/mpich-3.2/include
endif


COMPILE=	$(COMPILER) $(OPTION) -c 

f90source = modules.f90 pputil_yj.f90  splines.f90  math.f90  \
                equilibrium.f90 magnetic_coordinates.f90  table_in_mc.f90 transform.f90 fk_module.f90 \
               profiles.f90 gk_module.f90     \
               metric.f90  mapping.f90  load_gk.f90 \
                 field_lines_analyse.f90  \
                 deposit_gk.f90  deposit_fk.f90 \
              communication_connection.f90 filter.f90 derivatives_in_xyz.f90   \
             gk_polarization.f90 gyro_ring.f90 gyro_average.f90 poisson.f90 ampere.f90 \
            force.f90          boris.f90   push_fk_orbit.f90  \
              drift.f90 push_gk_orbit.f90 \
            push_gk_weight.f90   initial_weight.f90       \
            push_fk_weight.f90    $(source_term_src)     \
            diagnosis.f90  write_data_for_restarting.f90      $(field_solver_src) $(main_src)


#f90objs_tmp= $(f90source:.f90=.o)
#f90objs= $(f90objs_tmp:.F90=.o)
f90objs= $(f90source:.f90=.o)
#f77source=  pnpoly.for
#f77objs= $(f77source:.for=.o)
build_directory=build/

main_obj=$(main_src:.f90=.o)
$(program_name):    $(f90objs)  $(f77objs)   makefile
	 $(COMPILER)  $(OPTION)   $(f90objs) $(f77objs) $(lapack_location) $(blas_location) $(fftw_lib)  -o $(program_name)
%.o: %.f90 modules.f90
	$(COMPILE) -I$(fftw_include) $< -o $@

.PHONY : clean run tarfile create_tags to_cori to_edison on_tianhe
clean :
	 rm -f $(program_name)   $(f90objs) $(f77objs) *.mod  $(program_name) ms/*

create_tags:
	 etags   $(f90source) $(f77source)

tarfile:
	mkdir $(program_name)_version`date +"%Y-%m-%d"` && cp -r $(f90source) $(f77source) makefile input.nmlt $(program_name)_version`date +"%Y-%m-%d"` && tar -cf $(program_name)_version`date +"%Y-%m-%d"`.tar $(program_name)_version`date +"%Y-%m-%d"` && rm -r $(program_name)_version`date +"%Y-%m-%d"` 
merge_to_one_file:
	cat $(f90source) >one_file.f90
	$(COMPILER)  $(OPTION) one_file.f90   $(lapack_location) $(blas_location) -I$(fftw_include) $(fftw_lib) -lfftw3 -lm -o ./a.out
#sync_to_cluster:
#	 make clean && rsync -avz --delete ./ yj@202.127.204.22:/scratch/yj/TEK
to_cori:
	#rsync -avz --delete ./ yjhu@cori.nersc.gov:~/gk_ion_itg/
	rsync -avz --delete ./ yjhu@cori.nersc.gov:~/my_etg/
to_edison:
	rsync -avz --delete ./ yjhu@edison.nersc.gov:~/my_etg_edison/
run: $(program_name)
	date +%H:%M:%S
	#export OMP_NUM_THREADS=1 && mpiexec -n 4 ./$(program_name) -pc_type bjacobi -sub_pc_type ilu -ksp_type cg -ksp_monitor -ksp_max_it 100
	export OMP_NUM_THREADS=1 && mpiexec -n 4 ./$(program_name) 
	#export OMP_NUM_THREADS=1 && mpiexec -n 4 ./$(program_name)  -pc_type svd -pc_svd_monitor
	date +%H:%M:%S
on_tianhe:
	export which_machine=tianhe; make
to_tianhe:
	make clean && rsync -avz --delete --exclude=".git" ./ huyj@192.168.2.107:/THL8/home/huyj/tek/
to_3m:
	make clean && rsync -avz --delete ./ huyj@192.168.10.20:~/tek/
on_3m:
	export which_machine=3m; make

on_hfc:
	export which_machine=hfc; make

to_hfc:
	make clean && rsync -avz --delete --exclude=".git" -e "ssh -p 65062  -i ~/Downloads/hfeshell.nscc-hf.cn_rsa.txt" ./ hfcas_user26@hfeshell.nscc-hf.cn:~/tek/

