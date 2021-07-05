cmake_minimum_required(VERSION 3.16)
project(MFE1 LANGUAGES Fortran)

add_executable(go mfe_extents.f90 prob_procs.f90)

target_sources(go PUBLIC
  source/bc.f90
  source/bls_solver.f90
  source/common_io.f90
  source/elt_array.f90
  source/elt_lapl.f90
  source/elt_mfe.f90
  source/fpa.f90
  source/init_procs.f90
  source/mfe1.f90
  source/mfe_data.f90
  #source/mfe_extents.f90  example
  source/mfe_procs.f90
  source/mfe_solver.f90
  source/norm_data.f90
  source/precision.f90
  #source/prob_procs.f90  example
  source/system_io.f90
)
