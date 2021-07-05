cmake_minimum_required(VERSION 3.16)
project(MFE1 LANGUAGES Fortran)

add_executable(go)

target_sources(go PRIVATE
  precision.f90
  mfe_extents.f90
  mfe_data.f90
  norm_data.f90
  bc.f90
  system_io.f90
  common_io.f90
  elt_array.f90
  elt_lapl.f90
  prob_procs.f90
  elt_mfe.f90
  bls_solver.f90
  mfe_procs.f90
  fpa.f90
  mfe_solver.f90
  init_procs.f90
  mfe1.f90
)
