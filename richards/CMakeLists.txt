cmake_minimum_required(VERSION 3.16)
project(MFE1 LANGUAGES Fortran)

add_executable(go problem.f90)

target_sources(go PUBLIC
  source/kindpar.f90
  source/sysio.f90
  source/mfeconst.f90
  source/comio.f90
  source/mfedata.f90
  source/mfetypes.f90
  source/locarray.f90
  source/loclapl.f90
  source/locmfe.f90
  source/locglb.f90
  source/locmfeprocs.f90
  source/blssolve.f90
  source/bc.f90
  source/mfeprocs.f90
  source/norm.f90
  source/fpa.f90
  source/mfesolve.f90
  source/init.f90
  source/output.f90
  source/mfe1.f90
)
