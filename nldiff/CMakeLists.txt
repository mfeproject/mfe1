cmake_minimum_required(VERSION 3.16)
project(MFE1 LANGUAGES Fortran)

add_executable(go init.f90 output.f90 problem.f90)

target_sources(go PUBLIC
  src/kindpar.f90
  src/sysio.f90
  src/mfeconst.f90
  src/comio.f90
  src/mfedata.f90
  src/mfetypes.f90
  src/locarray.f90
  src/loclapl.f90
  src/locmfe.f90
  src/locglb.f90
  src/locmfeprocs.f90
  src/blssolve.f90
  src/bc.f90
  src/mfeprocs.f90
  src/norm.f90
  src/fpa.f90
  src/mfesolve.f90
  #src/init.f90
  #src/output.f90
  src/mfe1.f90
)
