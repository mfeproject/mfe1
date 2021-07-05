
module local_mfe_procs

  use local_global, only: gather_local_solution, free_local_arrays, assemble_vector, assemble_matrix, assemble_diagonal
  use local_mfe, only: preprocessor, res_mass_matrix, eval_mass_matrix, reg_rhs, eval_dfdy
  use problem_pde, only: pde_rhs
  public

end module local_mfe_procs
