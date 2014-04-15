function varargout = hsl_ma97_expert(varargin)
% HSL_MA97  Sparse Symmetric Indefinite Linear Solver.
%     hsl_ma97_expert is a direct Fortran mex interface. The first parameter
%     specifies the action to be performed. Inputs and outputs depend on the
%     action.
%
%     Complex matrices are assumed to be symmetric unless they are specified as
%     Hermitian by setting the value control.hermitian=true.
%
%     Real symmetric and complex Hermitian matrices are assumed to be indefinite
%     unless control.pos_def=true.
%
%     [handle, info] = hsl_ma97_expert('factor', A, control, P)
%        Performs the factorization of a symmetric matrix A and returns an
%        integer handle for the factorization. This is equivalent to calling
%        the Fortran routines ma97_analyse and ma97_factor.
%        The argument control is optional and is described below.
%        The argument P is optional and is a vector as returned by e.g.
%           symamd(A). If it is not present hsl_ma97 will find its own
%           fill-reducing permutation.
%        The argument info is optional and is described below.
%
%     [X, info, consist] = hsl_ma97_expert('solve', handle, B, control)
%        Solves the equation AX = B for X using a factorization previously
%        computed by either the 'factor' or 'backslash' actions.
%        The argument consist is optional and its presence invokes the Fredholm
%        alternative for singular systems. If rhs i is consistent, consist(i)
%        will be true and X(:,i) will contain a solution. If rhs i is
%        inconsistent consist(i) will be false and X(:,i) will satsify Ax=0 and
%        x'*b!=0.
%        The arguments control and info are optional and are described below.
%
%     [X, info, handle] = hsl_ma97_expert('backslash', A, B, control, P)
%        Combines the 'factor' and 'solve' actions in a more efficient fashion.
%        It solves the equation AX = B for X using a factorization it computes.
%        The arguments control and info are optional and are described below.
%        The argument P is optional and as described for the 'factor' option
%        above.
%        The argument handle is optional, and if supplied then may be used to
%        reference the factorization in future calls to 'solve'.
%
%     hsl_ma97_expert('destroy', handle)
%        Destroys the factorization referecenced by handle, freeing the memory
%        associated with it.
%
%     The optional argument CONTROL may have the following components set. If
%     they are not set then the stated default is used.
%     control.hermitian    - True or false. Determines if a complex matrix is
%                            treated as Hermitian (true) or symmetric (false).
%                            Default is false.
%     control.nemin        - Maximum number of columns in candidates for
%                            supernode amalgamation. Default is 32.
%     control.num_threads  - Number of threads on which to run. Default is the
%                            maximum available.
%     control.ordering     - If an explicit P is not supplied, determines the
%                            ordering algorithm employed:
%                                   1 : Approximate minimum Degree
%                                   2 : Minimum Degree
%                                   3 : METIS (requires METIS, see README)
%                                   4 : MA47 ordering
%                                   5 : Heurisitc choice between AMD and METIS
%                            Default is 5.
%     control.pos_def      - True or false. Determines if a matrix is treated
%                            as positive-definite (true) or indefinite (false).
%                            In the complex case only Hermitian matrices can
%                            be positive definite.
%                            Default is false.
%     control.scaling      - Determines if scaling is to be used with values:
%                                 <=0 : No scaling
%                                   1 : MC64 Bipartite matching based scaling
%                                   2 : MC77 (1 itr inf norm, 3 itr one norm)
%                                 >=4 : MC30
%                            Default is 0.
%     control.small        - Pivots of modulus less than this are treated as
%                            zero. Default is 1e-20.
%     control.u            - Initial relative pivot tolerance threshold. Default
%                            is 0.01.
%
%     The optional return value INFO will have some of the following components
%     set on exit.
%     info.matrix_rank        - Number of non-zero pivots.
%     info.num_delay          - Number of delayed pivots.
%     info.num_factor         - Number of entries in the factors (after
%                               supernode amalgamation and pivoting).
%     info.num_flops          - Number of floating point operations to form
%                               factors (after supernode amalgamation and
%                               pivoting).
%     info.num_neg            - Number of negative pivots in factors.
%     info.num_two            - Number of 2x2 pivots used in factorization.
%     info.order              - Ordering used. One of 'AMD', 'MD', 'MC47',
%                               'MeTiS' or 'user'.
%     info.analyse_time       - Wall clock time for Fortran ma97_analyse call
%     info.factor_time        - Wall clock time for Fortran ma97_factor call
%     info.factor_solve_time  - Wall clock time for Fortran ma97_factor_solve
%                               call
%     info.solve_time         - Wall clock time for Fortran ma97_solve call

error('hsl_ma97 must be compiled before use')
