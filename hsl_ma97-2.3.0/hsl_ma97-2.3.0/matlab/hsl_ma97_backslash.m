function [X, varargout] = ma97_backslash(A, B, varargin)
% HSL_MA97_BACKSLASH Sparse Symmetric Indefinite Linear Solve.
%     X = hsl_ma97_backslash(A, B) solves the equation AX=B for X by means of
%     a symmetric indefinite factorization PAP' = LDL'.  If A is
%     positive-definite then a PAP'=LL' factorization may be computed instead
%     by setting control.pos_def=true. P is generated automatically to reduce
%     fill-in. A is assumed to be symmetric and only the lower triangular part
%     is referenced. A must be square.
%     Complex matrices are assumed to be symmetric unless they are specified as
%     Hermitian by setting the value control.hermitian=true.
%
%     Usage: X = hsl_ma97_backslash(A, B)
%            [X, info, handle] = ma97_backslash(A, B, control, P)
%
%     control is a structure described below. P is a permuaton such as that
%     output from symamd(A). info is a structure described below. handle
%     provides a reference to the factorization so it may be reused at a later
%     date.
% 
%     control may have the following components set. If they are not set then
%     the stated default is used.
%     control.hermitian    - True or false. Determines if a complex matrix is
%                            treated as Hermitian (true) or symmetric (false).
%                            Default is false.
%     control.nemin        - Maximum number of columns in candidates for
%                            supernode amalgamation. Default is 8.
%     control.num_threads  - Number of threads on which to run. Default is the
%                            maximum available.
%     control.ordering     - If an explicit P is not supplied, determines the
%                            ordering algorithm employed:
%                                   1 : Approximate minimum Degree
%                                   2 : Minimum Degree
%                                   3 : METIS (requires METIS, see README)
%                                   4 : MA47 ordering
%                                   5 : Heuristic choice between AMD and METIS
%                                       for parallel computation.
%                                   6 : Heuristic choice between AMD and METIS
%                                       for serial computation.
%                                   7 : Matching based ordering using AMD.
%                                   8 : Matching based ordering using METIS.
%                            Default is 5.
%                            Note that for 7 and 8, MC64 scaling is calculated
%                            as a side-effect and can therefore be used for
%                            free.
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
%     control.u            - Relative pivot tolerance threshold. Default
%                            is 0.01.
%
%     On return, info will have the following components set.
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
%     info.factor_solve_time  - Wall clock time for Fortran ma97_factor_solve
%                               call
%
%     Note that if the optional return value handle is used, then memory will
%     be tied up storing the factorization. This may be recovered by calling
%     ma97_destroy.
%
%     Please cite HSL as:
%     [1] HSL, a collection of Fortran codes for large-scale scientific
%         computation. See http://www.hsl.rl.ac.uk/.
%
%     This code is described in:
%     [2] HSL_MA97: a bit-compatible multifrontal code for sparse symmetric
%         systems. J.D. Hogg and J.A. Scott. Technical Report RAL-TR-2011-024.
%
%     See also: ma97_destroy, ma97_factor, ma97_multiply, ma97_solve

optargin = size(varargin,2);
optargout = max(nargout,1)-1;
if(optargin>2)
   error ('Too many arguments')
end
if(optargout == 0)
   if(optargin == 0)
      X = hsl_ma97_expert('backslash', A, B);
   elseif(optargin == 1)
      X = hsl_ma97_expert('backslash', A, B, varargin{1});
   elseif(optargin == 2)
      X = hsl_ma97_expert('backslash', A, B, varargin{1}, varargin{2});
   end
elseif(optargout == 1)
   if(optargin == 0)
      [X, info] = hsl_ma97_expert('backslash', A, B);
   elseif(optargin == 1)
      [X, info] = hsl_ma97_expert('backslash', A, B, varargin{1});
   elseif(optargin == 2)
      [X, info] = ...
         hsl_ma97_expert('backslash', A, B, varargin{1}, varargin{2});
   end
   varargout(1) = {info};
elseif(optargout == 2)
   if(optargin == 0)
      [X, info, handle] = hsl_ma97_expert('backslash', A, B);
   elseif(optargin == 1)
      [X, info, handle] = ...
         hsl_ma97_expert('backslash', A, B, varargin{1});
   elseif(optargin == 2)
      [X, info, handle] = ...
         hsl_ma97_expert('backslash', A, B, varargin{1}, varargin{2});
   end
   varargout(1) = {info};
   varargout(2) = {handle};
else
    error ('Too many output arguments')
end
