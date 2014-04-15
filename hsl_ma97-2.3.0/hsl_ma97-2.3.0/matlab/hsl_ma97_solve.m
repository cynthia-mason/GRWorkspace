function [X, info, varargout] = hsl_ma97_solve(handle, B, varargin)
% HSL_MA97_SOLVE  Sparse Symmetric Indefinite Solve.
%     X = hsl_ma97_solve(handle, B) solves the equation AX=B for X given
%     precomputed factors associated with handle. The handle must have been
%     obtained by a prior call to hsl_ma97_factor of hsl_ma97_backslash.
%
%     Usage: X = hsl_ma97_solve(handle, B)
%            [X, info] = hsl_ma97_solve(handle, B, control)
%            [X, info, consist] = hsl_ma97_solve(handle, B, control)
%
%     The optional argument CONTROL may have the following components set. If
%     they are not set then the stated default is used.
%     control.num_threads  - Number of threads on which to run. Default is the
%                            maximum available.
%
%     The optional return value INFO will have some of the following components
%     set on exit.
%     info.solve_time         - Wall clock time for Fortran ma97_solve call
%
%     The presence of the optional return value CONSIST invokes the Fredholm
%     alternative. For each right-hand side i, consist(i) will be true if the
%     right-hand side is consistent, and false if it is inconsistent. For
%     each rhs X(:,i) will contain the same solution as before, however for each
%     inconsistent rhs X(:,size(B,2)+i) will contain a vector satisfying the
%     Fredholm alternative Ax=0 and x'*b!=0. Note that only singular systems
%     can have inconsistent rhs.
%
%     Please cite HSL as:
%     [1] HSL, a collection of Fortran codes for large-scale scientific
%         computation. See http://www.hsl.rl.ac.uk/.
%
%     This code is described in
%     [2] HSL_MA97: a bit-compatible multifrontal code for sparse symmetric
%         systems. J.D. Hogg and J.A. Scott. Technical Report RAL-TR-2011-024.
%
%     See also: ma97_backslash, ma97_destroy, ma97_factor, ma97_multiply

optargin = size(varargin,2);
if(nargout > 3)
   error 'Too many output arguments.';
end
if(optargin == 0)
   if(nargout <= 2)
      [X, info] = hsl_ma97_expert('solve', handle, B);
   else
      % Fortran can't interact with matlab logical as Mathworks don't define
      % how to, so we return int4 and convert on the matlab side.
      [X, info, varargout{1}] = hsl_ma97_expert('solve', handle, B);
      varargout{1} = logical(varargout{1});
   end
elseif(optargin == 1)
   if(nargout <= 2)
      [X, info] = hsl_ma97_expert('solve', handle, B, varargin{1});
   else
      % Fortran can't interact with matlab logical as Mathworks don't define
      % how to, so we return int4 and convert on the matlab side.
      [X, info, varargout{1}] = hsl_ma97_expert('solve', handle, B, varargin{1});
   end
else
   error ('Too many arguments')
end
