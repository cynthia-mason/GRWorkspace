function [Y, info] = hsl_ma97_multiply(handle, transpose, inverse, X, varargin)
% HSL_MA97_MULTIPLY  Multiply by factor of Sparse Symmetric Matrix
%     Y = hsl_ma97_multiply(handle, transpose, inverse, X) forms a product
%     of X with the lower triangular factor L depending on the values of
%     transpose and inverse:
%
%             transpose= false     true
%     inverse=false        Lx      L^Tx
%     inverse=true      L^{-1}x  L^{-T}x
%
%     Usage: Y = hsl_ma97_multiply(handle, transpose, inverse, X)
%            [Y, info] =  hsl_ma97_multiply(handle, transpose, inverse, X,
%                                           control)
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
%     Please cite HSL as:
%     [1] HSL, a collection of Fortran codes for large-scale scientific
%         computation. See http://www.hsl.rl.ac.uk/.
%
%     This code is described in
%     [2] HSL_MA97: a bit-compatible multifrontal code for sparse symmetric
%         systems. J.D. Hogg and J.A. Scott. Technical Report RAL-TR-2011-024.
%
%     See also: ma97_backslash, ma97_destroy, ma97_factor, ma97_solve

optargin = size(varargin, 2);
if(optargin == 0)
   % Note: +0 forces logical=>int conversion
   [Y, info] = hsl_ma97_expert('multiply', handle, transpose+0, inverse+0, X);
elseif(optargin == 1)
   % Note: +0 forces logical=>int conversion
   [Y, info] = hsl_ma97_expert('multiply', handle, transpose+0, inverse+0, X, varargin{1});
else
   error ('Too many arguments')
end
