function hsl_ma97_expert(handle)
% HSL_MA97_DESTROY  Free memory associated with factorization.
%     hsl_ma97_destroy(handle) will free all memory associated with handle.
%     The factorization may not be reused again.
%
%     Usage: hsl_ma97_destroy(handle)
%
%     Please cite HSL as:
%     [1] HSL, a collection of Fortran codes for large-scale scientific
%         computation. See http://www.hsl.rl.ac.uk/.
%
%     This code is described in
%     [2] HSL_MA97: a bit-compatible multifrontal code for sparse symmetric
%         systems. J.D. Hogg and J.A. Scott. Technical Report RAL-TR-2011-024.
%
%     See also: ma97_backslash, ma97_factor, ma97_multiply, ma97_solve

hsl_ma97_expert('destroy', handle)
