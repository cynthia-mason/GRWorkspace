function hsl_ma97_full_test

clear all;

A = gallery('poisson', 2);
x = [1; 2; 3; 4];
b = A*x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insufficient arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
   handle = hsl_ma97_factor();
   error('Unexpected success at insufficient inputs to hsl_ma97_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "A" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_factor')
   end
end

try
   X = hsl_ma97_backslash();
   error('Unexpected success at insufficient inputs to hsl_ma97_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "A" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_backslash')
   end
end

try
   X = hsl_ma97_backslash(A);
   error('Unexpected success at insufficient inputs to hsl_ma97_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "B" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_backslash')
   end
end

try
   X = hsl_ma97_solve();
   error('Unexpected success at insufficient inputs to hsl_ma97_solve')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "handle" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_solve')
   end
end

handle = hsl_ma97_factor(A);
try
   X = hsl_ma97_solve(handle);
   error('Unexpected success at insufficient inputs to hsl_ma97_solve')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "B" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_solve')
   end
end
hsl_ma97_destroy(handle);

try
   Y = hsl_ma97_multiply();
   error('Unexpected success at insufficient inputs to hsl_ma97_multiply')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "handle" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_multiply')
   end
end

handle = hsl_ma97_factor(A);
try
   Y = hsl_ma97_multiply(handle);
   error('Unexpected success at insufficient inputs to hsl_ma97_multiply')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "transpose" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_multiply')
   end
end
hsl_ma97_destroy(handle);

handle = hsl_ma97_factor(A);
try
   Y = hsl_ma97_multiply(handle, false);
   error('Unexpected success at insufficient inputs to hsl_ma97_multiply')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "inverse" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_multiply')
   end
end
hsl_ma97_destroy(handle);

handle = hsl_ma97_factor(A);
try
   Y = hsl_ma97_multiply(handle, false, false);
   error('Unexpected success at insufficient inputs to hsl_ma97_multiply')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "X" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_multiply')
   end
end
hsl_ma97_destroy(handle);

try
   hsl_ma97_destroy();
   error('Unexpected success at insufficient inputs to hsl_ma97_destroy')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Input argument "handle" is undefined');
   str2 = strtrim('Not enough input arguments.');
   if (size(strfind(str,str1),2)==0 && size(strfind(str,str2),2)==0)
      str
      error('Failure at insufficient inputs to hsl_ma97_destroy')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Too many input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

control.nemin=8;
P = symamd(A);
extra = 1;

try
   handle = hsl_ma97_factor(A, control, P, extra);
   error('Unexpected success at too many inputs to hsl_ma97_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many inputs to hsl_ma97_factor')
   end
end

handle = hsl_ma97_factor(A);
try
   X = hsl_ma97_solve(handle, b, control, extra);
   error('Unexpected success at too many inputs to hsl_ma97_solve')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many inputs to hsl_ma97_solve')
   end
end
hsl_ma97_destroy(handle);

handle = hsl_ma97_factor(A);
try
   X = hsl_ma97_multiply(handle, false, false, b, control, extra);
   error('Unexpected success at too many inputs to hsl_ma97_multiply')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many inputs to hsl_ma97_multiply')
   end
end
hsl_ma97_destroy(handle);

try
   X = hsl_ma97_backslash(A, b, control, P, extra);
   error('Unexpected success at too many inputs to hsl_ma97_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many inputs to hsl_ma97_backslash')
   end
end

try
   hsl_ma97_destroy(handle, extra);
   error('Unexpected success at too many inputs to hsl_ma97_destroy')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many input arguments.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many inputs to hsl_ma97_destroy')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Too many output arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
   [handle, info, extra] = hsl_ma97_factor(A, control, P);
   error('Unexpected success at too many outputs to hsl_ma97_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many output arguments.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many outputs to hsl_ma97_factor')
   end
end

handle = hsl_ma97_factor(A);
try
   [x, info, consist, extra] = hsl_ma97_solve(handle, b);
   error('Unexpected success at too many outputs to hsl_ma97_solve')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many output arguments.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many outputs to hsl_ma97_solve')
   end
end
hsl_ma97_destroy(handle);

handle = hsl_ma97_factor(A);
try
   [x, info, extra] = hsl_ma97_multiply(handle, false, false, b);
   error('Unexpected success at too many outputs to hsl_ma97_multiply')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many output arguments.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many outputs to hsl_ma97_multiply')
   end
end
hsl_ma97_destroy(handle);

try
   [x, info, handle, extra] = hsl_ma97_backslash(A, b);
   error('Unexpected success at too many outputs to hsl_ma97_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many output arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many outputs to hsl_ma97_backslash')
   end
end

try
   extra = hsl_ma97_destroy(handle);
   error('Unexpected success at too many outputs to hsl_ma97_destroy')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Too many output arguments');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at too many outputs to hsl_ma97_destroy')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Not a sparse matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Adense = full(A);

try
   [handle, info] = hsl_ma97_factor(Adense);
   error('Unexpected success at dense matrix to hsl_ma97_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Error in argument A. Expected sparse matrix.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at dense matrix to hsl_ma97_factor')
   end
end

try
   [handle, info] = hsl_ma97_backslash(Adense, b);
   error('Unexpected success at dense matrix to hsl_ma97_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Error in argument A. Expected sparse matrix.');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at dense matrix to hsl_ma97_backslash')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A not square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Arect = [A [1; 2; 3; 4;]];

try
   [handle, info] = hsl_ma97_factor(Arect);
   error('Unexpected success at rectangular matrix to hsl_ma97_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('The matrix must be square');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at rectangular matrix to hsl_ma97_factor')
   end
end

try
   x = hsl_ma97_backslash(Arect, b);
   error('Unexpected success at rectangular matrix to hsl_ma97_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('The matrix must be square');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at rectangular matrix to hsl_ma97_backslash')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A and b inconsistent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b5 = [1; 2; 3; 4; 5];

handle = hsl_ma97_factor(A);
try
   soln = hsl_ma97_solve(handle, b5);
   error('Unexpected success at inconsistent A, b to hsl_ma97_solve')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Dimensions of A and b inconsistent: A=    4x    4, b=    5x    1');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at inconsistent A, b to hsl_ma97_solve')
   end
end
hsl_ma97_destroy(handle);

handle = hsl_ma97_factor(A);
try
   soln = hsl_ma97_multiply(handle, false, false, b5);
   error('Unexpected success at inconsistent L, X to hsl_ma97_multiply')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Dimensions of L and X inconsistent: L=    4x    4, X=    5x    1');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at inconsistent L, X to hsl_ma97_multiply')
   end
end
hsl_ma97_destroy(handle);

try
   soln = hsl_ma97_backslash(A, b5);
   error('Unexpected success at inconsistent A, b to hsl_ma97_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Dimensions of A and b inconsistent: A=    4x    4, b=    5x    1');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at inconsistent A, b to hsl_ma97_backslash')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bad P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pbad = [1 2 2 3];

try
   handle = hsl_ma97_factor(A, control, Pbad);
   error('Unexpected success at bad P to hsl_ma97_factor')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Problem with P. P =    1   2   2   3');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at bad P to hsl_ma97_factor')
   end
end

try
   x = hsl_ma97_backslash(A, b, control, Pbad);
   error('Unexpected success at bad P to hsl_ma97_backslash')
catch
   errstr = lasterror();
   str = strtrim(errstr.message);
   str1 = strtrim('Problem with P. P =    1   2   2   3');
   if (size(strfind(str,str1),2)==0)
      str
      error('Failure at bad P to hsl_ma97_backslash')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('All tests succeeded.\n')
