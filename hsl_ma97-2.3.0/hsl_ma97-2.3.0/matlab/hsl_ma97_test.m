function hsl_ma97_test()
%
% Unit tests for hsl_ma97 matlab interface
%
fails = 0;

fprintf('Testing Poisson(2) real:\n')
A = gallery('poisson', 2);
fails = fails + test_with_matrix(A);

fprintf('Testing toy example real:\n')
A = sparse ([1 1 1 2 2 3 3 3 4 4], [2 3 4  1 3  1 2 3  1 4], [1.1 2.2 3.3, 1.1 4.4, 2.2 4.4 5.5, 3.3 6.6]);
fails = fails + test_with_matrix(A);

fprintf('Testing toy singular example real:\n')
A = sparse ([1 2 3 4 1 2 1 1], [1 1 1 1 2 2 3 4], [1 2 3 4 2 5 3 4]);
fails = fails + test_with_matrix(A);

fprintf('Testing toy example complex Hermitian:\n')
A = sparse ([1 1 1 2 2 3 3 3 4 4], [2 3 4  1 3  1 2 3  1 4], [1.1-i 2.2-i 3.3-i, 1.1+i 4.4-i, 2.2+i 4.4+i 5.5, 3.3+i 6.6]);
fails = fails + test_with_matrix(A);

fprintf('Testing toy singular example complex Hemitian:\n')
A = sparse ([1 2 3 4 1 2 1 1], [1 1 1 1 2 2 3 4], [1 2-i 3-i 4-i 2+i 5 3+i 4+i]);
fails = fails + test_with_matrix(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(fails == 0)
   fprintf('Test OK.\n')
else
   fprintf('Failed %i tests.\n', fails)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fails = test_with_matrix(A)
% Run through all tests with a specified matrix
fails = 0;

control.nemin = 8;
if(isreal(A))
   x = rand(size(A,1),1);
   x2 = rand(size(A,1));
else
   x = rand(size(A,1),1) + rand(size(A,1),1)*i;
   x2 = rand(size(A,1)) + rand(size(A,1))*i;
   if(A(1,2) == conj(A(2,1)))
      control.hermitian = true;
   end
end
b = A*x;
b2 = A*x2;

fprintf('   - mc68 ordering, separate calls\n')
[handleA, info] = hsl_ma97_factor(A, control);
soln = hsl_ma97_solve(handleA, b);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
hsl_ma97_destroy(handleA);

fprintf('   - symamd ordering, separate calls\n')
handleA = hsl_ma97_factor(A, control, symamd(A));
soln = hsl_ma97_solve(handleA, b);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
hsl_ma97_destroy(handleA);

fprintf('   - symamd ordering, backslash with handle\n')
[soln, info, handleA] = hsl_ma97_backslash(A, b, control, symamd(A));
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end

fprintf('   - subsequent solve\n')
[soln, info] = hsl_ma97_solve(handleA, b2, control);
res = norm(A*soln - b2, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b2, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end

fprintf('   - Fredholm solve')
b3 = b;
b3(4) = b3(4) - 1; % render inconsistent in singular case
B = [b b3];
[soln, info, consist] = hsl_ma97_solve(handleA, B, control);
res(1) = norm(A*soln(:,1) - B(:,1), inf) / ( norm(A, inf)*norm(soln(:,1), inf) + norm(B(:,1), inf) );
if(res(1) > 1e-14 && consist(1))
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
res(2) = norm(A*soln(:,2) - B(:,2), inf) / ( norm(A, inf)*norm(soln(:,2), inf) + norm(B(:,2), inf) );
if(res(2) > 1e-14 && consist(2))
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
fprintf(' consist = %d %d\n', consist(1), consist(2));
if(~consist(2))
   res = A*soln(:,2+size(B,2));
   if(norm(res, inf) > 1e-14)
      fprintf('fail Fredholm Ax = %d\n', res)
      fails = fails + 1;
   end
   res = soln(:,2+size(B,2))'*B(:,2);
   if(norm(res, inf) < 1e-14)
      fprintf('fail Fredholm x''*b = %d\n', res)
      fails = fails + 1;
   end
end

fprintf('   - L and L^T multiplies\n');
I = eye(size(A,1));
L = hsl_ma97_multiply(handleA, false, false, I);
LT = hsl_ma97_multiply(handleA, true, false, I);
res = norm(L'-LT, inf) / norm(A, inf);
if(res > 1e-14)
   fprintf('fail residual = %d\n', res);
   fails = fails + 1;
end

fprintf('   - L^-1 and L^-T multiplies\n');
I = eye(size(A,1));
L = hsl_ma97_multiply(handleA, false, true, I);
LT = hsl_ma97_multiply(handleA, true, true, I);
res = norm(L'-LT, inf) / norm(A, inf)^0.5;
if(res > 1e-14)
   fprintf('fail residual = %d\n', res);
   fails = fails + 1;
end

fprintf('   - sparse L^-1 and L^-T multiplies\n');
I = sparse(eye(size(A,1)));
for j = 1:size(A,1)
   L(:,j) = hsl_ma97_multiply(handleA, false, true, I(:,j));
end
LT = hsl_ma97_multiply(handleA, true, true, I);
res = norm(L'-LT, inf) / norm(A, inf)^0.5;
if(res > 1e-14)
   fprintf('fail residual = %d\n', res);
   fails = fails + 1;
end

%%% cleanup
hsl_ma97_destroy(handleA);

fprintf('   - simple backslash\n')
soln = hsl_ma97_backslash(A, b, control);
res = norm(A*soln - b, inf) / ( norm(A, inf)*norm(soln, inf) + norm(b, inf) );
if(res > 1e-14)
   fprintf('fail residual = %d\n', res)
   fails = fails + 1;
end
