function [M] = find_coeff_matrix(x,p,q)
%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2024. November 19. (2023a)
%

n = numel(p);
m = numel(q);

M = sym('m',[n,m]);
vars = M(:);

z = M*q - p;

c = sym([]);
for i = 1:n
    ci = coeffs(z(i),x);
    c = [c ; ci(:)];
end

[A,b] = equationsToMatrix(c,vars);

M = double(reshape(A\b,[n,m]));

end