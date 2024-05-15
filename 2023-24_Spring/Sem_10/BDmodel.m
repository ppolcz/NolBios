function [G,A] = BDmodel(nV,nE,args)
arguments
    nV,nE
    args.maxIter = 100;
end
%%
%   Bemenete:
%   nV       - csúcsok száma
%   nE       - az élek száma
%%

[i,j] = meshgrid(1:nV);
P = tril((i ./ j).^3,-1);

R = rand(nV,nV);

nE_ = @(s) sum(s*R(:) <= P(:));

s = [1,1];
while nE_(s(2)) > nE
    s(2) = s(2) * 2;
end
while nE_(s(1)) < nE
    s(1) = s(1) / 2;
end

actnE = nE_(mean(s));
it = 0;
while nE ~= actnE || it > args.maxIter

    if actnE < nE
        s = [s(1) mean(s)];
    else
        s = [mean(s) s(2)];
    end
    actnE = nE_(mean(s));
    it = it + 1;
end

if it > args.maxIter
    warning('Maximum number of iterations reached, terminating.')
end

A = mean(s)*R <= P;
A = A + A';
G = graph(A);
