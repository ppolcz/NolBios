function [rloc] = findroots(y,x)

idx = find(abs(diff(sign(y))));

x0 = x(idx);
x1 = x(idx+1);

y0 = y(idx);
y1 = y(idx+1);

rloc = x0 - y0 ./ (y1 - y0) .* (x1 - x0);

ldx = ismissing(rloc);
rloc(ldx) = (x0(ldx) + x1(ldx)) / 2;



