%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2024. November 05. (2023a)
%

syms a b real


phi = [
    1
    a
    b
    a^2
    a*b
    b^2
    ];


P = [
    2  0  0 0 0 0
    0  0 -3 0 0 1
    0 -3  0 1 0 0
    0  0  1 0 0 0
    0  0  0 0 0 0
    0  1  0 0 0 0
    ] / 2;

Motzkin = a^2*b + a*b^2 - 3*a*b + 1

simplify(phi.' * P * phi)
