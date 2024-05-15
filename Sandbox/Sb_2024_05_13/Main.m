%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2024. May 13. (2023a)
%

syms t u x y z real
r = [x;y;z];

B = [ 
     0 1 1 
     1 0 1 
     1 1 0 
    ];

Br = B*r;

dx = x*(Br(1) - r.'*B*r);
dy = y*(Br(2) - r.'*B*r);
dz = z*(Br(3) - r.'*B*r);

du = collect(simplify(subs(simplify(dx + dy),y,u-x)),[z,u])
dz = collect(simplify(subs(dz,y,u-x)),[z,u])

%%

f1 = matlabFunction(f_sym(1),'vars',[u;z]);
f2 = matlabFunction(f_sym(2),'vars',[u;z]);
f = matlabFunction(f_sym,'vars',{t,[u;z]});

[u_,z_] = meshgrid(linspace(0,1,101));

f1_ = f1(u_,z_);
f2_ = f2(u_,z_);

ldx = u_ + z_ > 1;
f1_(ldx) = NaN;
f2_(ldx) = NaN;
                
streamslice(u_,z_,f1_,f2_);