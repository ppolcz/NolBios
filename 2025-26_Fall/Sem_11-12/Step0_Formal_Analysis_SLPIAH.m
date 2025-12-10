syms beta nu omega real
syms S L P I A H D R real
syms tau_L tau_P tau_I tau_A tau_H p_I p_H p_D q_A real
syms N_p real
syms v real
%[text] The state dynamics. (The recovered compartment is only an accumulator, but it does not affect the overall dynamics, therefore, it is neglected for simplicity.)
dS = -beta*(P+I+q_A*A)*S/N_p - nu*S + omega*R;
dL = beta*(P+I+q_A*A)*S/N_p - tau_L*L;
dP = tau_L*L - tau_P*P;
dI = p_I*tau_P*P - tau_I*I;
dA = (1-p_I)*tau_P*P - tau_A*A;
dH = tau_I*p_H*I - tau_H*H;
dD = p_D*tau_H*H;
dR = tau_I*(1-p_H)*I + tau_A*A + (1-p_D)*tau_H*H + nu*S - omega*R;

u = beta;
x = [S ; L ; P ; I ; A ; H ; D ; R];
h = H;

F = [dS ; dL ; dP ; dI ; dA ; dH ; dD ; dR] %[output:140cbe5f]
f = subs(F,beta,0) %[output:59aa7822]
g = (F - f) / beta %[output:6e099020]
%%
%[text] ## The output function: $\\Phi\_1(x) = h(x) = \\mathbf{I}$
Phi1 = h %[output:3f07f5c1]
%[text] The first derivative: $\\Phi\_2(x) = L\_f h(x) + L\_g h(x) \\, u$, where $L\_g h(x) = 0$
Phi2 = Lie(1,f,Phi1,x) + Lie(1,g,Phi1,x)*u %[output:21ea26fd]
%[text] The second derivative: $\\Phi\_3(x) = L\_f^2 h(x) + L\_g L\_f h(x) \\, u$, where $L\_g L\_f h(x) = 0$
Phi3 = Lie(1,f,Phi2,x) + Lie(1,g,Phi2,x)*u %[output:21c08d1c]
%[text] The third derivative: $\\Phi\_4(x) = L\_f^3 h(x) + L\_g L\_f^2 h(x) \\, u$, where $L\_g L\_f^2 h(x) = 0$
Phi4 = Lie(1,f,Phi3,x) + Lie(1,g,Phi3,x)*u %[output:5398582d]
%[text] The fourth derivative: $\\Phi\_5(x) = L^4\_f h(x) + L\_g L\_f^3 h(x) \\, u$, where $L\_g L\_f^3 h(x) \\neq 0$, therefore, the relative degree with respect to the output $y = h(x) = \\textbf{I}$ is $r = 4\n$
Phi5 = Lie(1,f,Phi4,x) + Lie(1,g,Phi4,x)*u %[output:79b85040]
Phi5 = collect(Phi5,[x;u]) %[output:31ada6b4]
%[text] The feedback linearizing input: $u = \\frac{v - L\_f^4 h(x)}{L\_g L\_f^3 h(x)}$
u_FdbLin = Lie(1,g,Phi4,x) \ (v - Lie(1,f,Phi4,x)) %[output:94a72e00]
%[text] As the relative degree is 4, the fifth output transformation $\\Phi\_5(x)$ is not used in the final coordinates transformation, but $\\Phi(x) = \\pmatrix{ \\Phi\_\\xi(x) \\cr \\Phi\_\\eta(x) }$, where $\\Phi\_\\xi(x) = \\pmatrix{ \\Phi\_1(x) \\cr \\Phi\_2(x) \\cr \\Phi\_3(x) \\cr \\Phi\_4(x) }$
Phi_xi = [Phi1 ; Phi2 ; Phi3 ; Phi4] %[output:4f64b4a9]
T_x = jacobian(Phi_xi,x) % [2022_Csutak.etal, Eqs. (7)-(9)] %[output:998750c7]
simplify(T_x * x - Phi_xi) %[output:97041532]
%%
%[text] Then, we select $\\Phi\_\\eta(x) = \\pmatrix{ \\mathbf{S} + \\mathbf{L} \\cr \\mathbf{A} \\cr \\mathbf{D} \\cr \\mathbf{R}}$ and check whether $L\_g \\Phi\_\\eta(x) = 0$.
Phi_eta = [S+L ; A ; D ; R];
Lg_Phi_eta = Lie(1,g,Phi_eta,x) %[output:4935dbd6]
Phi = [ %[output:group:4171a549] %[output:2ca47b1b]
    Phi_xi %[output:2ca47b1b]
    Phi_eta %[output:2ca47b1b]
    ] %[output:group:4171a549] %[output:2ca47b1b]
J_Phi = jacobian(Phi,x) %[output:99b96db3]
det(J_Phi) %[output:3a1c200e]
%%
%[text] Find inverse transformation
r = 4;
n = 8;
xi = sym('xi',[r,1]);
eta = sym('eta',[n-r,1]);

z = [xi;eta] %[output:0d0da1cb]
iPhi_sol = solve(Phi - z, x,'Real',true) %[output:4ad90f13]
%%
%[text] ## Zero dynamics
%[text] $\\dot \\eta = \\dot \\mathbf{S} + \\dot \\mathbf{L}  = \\left\[ -\\alpha \\,\\mathbf{L} \\right\]\_{\\xi = 0} = 0$, therefore, the zero dynamics is stable (but not asymptotically stable).
q = subs(subs(Lie(1,f,Phi_eta,x),iPhi_sol),xi,xi*0) %[output:88b5c3e5]
Aq = jacobian(q,eta) %[output:6ed2f98c]
eig(Aq) %[output:27e0a488]
%%
function ret = Lie(r,f,h,x)
    if r == 1
        ret = jacobian(h,x)*f;
    else
        ret = Lie(1,f,Lie(r-1,f,h,x),x);
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:140cbe5f]
%   data: {"dataType":"symbolic","outputData":{"name":"F","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\nR\\,\\omega -S\\,\\nu -\\sigma_1 \\\\\n\\sigma_1 -L\\,\\tau_L \\\\\nL\\,\\tau_L -P\\,\\tau_P \\\\\nP\\,p_I \\,\\tau_P -\\textrm{I}\\,\\tau_I \\\\\n-A\\,\\tau_A -P\\,\\tau_P \\,{\\left(p_I -1\\right)}\\\\\n\\textrm{I}\\,p_H \\,\\tau_I -H\\,\\tau_H \\\\\nH\\,p_D \\,\\tau_H \\\\\nA\\,\\tau_A -R\\,\\omega +S\\,\\nu -H\\,\\tau_H \\,{\\left(p_D -1\\right)}-\\textrm{I}\\,\\tau_I \\,{\\left(p_H -1\\right)}\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{S\\,\\beta \\,{\\left(\\textrm{I}+P+A\\,q_A \\right)}}{N_p }\n\\end{array}"}}
%---
%[output:59aa7822]
%   data: {"dataType":"symbolic","outputData":{"name":"f","value":"\\left(\\begin{array}{c}\nR\\,\\omega -S\\,\\nu \\\\\n-L\\,\\tau_L \\\\\nL\\,\\tau_L -P\\,\\tau_P \\\\\nP\\,p_I \\,\\tau_P -\\textrm{I}\\,\\tau_I \\\\\n-A\\,\\tau_A -P\\,\\tau_P \\,{\\left(p_I -1\\right)}\\\\\n\\textrm{I}\\,p_H \\,\\tau_I -H\\,\\tau_H \\\\\nH\\,p_D \\,\\tau_H \\\\\nA\\,\\tau_A -R\\,\\omega +S\\,\\nu -H\\,\\tau_H \\,{\\left(p_D -1\\right)}-\\textrm{I}\\,\\tau_I \\,{\\left(p_H -1\\right)}\n\\end{array}\\right)"}}
%---
%[output:6e099020]
%   data: {"dataType":"symbolic","outputData":{"name":"g","value":"\\left(\\begin{array}{c}\n-\\frac{S\\,{\\left(\\textrm{I}+P+A\\,q_A \\right)}}{N_p }\\\\\n\\frac{S\\,{\\left(\\textrm{I}+P+A\\,q_A \\right)}}{N_p }\\\\\n0\\\\\n0\\\\\n0\\\\\n0\\\\\n0\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:3f07f5c1]
%   data: {"dataType":"symbolic","outputData":{"name":"Phi1","value":"H"}}
%---
%[output:21ea26fd]
%   data: {"dataType":"symbolic","outputData":{"name":"Phi2","value":"\\textrm{I}\\,p_H \\,\\tau_I -H\\,\\tau_H"}}
%---
%[output:21c08d1c]
%   data: {"dataType":"symbolic","outputData":{"name":"Phi3","value":"\\tau_H \\,{\\left(H\\,\\tau_H -\\textrm{I}\\,p_H \\,\\tau_I \\right)}-p_H \\,\\tau_I \\,{\\left(\\textrm{I}\\,\\tau_I -P\\,p_I \\,\\tau_P \\right)}"}}
%---
%[output:5398582d]
%   data: {"dataType":"symbolic","outputData":{"name":"Phi4","value":"{\\left(\\textrm{I}\\,\\tau_I -P\\,p_I \\,\\tau_P \\right)}\\,{\\left(p_H \\,{\\tau_I }^2 +p_H \\,\\tau_H \\,\\tau_I \\right)}-{\\tau_H }^2 \\,{\\left(H\\,\\tau_H -\\textrm{I}\\,p_H \\,\\tau_I \\right)}+p_H \\,p_I \\,\\tau_I \\,\\tau_P \\,{\\left(L\\,\\tau_L -P\\,\\tau_P \\right)}"}}
%---
%[output:79b85040]
%   data: {"dataType":"symbolic","outputData":{"name":"Phi5","value":"\\begin{array}{l}\n{\\tau_H }^3 \\,{\\left(H\\,\\tau_H -\\textrm{I}\\,p_H \\,\\tau_I \\right)}-{\\left(\\tau_I \\,\\sigma_1 +p_H \\,{\\tau_H }^2 \\,\\tau_I \\right)}\\,{\\left(\\textrm{I}\\,\\tau_I -P\\,p_I \\,\\tau_P \\right)}-{\\left(L\\,\\tau_L -P\\,\\tau_P \\right)}\\,{\\left(p_H \\,p_I \\,\\tau_I \\,{\\tau_P }^2 +p_I \\,\\sigma_1 \\,\\tau_P \\right)}-L\\,p_H \\,p_I \\,\\tau_I \\,{\\tau_L }^2 \\,\\tau_P +\\frac{S\\,\\beta \\,p_H \\,p_I \\,\\tau_I \\,\\tau_L \\,\\tau_P \\,{\\left(\\textrm{I}+P+A\\,q_A \\right)}}{N_p }\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =p_H \\,{\\tau_I }^2 +p_H \\,\\tau_H \\,\\tau_I \n\\end{array}"}}
%---
%[output:31ada6b4]
%   data: {"dataType":"symbolic","outputData":{"name":"Phi5","value":"\\begin{array}{l}\n\\sigma_1 \\,S\\,P\\,\\beta +\\sigma_1 \\,S\\,\\textrm{I}\\,\\beta +\\frac{p_H \\,p_I \\,q_A \\,\\tau_I \\,\\tau_L \\,\\tau_P }{N_p }\\,S\\,A\\,\\beta +{\\left(-\\tau_L \\,\\sigma_2 -p_H \\,p_I \\,\\tau_I \\,{\\tau_L }^2 \\,\\tau_P \\right)}\\,L+{\\left(\\tau_P \\,\\sigma_2 +p_I \\,\\tau_P \\,\\sigma_3 \\right)}\\,P+{\\left(-\\tau_I \\,\\sigma_3 -p_H \\,{\\tau_H }^3 \\,\\tau_I \\right)}\\,\\textrm{I}+{\\tau_H }^4 \\,H\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{p_H \\,p_I \\,\\tau_I \\,\\tau_L \\,\\tau_P }{N_p }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =p_H \\,p_I \\,\\tau_I \\,{\\tau_P }^2 +p_I \\,\\sigma_4 \\,\\tau_P \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =\\tau_I \\,\\sigma_4 +p_H \\,{\\tau_H }^2 \\,\\tau_I \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =p_H \\,{\\tau_I }^2 +p_H \\,\\tau_H \\,\\tau_I \n\\end{array}"}}
%---
%[output:94a72e00]
%   data: {"dataType":"symbolic","outputData":{"name":"u_FdbLin","value":"\\begin{array}{l}\n\\frac{N_p \\,{\\left(v-{\\tau_H }^3 \\,{\\left(H\\,\\tau_H -\\textrm{I}\\,p_H \\,\\tau_I \\right)}+{\\left(\\tau_I \\,\\sigma_1 +p_H \\,{\\tau_H }^2 \\,\\tau_I \\right)}\\,{\\left(\\textrm{I}\\,\\tau_I -P\\,p_I \\,\\tau_P \\right)}+{\\left(L\\,\\tau_L -P\\,\\tau_P \\right)}\\,{\\left(p_H \\,p_I \\,\\tau_I \\,{\\tau_P }^2 +p_I \\,\\sigma_1 \\,\\tau_P \\right)}+L\\,p_H \\,p_I \\,\\tau_I \\,{\\tau_L }^2 \\,\\tau_P \\right)}}{S\\,p_H \\,p_I \\,\\tau_I \\,\\tau_L \\,\\tau_P \\,{\\left(\\textrm{I}+P+A\\,q_A \\right)}}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =p_H \\,{\\tau_I }^2 +p_H \\,\\tau_H \\,\\tau_I \n\\end{array}"}}
%---
%[output:4f64b4a9]
%   data: {"dataType":"symbolic","outputData":{"name":"Phi_xi","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\nH\\\\\n\\textrm{I}\\,p_H \\,\\tau_I -H\\,\\tau_H \\\\\n\\tau_H \\,\\sigma_2 -p_H \\,\\tau_I \\,\\sigma_1 \\\\\n\\sigma_1 \\,{\\left(p_H \\,{\\tau_I }^2 +p_H \\,\\tau_H \\,\\tau_I \\right)}-{\\tau_H }^2 \\,\\sigma_2 +p_H \\,p_I \\,\\tau_I \\,\\tau_P \\,{\\left(L\\,\\tau_L -P\\,\\tau_P \\right)}\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\textrm{I}\\,\\tau_I -P\\,p_I \\,\\tau_P \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =H\\,\\tau_H -\\textrm{I}\\,p_H \\,\\tau_I \n\\end{array}"}}
%---
%[output:998750c7]
%   data: {"dataType":"symbolic","outputData":{"name":"T_x","value":"\\begin{array}{l}\n\\left(\\begin{array}{cccccccc}\n0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\\\\n0 & 0 & 0 & p_H \\,\\tau_I  & 0 & -\\tau_H  & 0 & 0\\\\\n0 & 0 & p_H \\,p_I \\,\\tau_I \\,\\tau_P  & -p_H \\,{\\tau_I }^2 -p_H \\,\\tau_H \\,\\tau_I  & 0 & {\\tau_H }^2  & 0 & 0\\\\\n0 & p_H \\,p_I \\,\\tau_I \\,\\tau_L \\,\\tau_P  & -p_H \\,p_I \\,\\tau_I \\,{\\tau_P }^2 -p_I \\,\\sigma_1 \\,\\tau_P  & \\tau_I \\,\\sigma_1 +p_H \\,{\\tau_H }^2 \\,\\tau_I  & 0 & -{\\tau_H }^3  & 0 & 0\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =p_H \\,{\\tau_I }^2 +p_H \\,\\tau_H \\,\\tau_I \n\\end{array}"}}
%---
%[output:97041532]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n0\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:4935dbd6]
%   data: {"dataType":"symbolic","outputData":{"name":"Lg_Phi_eta","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n0\\\\\n0\n\\end{array}\\right)"}}
%---
%[output:2ca47b1b]
%   data: {"dataType":"symbolic","outputData":{"name":"Phi","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\nH\\\\\n\\textrm{I}\\,p_H \\,\\tau_I -H\\,\\tau_H \\\\\n\\tau_H \\,\\sigma_2 -p_H \\,\\tau_I \\,\\sigma_1 \\\\\n\\sigma_1 \\,{\\left(p_H \\,{\\tau_I }^2 +p_H \\,\\tau_H \\,\\tau_I \\right)}-{\\tau_H }^2 \\,\\sigma_2 +p_H \\,p_I \\,\\tau_I \\,\\tau_P \\,{\\left(L\\,\\tau_L -P\\,\\tau_P \\right)}\\\\\nL+S\\\\\nA\\\\\n\\textrm{D}\\\\\nR\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\textrm{I}\\,\\tau_I -P\\,p_I \\,\\tau_P \\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =H\\,\\tau_H -\\textrm{I}\\,p_H \\,\\tau_I \n\\end{array}"}}
%---
%[output:99b96db3]
%   data: {"dataType":"symbolic","outputData":{"name":"J_Phi","value":"\\begin{array}{l}\n\\left(\\begin{array}{cccccccc}\n0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\\\\n0 & 0 & 0 & p_H \\,\\tau_I  & 0 & -\\tau_H  & 0 & 0\\\\\n0 & 0 & p_H \\,p_I \\,\\tau_I \\,\\tau_P  & -p_H \\,{\\tau_I }^2 -p_H \\,\\tau_H \\,\\tau_I  & 0 & {\\tau_H }^2  & 0 & 0\\\\\n0 & p_H \\,p_I \\,\\tau_I \\,\\tau_L \\,\\tau_P  & -p_H \\,p_I \\,\\tau_I \\,{\\tau_P }^2 -p_I \\,\\sigma_1 \\,\\tau_P  & \\tau_I \\,\\sigma_1 +p_H \\,{\\tau_H }^2 \\,\\tau_I  & 0 & -{\\tau_H }^3  & 0 & 0\\\\\n1 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\\\\n0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\\\\n0 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\\\\n0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =p_H \\,{\\tau_I }^2 +p_H \\,\\tau_H \\,\\tau_I \n\\end{array}"}}
%---
%[output:3a1c200e]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-{p_H }^3 \\,{p_I }^2 \\,{\\tau_I }^3 \\,\\tau_L \\,{\\tau_P }^2"}}
%---
%[output:0d0da1cb]
%   data: {"dataType":"symbolic","outputData":{"name":"z","value":"\\left(\\begin{array}{c}\n\\xi_1 \\\\\n\\xi_2 \\\\\n\\xi_3 \\\\\n\\xi_4 \\\\\n\\eta_1 \\\\\n\\eta_2 \\\\\n\\eta_3 \\\\\n\\eta_4 \n\\end{array}\\right)"}}
%---
%[output:4ad90f13]
%   data: {"dataType":"textualVariable","outputData":{"header":"struct with fields:","name":"iPhi_sol","value":"    S: -(xi4 + tau_H*xi3 + tau_I*xi3 + tau_P*xi3 + tau_H*tau_I*xi2 + tau_H*tau_P*xi2 + tau_I*tau_P*xi2 + tau_H*tau_I*tau_P*xi1 - eta1*p_H*p_I*tau_I*tau_L*tau_P)\/(p_H*p_I*tau_I*tau_L*tau_P)\n    L: (xi4 + tau_H*xi3 + tau_I*xi3 + tau_P*xi3 + tau_H*tau_I*xi2 + tau_H*tau_P*xi2 + tau_I*tau_P*xi2 + tau_H*tau_I*tau_P*xi1)\/(p_H*p_I*tau_I*tau_L*tau_P)\n    P: (xi3 + tau_H*xi2 + tau_I*xi2 + tau_H*tau_I*xi1)\/(p_H*p_I*tau_I*tau_P)\n    I: (xi2 + tau_H*xi1)\/(p_H*tau_I)\n    A: eta2\n    H: xi1\n    D: eta3\n    R: eta4\n"}}
%---
%[output:88b5c3e5]
%   data: {"dataType":"symbolic","outputData":{"name":"q","value":"\\left(\\begin{array}{c}\n\\eta_4 \\,\\omega -\\eta_1 \\,\\nu \\\\\n-\\eta_2 \\,\\tau_A \\\\\n0\\\\\n\\eta_1 \\,\\nu -\\eta_4 \\,\\omega +\\eta_2 \\,\\tau_A \n\\end{array}\\right)"}}
%---
%[output:6ed2f98c]
%   data: {"dataType":"symbolic","outputData":{"name":"Aq","value":"\\left(\\begin{array}{cccc}\n-\\nu  & 0 & 0 & \\omega \\\\\n0 & -\\tau_A  & 0 & 0\\\\\n0 & 0 & 0 & 0\\\\\n\\nu  & \\tau_A  & 0 & -\\omega \n\\end{array}\\right)"}}
%---
%[output:27e0a488]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"\\left(\\begin{array}{c}\n0\\\\\n0\\\\\n-\\tau_A \\\\\n-\\nu -\\omega \n\\end{array}\\right)"}}
%---
