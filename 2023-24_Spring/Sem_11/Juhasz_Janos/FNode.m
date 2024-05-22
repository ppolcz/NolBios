function dxdt=FNode(t,x, hip,Ihip,a,b,c)

dxdt=zeros(2,1);
if t<hip %lefojtó áram
    %dv=c*(x-x^3/3+y+I) 
    %dw=(1/c)*(a-x-b*y)
    dxdt(1)=;
else %normál egyenlet
    %dv=c*(x-x^3/3+y) 
    %dw=(1/c)*(a-x-b*y)
    dxdt(1)=;
end
dxdt(2)=;
