function dxdt=FNode(t,x, hip,Ihip,a,b,c)

dxdt=zeros(2,1);
if t<hip %lefojtó áram
    %dv=c*(x-x^3/3+y+I) 
    %dw=(1/c)*(a-x-b*y)
    dxdt(1)=c*(x(1)-x(1)^3/3+x(2)+Ihip);
else %normál egyenlet
    %dv=c*(x-x^3/3+y) 
    %dw=(1/c)*(a-x-b*y)
    dxdt(1)=c*(x(1)-x(1)^3/3+x(2));
end
dxdt(2)=(1/c)*(a-x(1)-b*x(2));
