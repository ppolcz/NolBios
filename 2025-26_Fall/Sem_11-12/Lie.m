function ret = Lie(r,f,h,x)
    if r == 1
        ret = jacobian(h,x)*f;
    else
        ret = Lie(1,f,Lie(r-1,f,h,x),x);
    end
end