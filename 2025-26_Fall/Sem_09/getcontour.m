function [C_out,H] = getcontour(x, y, z, alpha, varargin)
%% 

[C,handle] = contour(x,y,z,alpha);
delete(handle);

if nargout > 0 && isempty(C)
    C_out = C; 
    H = [];
    return 
end

N = size(C,2);

I = 1;
i = 1;
while I <= N
    C_out(i).level = C(1,I);
    C_out(i).nr = C(2,I);
    C_out(i).contour = C(:,I+1:I+C(2,I));
    
    H(i) = plot(C_out(i).contour(1,:), C_out(i).contour(2,:), varargin{:});

    I = I+C(2,I)+1;
    i = i+1;
end


end