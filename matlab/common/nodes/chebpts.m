function x = chebpts(n, I)
% chebpts in the interval I
if(nargin == 1)
    I = [-1 1];
end

if(n == 0)
    x = [];     
elseif(n == 1)
    x = 0; 
else
    m = n - 1;
    x = sin(pi*(-m:2:m)/(2*m)).';
end
x = (I(2) - I(1))/2 * x + (I(2) + I(1))/2;
end
