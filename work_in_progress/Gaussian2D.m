
matr = zeros(200,200);
gsize = size(matr);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(R,C, [50], [100,100]);


function val = gaussC(x, y, sigma, center)
xc = center(1);
yc = center(2);
exponent = ((x-xc).^2 + (y-yc).^2)./(1*sigma^2);
val       = (exp(-exponent));
end