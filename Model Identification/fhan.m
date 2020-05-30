function fh = fhan( x1, x2, r, h )

d = h*r;
d0 = h*d;
y = x1+h*x2;
a0 = sqrt(d^2+8*r*abs(y));

if abs(y)>d0
    a = x2+0.5*(a0-d)*sign(y);
else
    a = x2+y/h;
end

if abs(a)>d
    fh = -r*sign(a);
else
    fh = -r*a/d;
end

end

