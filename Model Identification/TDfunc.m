function [Xtrack, Xdiff] = TDfunc(x, r, h, h0, h1)

[rX, cX]=size(x);
Xtrack = zeros(rX, cX);  % tracking signal
XPredc = zeros(rX, cX);  % predicted signal
Xdiff = zeros(rX, cX);  % differential signal

for n = 1:cX
    for m = 1:rX-1
        XPredc(m,n) = x(m,n)+h1*Xdiff(m,n);
        fh = fhan( Xtrack(m,n)-XPredc(m,n), Xdiff(m,n), r, h0);
        Xtrack(m+1,n) = Xtrack(m,n)+h*Xdiff(m,n);
        Xdiff(m+1,n) = Xdiff(m,n)+h*fh;
    end
end
end

