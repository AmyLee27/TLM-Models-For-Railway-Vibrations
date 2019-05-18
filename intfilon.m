function  int = intfilon(x,s,f)

n = length(x) - 1;
sx = s.'*x;
S  = repmat(s.',1,length(x));

b = diff(f) ./ diff(x);
a = f(1:n) - b .* x(1:n);


ssx = sin(sx);
csx = cos(sx);
icon = ssx ./ S;
ilin = ( sx .* ssx + csx)  ./ S.^2;

%
con =[-a(1) -diff(a)  a(n)];
lin =[-b(1) -diff(b)  b(n)];
%
int = con * icon.' + lin * ilin.';
