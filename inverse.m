
function [s,dis]=inverse(k11,df)
%k11 - list of wave numbers on x direction;
%df - displacements of the beam on the transformed domain;
smax=1/(2*abs(k11(1)-k11(2))); sinc=1/(2*abs(k11(1)));
s=-2*pi*smax:2*pi*sinc:2*pi*smax;
u1=-ifftshift(df);
u2=ifft(u1);
u3=fftshift(u2);
dis=u3*length(k11)*abs(k11(1)-k11(2))/(2*pi);
%plot(s,dis);
%[hb,k11,x]=track(k1max,k1inc,k2max,k2inc,E,niu,del,m,hr,n,Eh,niuh,delh,mh,bon,b,EI,mb,de,v);
%[pk]=train(P,xl,k11);
%for j=1:length(k11),
%    df(j)=hb(j)*pk(j);
%end



