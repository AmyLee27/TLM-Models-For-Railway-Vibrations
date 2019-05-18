function  [f,recep] = time_to_recep(uh, dt, force, padding_multiple)

% padding for better frequency resolution
if(padding_multiple == 0)
    padding_length = 0;
else
    padding_length = length(uh)*10;
end
uh(length(uh):length(uh)+padding_length) = uh(end);
time = [1:length(uh)]*dt;


N = length(time);
At  =(time(end)-time(1))/(N-1);
fs = 1/At;
Af = fs/N;
ftwo = ((1:N)-1)*Af;
f = ftwo(1:N/2+1);

uh = detrend(uh);
u = gradient(uh,At);

nfft = size(u,1);
Uhelp = zeros(nfft,N);
U = zeros(nfft,round(N/2+1));
Uhelp = At*fft(u,[],2);     %factor At
U = Uhelp(:,1:N/2+1);

recep = abs(U)*force;
