Rref =1e3; 
A = 1;
c = 3e8;
fc = 9e9;
R1 = 930;
R2 = 1000;
R3 = 1050;
dR1 = 1000-R1;
dR2 = 1000-R2;
dR3 = 1000-R3;
Tp = 2e-6;
fs = 300e6;
B = 150e6;
Kr = B/Tp;
TpR = 3e-6;

t = -TpR:1/fs:TpR;

sref = A*rectpuls(t/TpR).*exp(1j*2*pi.*(fc.*t+1/2*Kr*t.^2));
s1 = A*rectpuls((t-(2*dR1/c))/Tp).*exp(1j*2*pi*(fc*(t-(2*dR1/c))+1/2*Kr*(t-(2*dR1/c)).^2));
s2 = A*rectpuls((t-(2*dR2/c))/Tp).*exp(1j*2*pi*(fc*(t-(2*dR2/c))+1/2*Kr*(t-(2*dR2/c)).^2));
s3 = A*rectpuls((t-(2*dR3/c))/Tp).*exp(1j*2*pi*(fc*(t-(2*dR3/c))+1/2*Kr*(t-(2*dR3/c)).^2));

% plot(t,sref,t,s1);

sif1 = conj(sref).*s1;
sif2 = conj(sref).*s2;
sif3 = conj(sref).*s3;
figure;
plot(t,sif1,t,sif2,t,sif3);

Fsif1 = fftshift(fft(sif1));
Fsif2 = fftshift(fft(sif2));
Fsif3 = fftshift(fft(sif3));
f = -fs/2:1/TpR/2:fs/2;
figure;
plot(f/2/Kr*c,abs(Fsif1),f/2/Kr*c,abs(Fsif2),f/2/Kr*c,abs(Fsif3));

Sef = exp(-1j*pi*f.^2/Kr);
Fsif11 = Fsif1.*Sef;
Fsif21 = Fsif2.*Sef;
Fsif31 = Fsif3.*Sef;
figure;
plot(f/2/Kr*c,abs(Fsif11),f/2/Kr*c,abs(Fsif21),f/2/Kr*c,abs(Fsif31));

sif11 = ifft(ifftshift(Fsif11));
sif21 = ifft(ifftshift(Fsif21));
sif31 = ifft(ifftshift(Fsif31));
figure;
plot(t,sif11,t,sif21,t,sif31);