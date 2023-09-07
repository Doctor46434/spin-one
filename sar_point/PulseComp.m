function [CompWave] = PulseComp(OriginWave)
% 距离向脉冲压缩
%% 参数
Tp = 20e-6;
B = 150e6;
Kr = B/Tp;
fs = 200e6;
c = 3e8;

Height = 3000;
WaveLength = 0.05;
AngleWaveWidth = 0.025;
Rmin = 15000;
PRF = 1000;
SpeedofFlight = 100;
R1 = 15200;
x1 = 0;
R2 = 15200;
x2 = 3;
R3 = 15250;
x3 = 0;


fc = c/WaveLength;

Rreal = sqrt(Height^2+R1^2);
Ls = AngleWaveWidth*Rreal;

Rreal2 = sqrt(Height^2+R2^2);
Ls2 = AngleWaveWidth*Rreal2;

Rreal3 = sqrt(Height^2+R3^2);
Ls3 = AngleWaveWidth*Rreal3;

MeterFreq = SpeedofFlight/PRF;
%% 算法
% Flight_x = -Ls3:MeterFreq:Ls2+x2;
Flight_x = -Ls:MeterFreq:Ls;
t = (-Tp*fs/2:Tp*fs/2)/fs;

RefSignal = exp(1j*pi*Kr.*(t).^2).*exp(1j*2*pi*fc.*(t));
CompWave = zeros(length(Flight_x),length(t));
for i = 1:length(Flight_x)
    CompWave(i,:) = ifftshift(ifft(fft(OriginWave(i,:)).*conj(fft(RefSignal))));
end

% 作图
% mesh(t*c/2+15700,Flight_x,(abs(CompWave)));
% 
% view(2);
% 
% xlabel("距离向/m");
% ylabel("方位向/m");
% title("单点回波距离向脉冲压缩图像");
end

