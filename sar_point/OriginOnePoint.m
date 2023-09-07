function [OriginReWave] = OriginOnePoint
% 产生单点目标的回波信号
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

fc = c/WaveLength;

Rreal = sqrt(Height^2+R1^2);
Ls = AngleWaveWidth*Rreal;

MeterFreq = SpeedofFlight/PRF;

%% 算法
Flight_x = -Ls:MeterFreq:Ls;
t = (-Tp*fs/2:Tp*fs/2)/fs;
ReWave = zeros(length(Flight_x),length(t));
for i = 1:length(Flight_x)
    if (Flight_x(i)>(x1-Ls)) && (Flight_x(i)<(x1+Ls))
        Rtemp = sqrt((Flight_x(i)-x1)^2+R1^2+Height^2);
        ReWave(i,:) = exp(1j*pi*Kr.*(t-2*Rtemp/c).^2).*exp(1j*2*pi*fc.*(t-2*Rtemp/c));
    end
end

OriginReWave = ReWave;
% 作图
% mesh(t,Flight_x,real(ReWave));
% 
% view(2);
% 
% xlabel("距离向");
% ylabel("方位向");
% title("单点目标回波原始图像");
end