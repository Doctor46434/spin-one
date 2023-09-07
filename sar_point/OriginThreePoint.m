function [OriginReWave] = OriginThreePoint
% 函数功能为生成三个点目标的回波信号
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

Rreal1 = sqrt(Height^2+R1^2);
Ls1 = AngleWaveWidth*Rreal1;

Rreal2 = sqrt(Height^2+R2^2);
Ls2 = AngleWaveWidth*Rreal2;

Rreal3 = sqrt(Height^2+R3^2);
Ls3 = AngleWaveWidth*Rreal3;

MeterFreq = SpeedofFlight/PRF;

%% 算法
Flight_x = -Ls3:MeterFreq:Ls2+x2;
t = (-Tp*fs/2:Tp*fs/2)/fs;
ReWave = zeros(length(Flight_x),length(t));
% 将三个点目标信号叠加
for i = 1:length(Flight_x)
    if (Flight_x(i)>(x1-Ls1)) && (Flight_x(i)<(x1+Ls1))
        Rtemp = sqrt((Flight_x(i)-x1)^2+R1^2+Height^2);
        ReWave(i,:) = ReWave(i,:) + exp(1j*pi*Kr.*(t-2*Rtemp/c).^2).*exp(1j*2*pi*fc.*(t-2*Rtemp/c));
    end

    if (Flight_x(i)>(x2-Ls2)) && (Flight_x(i)<(x2+Ls2))
        Rtemp = sqrt((Flight_x(i)-x2)^2+R2^2+Height^2);
        ReWave(i,:) = ReWave(i,:) + exp(1j*pi*Kr.*(t-2*Rtemp/c).^2).*exp(1j*2*pi*fc.*(t-2*Rtemp/c));
    end

    if (Flight_x(i)>(x3-Ls3)) && (Flight_x(i)<(x3+Ls3))
        Rtemp = sqrt((Flight_x(i)-x3)^2+R3^2+Height^2);
        ReWave(i,:) = ReWave(i,:) + exp(1j*pi*Kr.*(t-2*Rtemp/c).^2).*exp(1j*2*pi*fc.*(t-2*Rtemp/c));
    end
end

OriginReWave = ReWave;
% figure;
% mesh(t,Flight_x,real(ReWave));
% 
% view(2);
% 
% xlabel("距离向");
% ylabel("方位向");
% title("单点目标回波原始图像");

end