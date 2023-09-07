% 合成孔径rd成像
%% 三点回波的参数列表
H = 3000; %飞行高度
lambda = 0.05; %载波波长
A = 0.025; %方位向波束宽度
Rmin = 15000; %初始采样距离
Tp = 0.2e-6; %脉冲宽度
B = 150e6; %频带宽度
fs = 200e6; %采样率
PRF = 1000; %脉冲重复频率
v0 = 100; %飞机飞行速度
Rx = [15200 15200 15250];  %点目标垂直距离
Ry = [0 3 0]; %点目标沿航向坐标
%% 其余参数
c = 3e8; %光速
Kr = B/Tp; %LFM信号的调频率
RB = Rx;
Ls = A*Rmin;
N = Ls/v0;
Np = ceil(N*PRF);
Nr = 800;
tmin = 2*Rmin/c;
t = tmin:1/fs:tmin+(Nr-1)*1/fs;
tt = ones(Np,1)*t;
tv = (-Np/2:1:Np/2-1).'/PRF;

Fr=linspace(0,fs,length(t)+1);
Fr=Fr(1:end-1);
p=logical(Fr>fs/2);  %在-TP/2时是0，在0时刻，是B/2的频率，需要调整。
Fr(p)=Fr(p)-fs;    

fc = c/lambda;

fam = 2*v0/lambda;
Fa=(-PRF/2:PRF/Np:PRF/2-PRF/Np).';
Rs=ones(Np,1)*t*c/2;
%% 回波模型
sb = zeros(Np,Nr);
for i = 1:length(RB)
    Ri = sqrt(RB(i)^2 + (Ry(i)-v0.*tv).^2);
    tau = 2*Ri/c;
    si = rectpuls(tt-tau*ones(1,Nr),Tp).*exp(1j*pi*Kr.*(tt-tau*ones(1,Nr)).^2).*exp(1j*2*pi*fc.*(tt-tau*ones(1,Nr)));
    sb = sb + si;
end
figure
mesh(t*c/2,tv*v0,real(sb));
xlabel("距离向/m");
ylabel("方位向/m");
title("原始回波结果（实部）");
%% 距离脉压
F_sb =fft(sb,Nr,2);
figure;
mesh(t*c/2,tv*v0,real(F_sb));
temp = ones(Np,1)*((-Nr/2:1:Nr/2-1)/fs);
BaseFunc = rectpuls(temp,Tp).*exp(-1j*pi*Kr.*(tt).^2);

F_BaseFunc = fftshift(fft(BaseFunc,Nr,2));
% F_BaseFunc = exp(1i*pi*(ones(Np,1)*Fr).^2./Kr);
k = F_sb.*F_BaseFunc;
s1 = fftshift(ifft(k,Nr,2));
figure;
mesh(t*c/2,tv*v0,abs(s1));
xlabel("距离向/m");
ylabel("方位向/m");
title("距离向脉冲压缩结果");
%% 频域矫正距离弯曲
s2 = fftshift(fft(s1,Np,1),1);
s3 = fft(s2,Nr,2);
Fa=(-PRF/2:PRF/Np:PRF/2-PRF/Np).';
figure;
mesh(t*c/2,Fa,abs(s2));
xlabel("距离向/m");
ylabel("方位向/Hz");
title("变换到方位频域后结果");
H2=exp(1j*2*pi*Rs/c/fam/fam.*((Fa*ones(1,Nr)).^2).*(ones(Np,1)*Fr));
s4 = H2.*s3;
s5 = ifft(s4,Nr,2);
figure;
mesh(t*c/2,tv*v0,abs(s5));
xlabel("距离向/m");
ylabel("方位向/Hz");
title("频域矫正结果");
s6 = ifft(ifftshift(s5,1),Np,1);
figure;
mesh(t*c/2,tv*v0,abs(s6));

%% 方位向脉压
F_s1 = fft(s6);
tnew = tv*ones(1,Nr);
ABaseFunc = exp(-1j*pi*(-2*v0^2/lambda/Rmin).*(tnew).^2);
F_ABaseFunc = fft(ABaseFunc);
Fs6=F_s1.*F_ABaseFunc;
s6=ifft(Fs6,Np,1);

figure;
mesh(t*c/2,tv*v0,real(s6));
xlabel("距离向/m");
ylabel("方位向/m");
title("方位向脉冲压缩后结果");
figure;
imagesc(t*c/2,tv*v0,real(s6));
