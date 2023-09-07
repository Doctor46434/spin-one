%% rd成像算法 自己写的  五个点的成像
% 如果不提前对频率进行处理，则需要对fft的结果做fftshift，同样对ifft的 对象 做ifftshift，这是为了保证构建的频域函数与信号在坐标轴相对应。
% 频域轴的改变和fftshift是相对应的，改变则不需要shift，Tp/2仅造成延时，与其他无关。
% 作为例子，在距离向采用频域轴移动，方位向采用shift
% 在时域初始-Tp/2，需要在之后进行补偿，无论是频域利用fft性质还是单纯的坐标轴移动。否则会发生错位
% 代码没有考虑在方位向加根据合成孔径长度的窗，数据截取本身相当于矩形窗。
% 在编码之前写出全时间坐标轴比较方便，参考bp的代码
clc;
% close all;
% clear all;
%% 参数设置
c=3e8;
H=2000;  %高度
theta=0.025;
v=100;
lamda=0.05;fc=c/lamda;
R0=1.5e4;
Tp=0.2e-6;
B=150e6;
fs=300e6;Fs=fs;
PRF=500;
PRT=1/PRF;
Kr=B/Tp;
Ls=theta*R0;
Nr=800;%800是怎么来的？根据场景Rx宽度预设
%% 成像场景设置
Rx=[15100,15200,15300,15300];
Ry=[-3,0,3,0];
%% 坐标定义
Ta=Ls/v;  %方位向时间
Na=ceil(Ta/PRT);
ta=((-Na/2:1:Na/2-1)*PRT).';  % 纵向时间
Fa=(-PRF/2:PRF/Na:PRF/2-PRF/Na).';
% 起始时刻确定 从最小距离位置
trmin=2*R0/c;
tr=trmin:1/fs:trmin+(Nr-1)/fs;
% Fr的定义和做不做搬移是对应的 从-到正 则需要搬移。
%Fr=-fs/2:fs/Nr:fs/2-fs/Nr; %距离向频率范围

Fr=linspace(0,Fs,length(tr)+1);
Fr=Fr(1:end-1);
p=logical(Fr>Fs/2);  %在-TP/2时是0，在0时刻，是B/2的频率，需要调整。
Fr(p)=Fr(p)-Fs;    

t=ones(Na,1)*tr; %时域基准信号
t0=t-trmin;
df=10e6;
%% 回拨模型
s1=zeros(Na,Nr);
for ii=1:length(Rx);
    Rn=sqrt(Rx(ii)^2+(v*ta-Ry(ii)).^2);
    tao=2*Rn/c;
    tt=t-tao*ones(1,Nr); % 新时间
    k = rectpuls(tt,Tp);
    sb=rectpuls(tt,Tp).*exp(1i*pi*Kr*tt.^2).*(exp(-1i*4*pi/lamda*Rn)*ones(1,Nr)).*exp(j*2*pi*df*t0);%没有考虑方位向的窗
    s1=s1+sb;
end
figure
mesh(tr*c/2,ta*v,real(s1));

%% 距离向fft
% Fs1=fftshift(fft(s1,Nr,2),2); %
Fs1=fft(s1,Nr,2);
%% 距离向脉冲压缩
s2=zeros(Na,Nr);
H1=exp(1i*pi*(ones(Na,1)*Fr).^2./Kr);%.*exp(1i*pi*(ones(Na,1)*Fr)*Tp)
Fs2=Fs1.*H1;
% s3=ifft(ifftshift(PFs1,2),Nr,2); %
s2=ifft(Fs2,Nr,2); %
figure
mesh(tr*c/2,ta*v,abs(s2));
%% 方位向fft，变为距离多普勒域
Fs3=fftshift(fft(s2,Na,1),1);
%% 距离徙动校正 
% 采用第五章的频域内校正弯曲的做法 式（5.50）
Fs4=fft(Fs3,Nr,2);
Rs=ones(Na,1)*tr*c/2;  %场景中心线最近距离
fam=2*v/lamda;  %最大多普勒偏移
H2=exp(1j*2*pi*Rs/c/fam/fam.*((Fa*ones(1,Nr)).^2).*(ones(Na,1)*Fr)); % 二维去偶函数
Fs5=Fs4.*H2;
s3=ifft(Fs5,Nr,2);
figure;
s8=ifft(ifftshift(s3,1),Na,1);
imagesc(tr*c/2,ta*v,abs(s8));
title('校正之后');
figure
mesh(tr*c/2,ta*v,abs(s3));
%% 方位向脉冲压缩
% 传递距离逆fft之后的矩阵 （s3）
% 仿真为正侧视，theta0近似为0；
H3=exp(1i*2*pi*Rs/v.*sqrt(fam.^2-(Fa*ones(1,Nr)).^2)); %没有做泰勒展开近似
Fs6=s3.*H3;
s6=ifft(ifftshift(Fs6,1),Na,1);
%% 结果显示
figure
imagesc(tr*c/2,ta*v,abs(s6));
% E=tuxiangs(s6);

% s7=chazhi(s6,8);
% figure
% imagesc(db(s7));
% 
% figure
% surf(db(s7(7300:7600,3120:3210)));

return
figure;
mesh(tr*c/2,ta*v,abs(s6));
figure
imagesc(tr*c/2,ta*v,abs(s6));
xlabel('斜距坐标（m）');
ylabel('方位坐标（m）');
title('目标');
%% 方位相位
stemp=s6(:,197);
Stemp=fft(stemp);
ang=angle(Stemp);
ang=unwrap(ang);
% ang=phase(Stemp);
figure;plot(ang11); 


% ang11=ang;
% anghe=ang;
ang_cha=anghe-ang11-ang;
figure;plot(ang_cha);
figure;plot(ang);