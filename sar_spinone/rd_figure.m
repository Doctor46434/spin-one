%% rd�����㷨 �Լ�д��  �����ĳ���
% �������ǰ��Ƶ�ʽ��д�������Ҫ��fft�Ľ����fftshift��ͬ����ifft�� ���� ��ifftshift������Ϊ�˱�֤������Ƶ�������ź������������Ӧ��
% Ƶ����ĸı��fftshift�����Ӧ�ģ��ı�����Ҫshift��Tp/2�������ʱ���������޹ء�
% ��Ϊ���ӣ��ھ��������Ƶ�����ƶ�����λ�����shift
% ��ʱ���ʼ-Tp/2����Ҫ��֮����в�����������Ƶ������fft���ʻ��ǵ������������ƶ�������ᷢ����λ
% ����û�п����ڷ�λ��Ӹ��ݺϳɿ׾����ȵĴ������ݽ�ȡ�����൱�ھ��δ���
% �ڱ���֮ǰд��ȫʱ��������ȽϷ��㣬�ο�bp�Ĵ���
clc;
% close all;
% clear all;
%% ��������
c=3e8;
H=2000;  %�߶�
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
Nr=800;%800����ô���ģ����ݳ���Rx���Ԥ��
%% ���񳡾�����
Rx=[15100,15200,15300,15300];
Ry=[-3,0,3,0];
%% ���궨��
Ta=Ls/v;  %��λ��ʱ��
Na=ceil(Ta/PRT);
ta=((-Na/2:1:Na/2-1)*PRT).';  % ����ʱ��
Fa=(-PRF/2:PRF/Na:PRF/2-PRF/Na).';
% ��ʼʱ��ȷ�� ����С����λ��
trmin=2*R0/c;
tr=trmin:1/fs:trmin+(Nr-1)/fs;
% Fr�Ķ���������������Ƕ�Ӧ�� ��-���� ����Ҫ���ơ�
%Fr=-fs/2:fs/Nr:fs/2-fs/Nr; %������Ƶ�ʷ�Χ

Fr=linspace(0,Fs,length(tr)+1);
Fr=Fr(1:end-1);
p=logical(Fr>Fs/2);  %��-TP/2ʱ��0����0ʱ�̣���B/2��Ƶ�ʣ���Ҫ������
Fr(p)=Fr(p)-Fs;    

t=ones(Na,1)*tr; %ʱ���׼�ź�
t0=t-trmin;
df=10e6;
%% �ز�ģ��
s1=zeros(Na,Nr);
for ii=1:length(Rx);
    Rn=sqrt(Rx(ii)^2+(v*ta-Ry(ii)).^2);
    tao=2*Rn/c;
    tt=t-tao*ones(1,Nr); % ��ʱ��
    k = rectpuls(tt,Tp);
    sb=rectpuls(tt,Tp).*exp(1i*pi*Kr*tt.^2).*(exp(-1i*4*pi/lamda*Rn)*ones(1,Nr)).*exp(j*2*pi*df*t0);%û�п��Ƿ�λ��Ĵ�
    s1=s1+sb;
end
figure
mesh(tr*c/2,ta*v,real(s1));

%% ������fft
% Fs1=fftshift(fft(s1,Nr,2),2); %
Fs1=fft(s1,Nr,2);
%% ����������ѹ��
s2=zeros(Na,Nr);
H1=exp(1i*pi*(ones(Na,1)*Fr).^2./Kr);%.*exp(1i*pi*(ones(Na,1)*Fr)*Tp)
Fs2=Fs1.*H1;
% s3=ifft(ifftshift(PFs1,2),Nr,2); %
s2=ifft(Fs2,Nr,2); %
figure
mesh(tr*c/2,ta*v,abs(s2));
%% ��λ��fft����Ϊ�����������
Fs3=fftshift(fft(s2,Na,1),1);
%% �����㶯У�� 
% ���õ����µ�Ƶ����У������������ ʽ��5.50��
Fs4=fft(Fs3,Nr,2);
Rs=ones(Na,1)*tr*c/2;  %�����������������
fam=2*v/lamda;  %��������ƫ��
H2=exp(1j*2*pi*Rs/c/fam/fam.*((Fa*ones(1,Nr)).^2).*(ones(Na,1)*Fr)); % ��άȥż����
Fs5=Fs4.*H2;
s3=ifft(Fs5,Nr,2);
figure;
s8=ifft(ifftshift(s3,1),Na,1);
imagesc(tr*c/2,ta*v,abs(s8));
title('У��֮��');
figure
mesh(tr*c/2,ta*v,abs(s3));
%% ��λ������ѹ��
% ���ݾ�����fft֮��ľ��� ��s3��
% ����Ϊ�����ӣ�theta0����Ϊ0��
H3=exp(1i*2*pi*Rs/v.*sqrt(fam.^2-(Fa*ones(1,Nr)).^2)); %û����̩��չ������
Fs6=s3.*H3;
s6=ifft(ifftshift(Fs6,1),Na,1);
%% �����ʾ
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
xlabel('б�����꣨m��');
ylabel('��λ���꣨m��');
title('Ŀ��');
%% ��λ��λ
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