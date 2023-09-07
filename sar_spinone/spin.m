%% 参数
Rref =2e2; 
A = 1;
c = 3e8;
fc = 9e9;
omega = 10*pi;  %目标物体转速
K = 5;  %目标点数量
alpha = pi/6;  %观测角度
R = [1 1.5 2 3 3]*10;  %旋转半径
theta = [pi/2 0 pi/2 0 pi/2]; %初始相位
z = [0 0 20 20 50]; %初始z

Tp = 2e-6;
fs = 300e6;
B = 150e6;
Kr = B/Tp;
TpR = 3e-6;
PRF = 1000;

%% 构造回波

t = ones(PRF,1)*(-TpR:1/fs:TpR);   %%二维时间
tm = (0:PRF-1)/PRF;
s1 = zeros(PRF,2*TpR*fs+1);

for i=4
    R_delta = Rref + R(i)*cos(omega*tm+theta(i))*sin(alpha) - z(i)*cos(alpha);
    tt = t - 2*R_delta'/c*ones(1,2*TpR*fs+1);
    s1 = s1 + rectpuls(tt,Tp).*exp(1j*pi*Kr*tt.^2).*exp(1j*pi*fc*tt);
end

figure;
mesh(real(s1));
%% 脉压处理

t_ref = t - 2*Rref'/c*ones(1,2*TpR*fs+1);
s_ref = rectpuls(t_ref,TpR).*exp(1j*pi*Kr*t_ref.^2).*exp(1j*pi*fc*t_ref);

s2 = s1.*conj(s_ref);
Fs2 = fftshift((fft(s2,2*TpR*fs+1,2)),2);

f = -fs/2:1/TpR/2:fs/2;
figure;
mesh(real(Fs2));
figure;
imagesc(f/2/Kr*c,tm,db(abs(Fs2)));
ylabel("慢时间/s");
xlabel("相对REF处距离/m");

%% GRT方法确定参数

pre_theta = -pi:pi/180:pi;
pre_R = 5:0.1:30;
map = zeros(length(pre_R),length(pre_theta));
test = zeros(1000,1801);

for i = 1:length(pre_R)
    for j = 1:length(pre_theta)
        test = zeros(1000,1801);
        for k = 1:PRF
            map(i,j) = map(i,j) + abs(Fs2(k,round((-pre_R(i)*cos(omega*k/PRF+pre_theta(j))+10*sqrt(3))*4*Kr*TpR/c+TpR*fs+1)));
            test(k,round((-pre_R(i)*cos(omega*k/PRF+pre_theta(j))+10*sqrt(3))*3+TpR*fs+1)) = 1;
        end
    end
end

figure;
mesh(pre_theta,pre_R,map);


% imagesc(f/2/Kr*c,tm,db(real(Fs2)));

%% 梯度GRT方法确定参数
DFs2 = abs(Fs2(1:PRF,2:2*TpR*fs+1))-abs(Fs2(1:PRF,1:2*TpR*fs));
strat_theta = pi/12;
strat_R = 17;

point_theta = strat_theta;
point_R = strat_R;

point_Rmap = zeros(1,200);
point_thetamap = zeros(1,200);

alpha = 0.000008;
figure;
for i = 1:400
    DD_DR = 0;
    DD_Dtheta = 0;
    for m = 1:PRF
        DD_DR = DD_DR + DFs2(m,round((-point_R*cos(omega*k/PRF+point_theta)+10*sqrt(3))*4*Kr*TpR/c+TpR*fs+1))*4*Kr*TpR/c*(-1)*cos(omega*m/PRF+point_theta);
        DD_Dtheta = DD_Dtheta + DFs2(m,round((-point_R*cos(omega*k/PRF+point_theta)+10*sqrt(3))*4*Kr*TpR/c+TpR*fs+1))*4*Kr*TpR/c*point_R*sin(omega*m/PRF+point_theta);
    end
    map_dr = map(2:length(pre_R),1:length(pre_theta))-map(1:length(pre_R)-1,1:length(pre_theta));
    map_dtheta = map(1:length(pre_R),2:length(pre_theta))-map(1:length(pre_R),1:length(pre_theta)-1);
    DD_DR_real = map_dr(round(point_R*10)-49,round(point_theta/pi*180)+181);
    DD_Dtheta_real = map_dtheta(round(point_R*10)-49,round(point_theta/pi*180)+181);

    quiver(point_theta,point_R,alpha*DD_Dtheta_real*10*pi/180,alpha*DD_DR_real);
    hold on;
    point_R = point_R + alpha*DD_DR_real;
    point_theta = point_theta + alpha*DD_Dtheta_real*10*pi/180;
    % point_R = point_R + alpha*DD_DR;
    % point_theta = point_theta + alpha*DD_Dtheta;
    point_Rmap(i) = point_R;
    point_thetamap(i) = point_theta; 
end

% plot(point_Rmap);
% plot(point_thetamap);
%% 确定频率

New_x = (0:PRF-1)/PRF;
New_z = ones(PRF,1)*(-TpR:1/fs:TpR);
New_y = sum(abs(Fs2).*New_z,2);

figure;
plot(New_x,New_y);

New_F = -PRF/2:PRF/2-1;
Fre_y = fftshift(fft(New_y));
plot(New_F,db(real(Fre_y)));