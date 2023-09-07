% 生成Chirp信号并绘制幅度相位谱
%% 参数
Tp = 20e-6;
B = 150e6;
Kr = B/Tp;
fc = 9e9;
fs = 200e6;

%% 产生基带Chirp信号
t = (-Tp*fs/2:Tp*fs/2)/fs;
BaseChirpSignal = exp(1j*pi*Kr.*t.^2);
BaseChirpSignal_fft = fftshift(fft(BaseChirpSignal));
f = (-Tp*fs/2:Tp*fs/2)/Tp;
figure;
plot(f,abs(BaseChirpSignal_fft));
xlabel("频率/Hz");
title("Chirp幅度谱");
figure;
plot(t,angle(BaseChirpSignal_fft));
xlabel("频率/Hz");
title("Chirp相位谱");

