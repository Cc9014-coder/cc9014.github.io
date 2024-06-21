%%通过解析信号和对数变换的方法求取群时延：
%这种方法首先计算信号的解析信号，然后对解析信号的幅度谱取对数，并通过差分来估计群时延。
%通过取对数，可以减少幅度接近零时的数值问题，并且对数变换有助于平滑处理。
%由于直接对幅度谱进行操作，而不是直接处理相位，所以这种方法可以避免因相位卷绕引起的问题
clc;clear;
path = 'F:\文献\相位解卷绕\1990数据';
List=dir('**/*.xls');%设置路径
filename={List.name}';
n=length(List);
fni=List(1).name;
fid=xlsread(fni);
t=fid(:,1);%时间序列
v=fid(:,4);%速度时间序列 单位cm/sec
s=fid(:,6);%位移时间序列 单位 cm
a=diff(v)./(t(2)-t(1));
a=[0;a];
% figure(1);plot(t,a);xlabel('时间 s');ylabel('加速度 cm/s^2');title('加速度时程曲线');
%%
%------------------------------计算持时------------------------------------
f=a.^2;%公式要求需去加速度的平方
for i=1:(length(a)-1)
    aint(i,:)=(t(2)-t(1))*(f(i+1)+f(i))/2;%梯形积分公式
end
arias=cumsum(aint);%能量持时的依次累加
arias2=arias(end);%能量持时的总和

T=t(arias>=0.001*arias2 & arias<=0.991*arias2);
I=[T(1),T(end)];
t_1=T(1);t_2=T(end);td=t_2-t_1;
a_=a(find(t==t_1):find(t==t_2),1);
N=length(a_);
%%
k=1;
nfft=2^nextpow2(N)*k;%取大于并最接近2的幂次方的长度
dt=td/nfft;df=1/td;
%--------------------------离散傅里叶变换----------------------------------
f = 0:df:(nfft/2-1)*df;%定义反应谱的离散频率向量
%对生成的高频地震波做傅里叶变换
AF = fft(a_,nfft);

%----------------------------计算群时延------------------------------------
% Calculate the Hilbert transform of h(t)
H_hilbert = hilbert(a_);

% Calculate the Fourier transform of the Hilbert transform
G = fft(H_hilbert,nfft);

% 因为MATLAB中FFT的变换矩阵不是一个酉矩阵(Unitary Matrix)，
% 该阵除以1/sqrt(N)就是个酉矩阵。故经过变换后对信号有放大作用，
% 所以要在fft处理后结果除以N/2来，修正此“放大”作用。

GB = 2*abs(G(1:nfft/2))/nfft;%幅值谱（真实傅里叶幅值谱）
Gan=angle(G(1:nfft/2));%主值相位谱 限定范围[-pi,pi]

% 对幅度谱进行对数变换前，加上一个小的常数以避免对数为负无穷
epsilon = 1e-6;
logH = log(G + epsilon) ;
%-------------------------------平滑处理-----------------------------------
% 初始化权重函数
weight_length = 3; % 权重函数的长度，必须是奇数
weights = zeros(1, weight_length);
half_length = round((weight_length + 1) / 2);
% 使用冒号表达式创建权重函数，注意使用方括号
weights(half_length:end) = 1:half_length;
weights(1:half_length) = half_length:-1:1;

% 计算权重函数的总和用于归一化
weights_sum = sum(weights);

% 初始化平滑后的分子和分母数组
smoothed_numerator = zeros(1, nfft/2);
smoothed_denominator = zeros(1, nfft/2);

% 应用三角加权函数进行平滑
for k = 1:nfft/2
    % 局部窗口，中心在 k
    window = (1:weight_length) - round(weight_length / 2) + k - 1;
    
    % 确保窗口内的索引不超出数组边界
    window(window < 1) = 1;
    window(window > nfft/2) = nfft/2;
    
    % 计算加权和
    weight_sum = weights * exp(1i * Gan(window));
    smoothed_numerator(k) = sum(weight_sum  .* GB(window));
    smoothed_denominator(k) = sum(abs(weight_sum) .^ 2);
end

% 避免除零错误
epsilon = 1e-6;
smoothed_denominator(smoothed_denominator == 0) = epsilon;
% 计算群时延
log_smoothed_numerator = log(smoothed_numerator + epsilon);
group_delay = -imag(diff(log_smoothed_numerator)) ./ (2 * pi * diff(f));

% % 计算相位导数，不需要展开相位
% group_delay = -imag(diff(logH)) ./ (2 * pi * diff(f) );

% 修正负数问题，群时延应该是非负的
group_delay(group_delay < 0) = NaN;
% 应用频率截止值
fm=10;
group_delay(f>fm)=NaN;
f(f>fm)=NaN;
figure(1);plot(group_delay,f,'o');ylabel('Frequency (Hz)');
xlabel('Group Delay (s)');
title('Boore方法');