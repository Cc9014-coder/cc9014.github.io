%%基于解析信号FFT G(f)的方法：
%这种方法直接使用原始信号的FFT结果和其共轭来计算群时延，通过分子中的实部和虚部乘积之和，
%以及分母中的幅度平方来实现。
%这种方法在计算时不直接处理相位，而是通过实部和虚部的乘积来间接利用相位信息，
%因此也可以避免相位卷绕的问题。
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
%%
td=t(end)-t(1);
N=length(a);
k=1;
nfft=2^nextpow2(N)*k;%取大于并最接近2的幂次方的长度
dt=td/nfft;df=1/td;
%--------------------------离散傅里叶变换----------------------------------
f = 0:df:(nfft/2-1)*df;%定义反应谱的离散频率向量
%对生成的高频地震波做傅里叶变换
AF = fft(a,nfft);
% Calculate the Hilbert transform of h(t)
H_hilbert = hilbert(a);
% Calculate the Fourier transform of the Hilbert transform
G = fft(H_hilbert,nfft);
%----------------------------计算群时延------------------------------------

%计算群时延
tf = zeros(1, nfft/2);
numerator = real(AF) .* real(conj(G)) + imag(AF) .* imag(conj(G));

denominator = (AF .* conj(AF)) ;

%-------------------------------平滑处理-----------------------------------
% 定义权重函数的长度，必须是奇数
weight_length = 3; % 例如，权重函数的长度为 9

% 创建三角加权函数（权重数组）
weights = zeros(1, weight_length);
half_length = round((weight_length + 1) / 2);
weights(half_length:end) = 1:half_length;
weights(1:half_length) = half_length:-1:1;

% 现在 weights 包含了一个简单的三角加权函数

% 接下来，您可以使用这个权重函数对分子和分母进行平滑处理
% 假设 AF 和 G 是已经计算好的 FFT 结果
% 我们只对正频率部分进行平滑处理

% 初始化平滑后的分子和分母数组
smoothed_numerator = zeros(1, nfft/2);
smoothed_denominator = zeros(1, nfft/2);

% 应用三角加权函数进行平滑
for k = 1:nfft/2
    % 局部窗口，中心在 k
    window = (1:weight_length) - round(weight_length / 2) + k;
    
    % 确保窗口内的索引不超出数组边界
    window(window < 1) = 1;
    window(window > nfft/2) = nfft/2;
    
    % 计算加权和
    weight_sum = weights * exp(1i * angle(AF(window)));
    smoothed_numerator(k) = sum(weight_sum .* exp(1i * angle(G(window))));
    smoothed_denominator(k) = sum(abs(weight_sum) .^ 2);
end

% 确保分母不为零
smoothed_denominator(smoothed_denominator == 0) = 1; % 使用一个小的值避免除零错误

tf = - 2 * pi * smoothed_numerator ./ smoothed_denominator;

% % 修正负数问题，群时延应该是非负的
% tf(tf < 0) = NaN;
% % 应用频率截止值
% fm=60;
% tf(f>fm)=NaN;
% % f(f>fm)=NaN;

figure(3);plot(tf(1:nfft/2),f,'o');ylabel('Frequency (Hz)');
xlabel('Group Delay (s)');
title('Boore方法');
