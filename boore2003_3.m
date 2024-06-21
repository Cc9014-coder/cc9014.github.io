%%���ڽ����ź�FFT G(f)�ķ�����
%���ַ���ֱ��ʹ��ԭʼ�źŵ�FFT������乲��������Ⱥʱ�ӣ�ͨ�������е�ʵ�����鲿�˻�֮�ͣ�
%�Լ���ĸ�еķ���ƽ����ʵ�֡�
%���ַ����ڼ���ʱ��ֱ�Ӵ�����λ������ͨ��ʵ�����鲿�ĳ˻������������λ��Ϣ��
%���Ҳ���Ա�����λ���Ƶ����⡣
clc;clear;
path = 'F:\����\��λ�����\1990����';
List=dir('**/*.xls');%����·��
filename={List.name}';
n=length(List);
fni=List(1).name;
fid=xlsread(fni);
t=fid(:,1);%ʱ������
v=fid(:,4);%�ٶ�ʱ������ ��λcm/sec
s=fid(:,6);%λ��ʱ������ ��λ cm
a=diff(v)./(t(2)-t(1));
a=[0;a];
%%
td=t(end)-t(1);
N=length(a);
k=1;
nfft=2^nextpow2(N)*k;%ȡ���ڲ���ӽ�2���ݴη��ĳ���
dt=td/nfft;df=1/td;
%--------------------------��ɢ����Ҷ�任----------------------------------
f = 0:df:(nfft/2-1)*df;%���巴Ӧ�׵���ɢƵ������
%�����ɵĸ�Ƶ����������Ҷ�任
AF = fft(a,nfft);
% Calculate the Hilbert transform of h(t)
H_hilbert = hilbert(a);
% Calculate the Fourier transform of the Hilbert transform
G = fft(H_hilbert,nfft);
%----------------------------����Ⱥʱ��------------------------------------

%����Ⱥʱ��
tf = zeros(1, nfft/2);
numerator = real(AF) .* real(conj(G)) + imag(AF) .* imag(conj(G));

denominator = (AF .* conj(AF)) ;

%-------------------------------ƽ������-----------------------------------
% ����Ȩ�غ����ĳ��ȣ�����������
weight_length = 3; % ���磬Ȩ�غ����ĳ���Ϊ 9

% �������Ǽ�Ȩ������Ȩ�����飩
weights = zeros(1, weight_length);
half_length = round((weight_length + 1) / 2);
weights(half_length:end) = 1:half_length;
weights(1:half_length) = half_length:-1:1;

% ���� weights ������һ���򵥵����Ǽ�Ȩ����

% ��������������ʹ�����Ȩ�غ����Է��Ӻͷ�ĸ����ƽ������
% ���� AF �� G ���Ѿ�����õ� FFT ���
% ����ֻ����Ƶ�ʲ��ֽ���ƽ������

% ��ʼ��ƽ����ķ��Ӻͷ�ĸ����
smoothed_numerator = zeros(1, nfft/2);
smoothed_denominator = zeros(1, nfft/2);

% Ӧ�����Ǽ�Ȩ��������ƽ��
for k = 1:nfft/2
    % �ֲ����ڣ������� k
    window = (1:weight_length) - round(weight_length / 2) + k;
    
    % ȷ�������ڵ���������������߽�
    window(window < 1) = 1;
    window(window > nfft/2) = nfft/2;
    
    % �����Ȩ��
    weight_sum = weights * exp(1i * angle(AF(window)));
    smoothed_numerator(k) = sum(weight_sum .* exp(1i * angle(G(window))));
    smoothed_denominator(k) = sum(abs(weight_sum) .^ 2);
end

% ȷ����ĸ��Ϊ��
smoothed_denominator(smoothed_denominator == 0) = 1; % ʹ��һ��С��ֵ����������

tf = - 2 * pi * smoothed_numerator ./ smoothed_denominator;

% % �����������⣬Ⱥʱ��Ӧ���ǷǸ���
% tf(tf < 0) = NaN;
% % Ӧ��Ƶ�ʽ�ֵֹ
% fm=60;
% tf(f>fm)=NaN;
% % f(f>fm)=NaN;

figure(3);plot(tf(1:nfft/2),f,'o');ylabel('Frequency (Hz)');
xlabel('Group Delay (s)');
title('Boore����');
