%%ͨ�������źźͶ����任�ķ�����ȡȺʱ�ӣ�
%���ַ������ȼ����źŵĽ����źţ�Ȼ��Խ����źŵķ�����ȡ��������ͨ�����������Ⱥʱ�ӡ�
%ͨ��ȡ���������Լ��ٷ��Ƚӽ���ʱ����ֵ���⣬���Ҷ����任������ƽ������
%����ֱ�ӶԷ����׽��в�����������ֱ�Ӵ�����λ���������ַ������Ա�������λ�������������
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
% figure(1);plot(t,a);xlabel('ʱ�� s');ylabel('���ٶ� cm/s^2');title('���ٶ�ʱ������');
%%
%------------------------------�����ʱ------------------------------------
f=a.^2;%��ʽҪ����ȥ���ٶȵ�ƽ��
for i=1:(length(a)-1)
    aint(i,:)=(t(2)-t(1))*(f(i+1)+f(i))/2;%���λ��ֹ�ʽ
end
arias=cumsum(aint);%������ʱ�������ۼ�
arias2=arias(end);%������ʱ���ܺ�

T=t(arias>=0.001*arias2 & arias<=0.991*arias2);
I=[T(1),T(end)];
t_1=T(1);t_2=T(end);td=t_2-t_1;
a_=a(find(t==t_1):find(t==t_2),1);
N=length(a_);
%%
k=1;
nfft=2^nextpow2(N)*k;%ȡ���ڲ���ӽ�2���ݴη��ĳ���
dt=td/nfft;df=1/td;
%--------------------------��ɢ����Ҷ�任----------------------------------
f = 0:df:(nfft/2-1)*df;%���巴Ӧ�׵���ɢƵ������
%�����ɵĸ�Ƶ����������Ҷ�任
AF = fft(a_,nfft);

%----------------------------����Ⱥʱ��------------------------------------
% Calculate the Hilbert transform of h(t)
H_hilbert = hilbert(a_);

% Calculate the Fourier transform of the Hilbert transform
G = fft(H_hilbert,nfft);

% ��ΪMATLAB��FFT�ı任������һ���Ͼ���(Unitary Matrix)��
% �������1/sqrt(N)���Ǹ��Ͼ��󡣹ʾ����任����ź��зŴ����ã�
% ����Ҫ��fft�����������N/2���������ˡ��Ŵ����á�

GB = 2*abs(G(1:nfft/2))/nfft;%��ֵ�ף���ʵ����Ҷ��ֵ�ף�
Gan=angle(G(1:nfft/2));%��ֵ��λ�� �޶���Χ[-pi,pi]

% �Է����׽��ж����任ǰ������һ��С�ĳ����Ա������Ϊ������
epsilon = 1e-6;
logH = log(G + epsilon) ;
%-------------------------------ƽ������-----------------------------------
% ��ʼ��Ȩ�غ���
weight_length = 3; % Ȩ�غ����ĳ��ȣ�����������
weights = zeros(1, weight_length);
half_length = round((weight_length + 1) / 2);
% ʹ��ð�ű��ʽ����Ȩ�غ�����ע��ʹ�÷�����
weights(half_length:end) = 1:half_length;
weights(1:half_length) = half_length:-1:1;

% ����Ȩ�غ������ܺ����ڹ�һ��
weights_sum = sum(weights);

% ��ʼ��ƽ����ķ��Ӻͷ�ĸ����
smoothed_numerator = zeros(1, nfft/2);
smoothed_denominator = zeros(1, nfft/2);

% Ӧ�����Ǽ�Ȩ��������ƽ��
for k = 1:nfft/2
    % �ֲ����ڣ������� k
    window = (1:weight_length) - round(weight_length / 2) + k - 1;
    
    % ȷ�������ڵ���������������߽�
    window(window < 1) = 1;
    window(window > nfft/2) = nfft/2;
    
    % �����Ȩ��
    weight_sum = weights * exp(1i * Gan(window));
    smoothed_numerator(k) = sum(weight_sum  .* GB(window));
    smoothed_denominator(k) = sum(abs(weight_sum) .^ 2);
end

% ����������
epsilon = 1e-6;
smoothed_denominator(smoothed_denominator == 0) = epsilon;
% ����Ⱥʱ��
log_smoothed_numerator = log(smoothed_numerator + epsilon);
group_delay = -imag(diff(log_smoothed_numerator)) ./ (2 * pi * diff(f));

% % ������λ����������Ҫչ����λ
% group_delay = -imag(diff(logH)) ./ (2 * pi * diff(f) );

% �����������⣬Ⱥʱ��Ӧ���ǷǸ���
group_delay(group_delay < 0) = NaN;
% Ӧ��Ƶ�ʽ�ֵֹ
fm=10;
group_delay(f>fm)=NaN;
f(f>fm)=NaN;
figure(1);plot(group_delay,f,'o');ylabel('Frequency (Hz)');
xlabel('Group Delay (s)');
title('Boore����');