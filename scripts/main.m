clc;
clear;
close all;


chan_num_tar = 2;
chan_num_ref = 1;
chan_num_total = chan_num_tar + chan_num_ref;
total_duration = 9;

% 
xT = 0.0;
yT = 0.0;
%sur1
xR1 = 2.6;
yR1 = 0.0;
%sur2
xR2 = 2*cos(pi/7);
yR2 = -2*sin(pi/7);
%
xTar = 2.3;
yTar = -2*sin(pi/7);

f_c = 0.5e9;
f_s = 1e6;

CIT = 0.1; %
CIT_region = CIT * f_s;


N_slide = 10;
T_slide = CIT / N_slide;% 0.1s N_slide 

max_dop = 1000;
step_dop = 1/CIT;
array_Doppler_frequency = -max_dop:step_dop:max_dop;


filename_data = "D:\Github\Passive-Handwriting-Tracking\data\1219_jcq\8_1.bin";


fid1=fopen(filename_data,'rb');  %rb - (1 Btye = 8 bits)
fseek(fid1,0,'eof');  %eof  
fsize = ftell(fid1);  %

total_samplelen = fsize/4;  % Byte/8 * 2 IQ

sample_1channel = 1996; %unit   unitIQ unitchannel byte

samplelen = sample_1channel * (chan_num_total); %block = channel*
%sampleBytes = samplelen*4;
loopNum = floor(total_samplelen/samplelen); % = (/block)

samples_per_chan = sample_1channel * loopNum ;%channel:*()
%DurationTruth = samples_per_chan/f_s;
fclose(fid1);

%% Read Calibration File
fid_cali=fopen(filename_data,'rb');%byte
fseek(fid_cali,0,'bof');%bof 
chan_seq_cali = zeros(chan_num_total,samples_per_chan);% 3xX 
for loopIdx = 1:loopNum
    [data,~]=fread(fid_cali,[2,sample_1channel*chan_num_total],'int16','l');% 2xXIQ int16 16
    for idx_chan = 0:chan_num_total-1
        frame_start = sample_1channel*(0 + idx_chan)+1;
        frame_end = sample_1channel*(1 + idx_chan);
        chan_seq_cali(idx_chan+1,(loopIdx-1)*sample_1channel+1:loopIdx*sample_1channel) = ...%
            data(1,frame_start:frame_end) + 1i*data(2,frame_start:frame_end);%IQIQ
    end
end
fclose(fid_cali);

%% 
offset = 1000;
ChannelIdx = 1;
EndSearch = samples_per_chan;
StartSearch = 1+offset;

%% 
data_count = 2184*1;
S_ref = chan_seq_cali(3,1:data_count);
S_tar = chan_seq_cali(1,1:data_count);
h_temp = xcorr(S_ref,S_tar);
[row_h,col_h] = size(h_temp);
col_h = col_h +1;
% N = floor(CIT*f_s)/100-1;
% N = 1000;
% col_max = find(max(abs(h_temp(col_h/2-N:col_h/2+N))) == abs(h_temp(col_h/2-N:col_h/2+N)));%
% array_sample_shift = col_max - N -1;

% figure(3)
% plot(abs(h_temp));

%
S_tar=chan_seq_cali(1,:);
S_ref1=chan_seq_cali(2,:);
S_ref2=chan_seq_cali(3,:);

N = floor(CIT*f_s)/1000-1;
col_max = find(max(abs(h_temp(col_h/2-N:col_h/2+N))) == abs(h_temp(col_h/2-N:col_h/2+N)));
array_sample_shift = col_max -N -1;
if array_sample_shift > 0
    data_cor1(1,:) = chan_seq_cali(1,1+array_sample_shift:end);
    data_cor1(2,:) = chan_seq_cali(2,1:end-array_sample_shift);
    data_cor2(1,:) = chan_seq_cali(1,1+array_sample_shift:end);
    data_cor2(2,:) = chan_seq_cali(3,1:end-array_sample_shift);
else
    data_cor1(1,:) = chan_seq_cali(1,1:end+array_sample_shift);
    data_cor1(2,:) = chan_seq_cali(2,1-array_sample_shift:end);
    data_cor2(1,:) = chan_seq_cali(1,1:end+array_sample_shift);
    data_cor2(2,:) = chan_seq_cali(3,1-array_sample_shift:end);
end


%
array_start_time = 0:T_slide:total_duration-CIT;
A_TD = zeros(length(array_start_time),length(array_Doppler_frequency));
ref_array = data_cor1(1,:); % ref_array
tar_array = data_cor1(2,:); % tar_array1

A_TD2 = zeros(length(array_start_time),length(array_Doppler_frequency));
ref_array2 = data_cor2(1,:); % ref_array2,ref_array
tar_array2 = data_cor2(2,:); % tar_array22

%
ref_integer=zeros(length(array_start_time),CIT_region);
tar_integer=zeros(length(array_start_time),CIT_region);%fft

ref_integer2=zeros(length(array_start_time),CIT_region);
tar_integer2=zeros(length(array_start_time),CIT_region);%fft

for i= 1:length(array_start_time)-2
    ref_integer(i,:)=ref_array((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer(i,:)=tar_array((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    ref_integer2(i,:)=ref_array2((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer2(i,:)=tar_array2((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end                                                    %fft

%FFT
% tar_integer=ClutterCancellation_Doppler(tar_integer,ref_integer);
% tar_integer2=ClutterCancellation_Doppler(tar_integer2,ref_integer);

for i= 1:length(array_start_time)-2
    final=fftshift(fft(tar_integer(i,:).*conj(ref_integer(i,:)),CIT_region));
    A_TD(i,:) = final(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
    final2=fftshift(fft(tar_integer2(i,:).*conj(ref_integer2(i,:)),CIT_region));
    A_TD2(i,:) = final2(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
end

%% CAF
% -Inf
A_TD = A_TD(:, 1:end-2);
A_TD2 = A_TD2(:, 1:end-2);
% 
A_TD(:, 101) = [];
A_TD2(:, 101) = [];

figure();
plot(A_TD);
plot(A_TD2);

%% CAF
% 
[maxValues, columnIndices] = max(abs(A_TD), [], 2);
[maxValues2, columnIndices2] = max(abs(A_TD2), [], 2);
% columnIndices ,
maxDF1 = (columnIndices - 101) * step_dop;
maxDF2 = (columnIndices2 - 101) * step_dop;

%% 

%%
fai_sur1 = zeros(1,length(array_start_time)-1);
fai_sur2 = zeros(1,length(array_start_time)-1);
fai_tx = zeros(1,length(array_start_time)-1);

fai_sur1(1) = atan((yTar)/ (xTar - xR1));
fai_sur2(1) = atan((yTar - yR2)/ (xTar - xR2));
fai_tx(1) = atan(yTar/xTar);

%
xtar = zeros(1,length(array_start_time)-1);
ytar = zeros(1,length(array_start_time)-1);
xtar(1) = xTar;
ytar(1) = yTar;
%
v_xy = zeros(1,2);
%
fd = zeros(1,2);
%F,ICCC
F = zeros(2,2);
%
fc = 60.48e9;
%
c = 3e8;
% 
A = eye(2);  % 
H = eye(2);  % 
Q = eye(2) * 1e-4;  % 
R = eye(2) * 1e-2;  % 
P = eye(2);  % 

vxArray = [0.1];
vyArray = [0.1];

for i = 1:1:length(array_start_time)-2
    fd = [maxDF1(i); maxDF2(i)];
    F = -2*fc/c*[cos((fai_sur1(i) - fai_tx(i))/2) * cos((fai_sur1(i) + fai_tx(i))/2) , ...
                 cos((fai_sur1(i) - fai_tx(i))/2) * sin((fai_sur1(i) + fai_tx(i))/2); ...
                 cos((fai_sur2(i) - fai_tx(i))/2) * cos((fai_sur2(i) + fai_tx(i))/2), ...
                 cos((fai_sur2(i) - fai_tx(i))/2) * sin((fai_sur2(i) + fai_tx(i))/2)];

    % 
    x_hat_minus = A * [vxArray(end); vyArray(end)];
    P_minus = A * P * A' + Q;
    % 
    K = P_minus * H' / (H * P_minus * H' + R);
    x_hat = x_hat_minus + K * (fd - H * x_hat_minus);
    P = (eye(2) - K * H) * P_minus;

    % 
    vMax = 0.5;
    for axis = [1, 2]
        if abs(x_hat(axis)) > vMax
            x_hat(axis) = sign(x_hat(1)) * vMax;
        end
    end

    vxArray = [vxArray, x_hat(1)];
    vyArray = [vyArray, x_hat(2)];

    xtar(i+1) = xtar(i) + x_hat(1) * T_slide;
    ytar(i+1) = ytar(i) + x_hat(2) * T_slide;

    % 
    fai_sur1(i+1) = atan((ytar(i+1) - yR1) / (xtar(i+1) - xR1));
    fai_sur2(i+1) = atan((ytar(i+1) - yR2) / (xtar(i+1) - xR2));
    fai_tx(i+1) = atan(ytar(i+1) / xtar(i+1));
end

% % 
% windowSize = 40;
% % vxArray
% vxArray = medfilt1(vxArray, windowSize);
% vyArray = medfilt1(vyArray, windowSize);

figure(2);
plot(vxArray, '-');
hold on;
plot(vyArray, '-');
hold off;
legend;

%% 
figure(3);
%  plot 
plot(xtar, ytar, 'x-', 'LineWidth', 2, 'MarkerSize', 4);

% 
xlabel('X');
ylabel('Y');
title('');

% 
title('');
xlabel('X');
ylabel('Y');

% 
grid on;


