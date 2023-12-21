%% 
clear ;
close all;
clc;


chan_num_tar = 2;
chan_num_ref = 1;
chan_num_total = chan_num_tar+chan_num_ref;
total_duration = 4;

f_c = 0.5e9;
f_s = 1e6;
%
c = 3e8;
%

CIT = 0.1; %
CIT_region = CIT*f_s;


N_slide = 10;
T_slide = CIT / N_slide;% 0.1s N_slide 

max_dop = 1000;
step_dop = 1/CIT;
array_Doppler_frequency = -max_dop:step_dop:max_dop;

filename_data = "D:\Github\Passive-Handwriting-Tracking\data\0413\0_small.bin";
% filename_data = "D:\Github\Passive-Handwriting-Tracking\data\0413\3.bin";
% filename_data = "D:\Github\Passive-Handwriting-Tracking\data\0413\7.bin";
% filename_data = "D:\Github\Passive-Handwriting-Tracking\data\0413\8.bin";

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

%% 
data_count = 2184*1;
S_ref = chan_seq_cali(3,1:data_count);
S_tar = chan_seq_cali(1,1:data_count);
h_temp = xcorr(S_ref,S_tar);
[row_h,col_h] = size(h_temp);
col_h = col_h +1;


%
S_tar=chan_seq_cali(1,:);
S_ref1=chan_seq_cali(2,:);
S_ref2=chan_seq_cali(3,:);

N = floor(CIT*f_s)/1000-1;
col_max = find(max(abs(h_temp(col_h/2-N:col_h/2+N))) == abs(h_temp(col_h/2-N:col_h/2+N)));
array_sample_shift = col_max -N -1;
if array_sample_shift>0
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
ref_array=data_cor1(1,:); % ref_array
tar_array=data_cor1(2,:); % tar_array1

A_TD2 = zeros(length(array_start_time),length(array_Doppler_frequency));
ref_array2=data_cor2(1,:); % ref_array2,ref_array
tar_array2=data_cor2(2,:); % tar_array22

%
ref_integer=zeros(length(array_start_time),CIT_region);
tar_integer=zeros(length(array_start_time),CIT_region);%fft

ref_integer2=zeros(length(array_start_time),CIT_region);
tar_integer2=zeros(length(array_start_time),CIT_region);%fft

for i= 1:length(array_start_time)-2
    ref_integer(i,:)=ref_array((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer(i,:)=tar_array((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end                                                    %fft

for i= 1:length(array_start_time)-2
    ref_integer2(i,:)=ref_array2((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer2(i,:)=tar_array2((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end  

%FFT
% tar_integer=ClutterCancellation_Doppler(tar_integer,ref_integer);
% tar_integer2=ClutterCancellation_Doppler(tar_integer2,ref_integer);

for i= 1:length(array_start_time)-2
    final=fftshift(fft(tar_integer(i,:).*conj(ref_integer(i,:)),CIT_region));
    A_TD(i,:) = final(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
    final2=fftshift(fft(tar_integer2(i,:).*conj(ref_integer2(i,:)),CIT_region));
    A_TD2(i,:) = final2(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
end



%% CAF

%%refsur1CAF
thres_A_TRD = -30;
fig1 = figure(1);
set(fig1,'position',[50,50,900,600]);
plot_A_DT = abs(A_TD');
plot_A_DT = mag2db(plot_A_DT/max(max(plot_A_DT)));
h1 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT);
xlim([array_start_time(1),array_start_time(end)]);
ylim([-100, 100]);
% y
set(gca, 'YDir', 'normal');
set(gcf,'unit','centimeters','position',[5 3 30 15]);
set(get(gca,'XLabel'),'FontSize',22);
set(get(gca,'YLabel'),'FontSize',22);
colorbar;
title('ref and sur1')
xlabel('Time (s)')
ylabel('Doppler frequency (Hz)')
colormap('jet');
clim([thres_A_TRD,0]);


%%refsur2CAF 
thres_A_TRD = -30;
fig2 = figure(2);
set(fig2,'position',[50,50,900,600]);
plot_A_DT2 = abs(A_TD2');
plot_A_DT2 = mag2db(plot_A_DT2/max(max(plot_A_DT2)));
h2 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT2);
xlim([array_start_time(1),array_start_time(end)]);
ylim([-100, 100]);
% y
set(gca, 'YDir', 'normal');
set(gcf,'unit','centimeters','position',[5 3 30 15]);
set(get(gca,'XLabel'),'FontSize',22);
set(get(gca,'YLabel'),'FontSize',22);
colorbar;
xlabel('Time (s)')
ylabel('Doppler frequency (Hz)')
title('ref and sur2')
colormap('jet');
clim([thres_A_TRD,0]);

%% CAF

% -Inf
A_TD = A_TD(1:end-2, :);
A_TD2 = A_TD2(1:end-2,:);
% 
A_TD(:, 101) = [];
A_TD2(:, 101) = [];
%%delete_index
delete_index = 51;
A_TD = A_TD(delete_index:end,:);
A_TD2 = A_TD2(delete_index:end,:);
%%


%refsur1CAF
thres_A_TRD = -30;
fig3 = figure(3);
set(fig3,'position',[50,50,900,600]);
plot_A_DT = abs(A_TD');
plot_A_DT = mag2db(plot_A_DT/max(max(plot_A_DT)));
h1 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT);
xlim([array_start_time(1),array_start_time(end)]);
ylim([-100, 100]);
% y
set(gca, 'YDir', 'normal');
set(gcf,'unit','centimeters','position',[5 3 30 15]);
set(get(gca,'XLabel'),'FontSize',22);
set(get(gca,'YLabel'),'FontSize',22);
colorbar;
title('ref and sur1')
xlabel('Time (s)')
ylabel('Doppler frequency (Hz)')
colormap('jet');
clim([thres_A_TRD,0]);


%refsur2CAF 
thres_A_TRD = -30;
fig4 = figure(4);
set(fig4,'position',[50,50,900,600]);
plot_A_DT2 = abs(A_TD2');
plot_A_DT2 = mag2db(plot_A_DT2/max(max(plot_A_DT2)));
h2 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT2);
xlim([array_start_time(1),array_start_time(end)]);
ylim([-100, 100]);
% y
set(gca, 'YDir', 'normal');
set(gcf,'unit','centimeters','position',[5 3 30 15]);
set(get(gca,'XLabel'),'FontSize',22);
set(get(gca,'YLabel'),'FontSize',22);
colorbar;
xlabel('Time (s)')
ylabel('Doppler frequency (Hz)')
title('ref and sur2')
colormap('jet');
clim([thres_A_TRD,0]);
close('all')

%% CAF
% 
[maxValues, columnIndices] = max(abs(A_TD), [], 2);
[maxValues2, columnIndices2] = max(abs(A_TD2), [], 2);
% columnIndices ,
maxDF1 = (columnIndices - 101) * step_dop;
maxDF2 = (columnIndices2 - 101) * step_dop;

%% 
% 对maxDF做积分
cumsumMaxDF = zeros(2,339);
cumsumMaxDF(1,:) = cumsum(maxDF1');
cumsumMaxDF(2,:) = cumsum(maxDF2');

% save('D:\Github\Passive-Handwriting-Tracking\data\train\data.mat', 'cumsumMaxDF');

% 对cumsum做卡尔曼平滑
cumsumMaxDF_KS = zeros(2,339);
cumsumMaxDF_KS(1,:) = KalmanSmoother(cumsumMaxDF(1,:));
cumsumMaxDF_KS(2,:) = KalmanSmoother(cumsumMaxDF(2,:));

figure(7);
plot(cumsumMaxDF(1,:), '-', 'LineWidth', 2, 'DisplayName', 'Channel 1');
hold on;
plot(cumsumMaxDF(2,:), '-', 'LineWidth', 2, 'DisplayName', 'Channel 2');
hold on;
plot(cumsumMaxDF_KS(1,:), '-', 'LineWidth', 2, 'DisplayName', 'Channel 1 KS');
hold on;
plot(cumsumMaxDF_KS(2,:), '-', 'LineWidth', 2, 'DisplayName', 'Channel 2 KS');
hold off;

legend('show');
grid on;


% 对cumsum做微分
maxDF1_KS = diff(cumsumMaxDF_KS(1,:));
maxDF2_KS = diff(cumsumMaxDF_KS(2,:));

figure(8);
plot(maxDF1, '-','LineWidth', 2, 'DisplayName', 'Channel 1 original');
hold on;
plot(maxDF2, '-','LineWidth', 2, 'DisplayName', 'Channel 2 original');
hold on;
plot(maxDF1_KS, '-','LineWidth', 2, 'DisplayName', 'Channel 1 KS undirect');
hold on;
plot(maxDF2_KS, '-','LineWidth', 2, 'DisplayName', 'Channel 2 KS undirect');
hold on;

% 对比直接做卡尔曼平滑
maxDF_KS = zeros(2,339);
maxDF_KS(1,:) = KalmanSmoother(maxDF1);
maxDF_KS(2,:) = KalmanSmoother(maxDF2);

plot(maxDF_KS(1,:), '-','LineWidth', 2, 'DisplayName', 'Channel 1 KS direct');
hold on;
plot(maxDF_KS(2,:), '-','LineWidth', 2, 'DisplayName', 'Channel 2 KS direct');
hold off;

ylim([-150, 150]);
yticks(-150:10:150);
legend('show');
grid on;

%% 
% 
xT = 0;
yT = 0;
%sur1
xR1 = 3;
yR1 = 0;
%sur2--sur2tx45
xR2 = sqrt(2);
yR2 = -sqrt(2);
%
xTar = 2;
yTar = -sqrt(2);


fai_sur1 = zeros(1,339);
fai_sur2 = zeros(1,339);
fai_tx = zeros(1,339);

fai_sur1(1) = atan((yTar)/ (xTar - xR1));
fai_sur2(1) = atan((yTar - yR2)/ (xTar - xR2));
fai_tx(1) = atan(yTar/xTar);

fai_sur1(1) = atan((yTar)/ (xTar - xR1));
fai_sur2(1) = atan((yTar - yR2)/ (xTar - xR2));
fai_tx(1) = atan(yTar/xTar);

%
xtar = zeros(1,339);
ytar = zeros(1,339);
xtar(1) = xTar;
ytar(1) = yTar;

%xy
v_xy = zeros(339,2);
vMax = 0.5;

%
fd = zeros(1,2);
%F,ICCC
F = zeros(2,2);


v_xy(1,:) = [0.1, 0.1];

for i = 1:1:338
    fd = [maxDF1_KS(i); maxDF2_KS(i)];
    F = -2*f_c/c*[cos((fai_sur1(i) - fai_tx(i))/2) * cos((fai_sur1(i) + fai_tx(i))/2) , ...
                 cos((fai_sur1(i) - fai_tx(i))/2) * sin((fai_sur1(i) + fai_tx(i))/2); ...
                 cos((fai_sur2(i) - fai_tx(i))/2) * cos((fai_sur2(i) + fai_tx(i))/2), ...
                 cos((fai_sur2(i) - fai_tx(i))/2) * sin((fai_sur2(i) + fai_tx(i))/2)];
    v_xy(i,:) = F \ fd;


    %
    % if (abs(v_xy(i,1)) > vMax) || (abs(v_xy(i,2)) > vMax)
    %     v_xy(i,:) = (v_xy(i-1,:) + v_xy(i-2,:) + v_xy(i-3,:))/3;
    %     v_xy(i,:) = sign(x_hat(1)) * vMax;
    % end

    if abs(v_xy(i,1)) > vMax
        v_xy(i,1) = sign(v_xy(i,1)) * vMax;
    end
    if abs(v_xy(i,2)) > vMax
        v_xy(i,2) = sign(v_xy(i,2)) * vMax;
    end

    xtar(i+1) = xtar(i) + v_xy(i,1) * T_slide;
    ytar(i+1) = ytar(i) + v_xy(i,2) * T_slide;

    %更新下一时刻的AOA角度以便下一次进行迭代
    fai_sur1(i+1) = atan((ytar(i+1) - yR1) / (xtar(i+1) - xR1));
    fai_sur2(i+1) = atan((ytar(i+1) - yR2) / (xtar(i+1) - xR2));
    fai_tx(i+1) = atan(ytar(i+1) / xtar(i+1));
end


figure(9);
plot(v_xy(:,1),'LineWidth', 2, 'DisplayName', 'vx');
hold on;
plot(v_xy(:,2),'LineWidth', 2, 'DisplayName', 'vy');
hold off;
legend('show');

%% 
% 
figure(10);
plot(xtar, ytar, '-','Color', [1, 0.5, 0],'LineWidth', 3,'DisplayName', 'Estimated Trajectory');
% 
hold on;
plot(xtar(1), ytar(1), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'Initial Point');
plot(xtar(end), ytar(end), 'x', 'MarkerSize', 6, 'MarkerEdgeColor', 'b',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'End Point');
hold off;

legend('show','Location', 'northeast');
axis equal;
% 
grid on;

