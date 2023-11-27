%% ��������
clear ;
close all;
clc;
%addpath('fun');



chan_num_tar = 2;
chan_num_ref = 1;
total_duration = 9;

f_c = 0.5e9;
f_s = 1e6;

CIT = 0.1; %��ģ���������������й�
CIT_region = CIT*f_s;


N_slide = 10;
T_slide = CIT / N_slide;%�� 0.1s �������и��N_slide �ݣ������Ĳ���

max_dop = 1000;
step_dop = 1/CIT;
array_Doppler_frequency = -max_dop:step_dop:max_dop;


chan_num_total = chan_num_tar+chan_num_ref;


filename_data = 'E:\Desktop\Project\data\����6��\6_1.bin';
% filename_data='E:\0618_b\0618_b\b10_'+string(p)+'.bin';

fid1=fopen(filename_data,'rb');  %rb - �����ƴ��ļ�(1 Btye = 8 bits)
fseek(fid1,0,'eof');  %eof  ���ļ�ĩβ
fsize = ftell(fid1);  %��ȡ�ļ�����

total_samplelen = fsize/4;  % Byte/8 * 2�� IQ����

sample_1channel = 1996; %unit ��С  unit��λ��IQ���� ÿһ��unitֻ����һ��channel ��λӦ����byte

samplelen = sample_1channel * (chan_num_total); %һ��block �ĳ���= channel����*��λ����
%sampleBytes = samplelen*4;
loopNum = floor(total_samplelen/samplelen); %ѭ������ = ����ȡ��(�ܹ�����/һ��block����)

samples_per_chan = sample_1channel * loopNum ;%һ��channel��������:��λ����*ѭ������(���������˶��ٴ�)�����ں��洴��һ���������������
%DurationTruth = samples_per_chan/f_s;
fclose(fid1);


%% Read Calibration File
fid_cali=fopen(filename_data,'rb');%������byte���ļ�
fseek(fid_cali,0,'bof');%bof ���ļ�ͷ
chan_seq_cali = zeros(chan_num_total,samples_per_chan);% 3xX ����
for loopIdx = 1:loopNum
    [data,~]=fread(fid_cali,[2,sample_1channel*chan_num_total],'int16','l');% ���Ŵ���������������2xX��IQ�ֽ� int16 16�ֽ���
    for idx_chan = 0:chan_num_total-1
        frame_start = sample_1channel*(0 + idx_chan)+1;
        frame_end = sample_1channel*(1 + idx_chan);
        chan_seq_cali(idx_chan+1,(loopIdx-1)*sample_1channel+1:loopIdx*sample_1channel) = ...%������д������
            data(1,frame_start:frame_end) + 1i*data(2,frame_start:frame_end);%ȸʳ��IQ�ֽ⣬��֮ǰ��IQ�ֽ����
    end
end
fclose(fid_cali);

%% ��ͼ�������ź��Ƿ���ȷ
offset = 1000;
ChannelIdx = 1;
EndSearch = samples_per_chan;
StartSearch = 1+offset;
chan_num_total = 3;

% figure(1)
% for idxChannel = 1:chan_num_total
%     subplot(chan_num_total,1,idxChannel)
%     plot(real(chan_seq_cali(idxChannel,StartSearch:EndSearch)));
%     hold on;
%     plot(imag(chan_seq_cali(idxChannel,StartSearch:EndSearch)));
%     hold off;
% end
% sgtitle('ʵ�ʽ��յ����ź�');% ����
% figure(2)
% 
% for idxChannel = 1:chan_num_total
%     subplot(chan_num_total,1,idxChannel)
%     pwelch(chan_seq_cali(idxChannel,StartSearch:EndSearch),[],[],[],f_s,'centered','psd');
% end
% sgtitle('psd');%����
% 
%% ���ݴ���
data_count = 2184*1;
S_ref = chan_seq_cali(3,1:data_count);
S_tar = chan_seq_cali(1,1:data_count);
h_temp = xcorr(S_ref,S_tar);
[row_h,col_h] = size(h_temp);
col_h = col_h +1;
% N = floor(CIT*f_s)/100-1;
% N = 1000;
% col_max = find(max(abs(h_temp(col_h/2-N:col_h/2+N))) == abs(h_temp(col_h/2-N:col_h/2+N)));%ģ����������Ԥ����һ��
% array_sample_shift = col_max - N -1;

% figure(3)
% plot(abs(h_temp));

%��ȡ�ŵ�
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


%��������
array_start_time = 0:T_slide:total_duration-CIT;
% A_TD = zeros(length(array_start_time),length(array_Doppler_frequency));
A_TD = zeros(length(array_start_time),length(-60:10:60));
ref_array=data_cor1(1,:); % ref_array����װ�Ųο�����������
tar_array=data_cor1(2,:); % tar_array����װ�ż���1����������

% A_TD2 = zeros(length(array_start_time),length(array_Doppler_frequency));
A_TD2 = zeros(length(array_start_time),length(-60:10:60));
ref_array2=data_cor2(1,:); % ref_array2����װ�Ųο�����������,��ref_array���������һģһ��
tar_array2=data_cor2(2,:); % tar_array2����װ�ż���2����������

%���ݴ���
ref_integer=zeros(length(array_start_time),CIT_region);
tar_integer=zeros(length(array_start_time),CIT_region);%׼��fft�����

ref_integer2=zeros(length(array_start_time),CIT_region);
tar_integer2=zeros(length(array_start_time),CIT_region);%׼��fft�����

for i= 1:length(array_start_time)-2
    ref_integer(i,:)=ref_array((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer(i,:)=tar_array((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end                                                    %���ಢ��дfft�����

for i= 1:length(array_start_time)-2
    ref_integer2(i,:)=ref_array2((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer2(i,:)=tar_array2((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end  

%FFT
% tar_integer=ClutterCancellation_Doppler(tar_integer,ref_integer);
% tar_integer2=ClutterCancellation_Doppler(tar_integer2,ref_integer);

for i= 1:length(array_start_time)-2
    final=fftshift(fft(tar_integer(i,:).*conj(ref_integer(i,:)),CIT_region));
    % A_TD(i,:) = final(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
    A_TD(i,:) = final(CIT_region/2+1-6:CIT_region/2+1+6);
     final2=fftshift(fft(tar_integer2(i,:).*conj(ref_integer2(i,:)),CIT_region));
    % A_TD2(i,:) = final2(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
    A_TD2(i,:) = final2(CIT_region/2+1-6:CIT_region/2+1+6);
end



%PLOT
% 
% A_TD2(95,:)=0;
% A_TD(95,:)=0;


thres_A_TRD = -30;
fig1 = figure(1);
set(fig1,'position',[50,50,900,600]);
plot_A_DT = abs(A_TD');
plot_A_DT = mag2db(plot_A_DT/max(max(plot_A_DT)));
% h1 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT);
h1 = imagesc(array_start_time,-60:10:60,plot_A_DT);
xlim([array_start_time(1),array_start_time(end)]);
% ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
ylim([-60,60]);
% set(gcf,'unit','centimeters','position',[5 3 30 15]);
% set(get(gca,'XLabel'),'FontSize',22);
% set(get(gca,'YLabel'),'FontSize',22);
colorbar;
% title('1 and 2')
xlabel('Time (s)')
ylabel('Doppler frequency (Hz)')
colormap('jet');
caxis([thres_A_TRD,0]);
% saveas(gcf, 'E:\0617\5s_'+string(p)+'-1fft.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5a.jpg', 'jpg')

% cfar2D = phased.CFARDetector2D('GuardBandSize',5,'TrainingBandSize',5,...
%   'ProbabilityFalseAlarm',0.5);
% resp=plot_A_DT;
% rngGrid=array_Doppler_frequency.';
% dopGrid=array_start_time.';
% rangeIndx(1)=229;
% rangeIndx(2)=287;
% dopplerIndx(1)=11;
% dopplerIndx(2)=150;
% [columnInds,rowInds] = meshgrid(dopplerIndx(1):dopplerIndx(2),...
%   rangeIndx(1):rangeIndx(2));
% CUTIdx = [rowInds(:) columnInds(:)]';
% detections = cfar2D(resp,CUTIdx);
% fig2 = figure(2);
% helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections)
% %saveas(gcf, 'E:\0617\5s_'+string(p)+'-2cfar.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5c.jpg', 'jpg')
% 
% 
% 
% 
% 
% 
% 
% 
% thres_A_TRD = -30;
% fig3 = figure(3);
% set(fig3,'position',[50,50,900,600]);
% plot_A_DT2 = abs(A_TD2');
% plot_A_DT2 = mag2db(plot_A_DT2/max(max(plot_A_DT2)));
% h2 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT2);
% xlim([array_start_time(1),array_start_time(end)]);
% ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
% set(gcf,'unit','centimeters','position',[5 3 30 15]);
% set(get(gca,'XLabel'),'FontSize',22);
% set(get(gca,'YLabel'),'FontSize',22);
% colorbar;
% xlabel('Time (s)')
% ylabel('Doppler frequency (Hz)')
% title('1 and 3')
% colormap('jet');
% caxis([thres_A_TRD,0]);
%saveas(gcf, 'E:\0617\5s_'+string(p)+'-3fft.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5b.jpg', 'jpg')
% 
% 

% cfar2D = phased.CFARDetector2D('GuardBandSize',5,'TrainingBandSize',5,...
%   'ProbabilityFalseAlarm',0.5);
% resp=plot_A_DT2;
% rngGrid=array_Doppler_frequency.';
% dopGrid=array_start_time.';
% rangeIndx(1)=229;
% rangeIndx(2)=287;
% dopplerIndx(1)=11;
% dopplerIndx(2)=150;
% [columnInds,rowInds] = meshgrid(dopplerIndx(1):dopplerIndx(2),...
%   rangeIndx(1):rangeIndx(2));
% CUTIdx = [rowInds(:) columnInds(:)]';
% detections = cfar2D(resp,CUTIdx);
% fig4 = figure(4);
% helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections)
% saveas(gcf, 'E:\0617\5s_'+string(p)+'-4cfar.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5f.jpg', 'jpg')
% close all;



