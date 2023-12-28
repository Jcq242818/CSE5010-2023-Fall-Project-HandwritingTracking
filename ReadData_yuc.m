clear ;
close all;
clc;
%addpath('fun');

chan_num_tar = 2;
chan_num_ref = 1;
total_duration=9;

f_c = 0.5e9;
f_s = 4e6;

CIT = 0.1; %跟模糊函数的数据量有关
CIT_region = CIT*f_s;

N_slide = 10;
T_slide = CIT / N_slide;%将 0.1s 的数据切割成N_slide 份，滑动的步长

max_dop = 1000;
step_dop = 1/CIT;
array_Doppler_frequency = [-max_dop:step_dop:max_dop];


chan_num_total = chan_num_tar+chan_num_ref;
for p = 1
%%
filename_data ='/home/lasso/Desktop/mmHandwriting_Jcq/data/1226/push_0.bin';


fid1=fopen(filename_data,'rb');  %rb - 二级制打开文件(1 Btye = 8 bits)
fseek(fid1,0,'eof');  %eof  打开文件末尾
fsize = ftell(fid1);  %获取文件长度

total_samplelen = fsize/4;  % /8 * 2， IQ分离

sample_1channel = 1996; %unit 大小  unit单位内IQ交叠 每一个unit只属于一个channel 单位应该是byte

samplelen = sample_1channel * (chan_num_total); %一个block 的长度= channel数量*单位长度
%sampleBytes = samplelen*4;
loopNum = floor(total_samplelen/samplelen); %循环次数 = 向下取整(总共长度/一个block长度)

samples_per_chan = sample_1channel * loopNum ;%一个channel采样点数:单位长度*循环次数(反复采样了多少次)，用于后面创建一个矩阵来存放数据
%DurationTruth = samples_per_chan/f_s;
fclose(fid1);


%% Read Calibration File
fid_cali=fopen(filename_data,'rb');%二进制byte打开文件
fseek(fid_cali,0,'bof');%bof 打开文件头
chan_seq_cali = zeros(chan_num_total,samples_per_chan);% 3xX 矩阵
for loopIdx = 1:loopNum
    [data,~]=fread(fid_cali,[2,sample_1channel*chan_num_total],'int16','l');% 竖着从左至右填满数组2xX，IQ分解?  int16 16字节数
    for idx_chan = 0:chan_num_total-1
        frame_start = sample_1channel*(0 + idx_chan)+1;
        frame_end = sample_1channel*(1 + idx_chan);
        chan_seq_cali(idx_chan+1,(loopIdx-1)*sample_1channel+1:loopIdx*sample_1channel) = ...%分三组写入数据
            data(1,frame_start:frame_end) + 1i*data(2,frame_start:frame_end);%雀食是IQ分解
    end
end
fclose(fid_cali);

%% plot 
offset = 1000;
ChannelIdx = 1;
EndSearch = samples_per_chan;
StartSearch = 1+offset;


% %提取信道
S_ref=chan_seq_cali(3,:);
S_tar1=chan_seq_cali(2,:);
% S_tar2=chan_seq_cali(3,:);
% 
% 
% 
%滑动窗口
array_start_time = 0:T_slide:total_duration-CIT;
A_TD = zeros(length(array_start_time),length(array_Doppler_frequency));
A_TD2 = zeros(length(array_start_time),length(array_Doppler_frequency));


%数据处理
ref_integer=zeros(length(array_start_time),CIT_region);
tar_integer=zeros(length(array_start_time),CIT_region);%准备fft数组表

% ref_integer2=zeros(length(array_start_time),CIT_region);
% tar_integer2=zeros(length(array_start_time),CIT_region);%准备fft数组表

for i= 1:length(array_start_time)-2
    ref_integer(i,:)=S_ref((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer(i,:)=S_tar1((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end                                                    %分类并填写fft数组表

% for i= 1:length(array_start_time)-2
%     ref_integer2(i,:)=S_ref((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
%     tar_integer2(i,:)=S_tar2((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
% end  

%FFT
% tar_integer=ClutterCancellation_Doppler(tar_integer,ref_integer);
% tar_integer2=ClutterCancellation_Doppler(tar_integer2,ref_integer);
for i= 1:length(array_start_time)-2
    final=fftshift(fft(tar_integer(i,:).*conj(ref_integer(i,:)),CIT_region));
    A_TD(i,:) = final(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
    %  final2=fftshift(fft(tar_integer2(i,:).*conj(ref_integer2(i,:)),CIT_region));
    % A_TD2(i,:) = final2(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
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
h1 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT);
xlim([array_start_time(1),array_start_time(end)]);
ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
set(gcf,'unit','centimeters','position',[5 3 30 15]);
set(get(gca,'XLabel'),'FontSize',22);
set(get(gca,'YLabel'),'FontSize',22);
colorbar;
title('1 and 2')
xlabel('Time (s)')
ylabel('Doppler frequency (Hz)')
colormap('jet');
caxis([thres_A_TRD,0]);
% saveas(gcf, 'E:\0617\5s_'+string(p)+'-1fft.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5a.jpg', 'jpg')
% 
% cfar2D = phased.CFARDetector2D('GuardBandSize',5,'TrainingBandSize',5,...
%   'ProbabilityFalseAlarm',0.5);
% resp=plot_A_DT;
% rngGrid=array_Doppler_frequency.';
% dopGrid=array_start_time.';
% rangeIndx(1)=20;
% rangeIndx(2)=180;
% dopplerIndx(1)=11;
% dopplerIndx(2)=800;
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

end

%% psd and real plot
figure
idxNum = 1e6;
for idx = 2:chan_num_total
    % subplot(chan_num_total,1,idx)
    pwelch(chan_seq_cali(idx,1:idxNum),[],[],[],f_s,'center','psd');hold on;
end

hold off;
figure
for idx = 2:chan_num_total
    % subplot(chan_num_total,1,idx)
    plot(real(chan_seq_cali(idx,1:idxNum)));hold on;
end
hold off;

legend;