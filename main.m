%% 参数定义
clear ;
close all;
clc;
%addpath('fun');



chan_num_tar = 2;
chan_num_ref = 1;
total_duration = 4;

f_c = 0.5e9;
f_s = 1e6;

CIT = 0.1; %跟模糊函数的数据量有关
CIT_region = CIT*f_s;


N_slide = 10;
T_slide = CIT / N_slide;%将 0.1s 的数据切割成N_slide 份，滑动的步长

max_dop = 1000;
step_dop = 1/CIT;
array_Doppler_frequency = -max_dop:step_dop:max_dop;


chan_num_total = chan_num_tar+chan_num_ref;


filename_data = "E:\Desktop\Project\data\0413\3_small.bin";


fid1=fopen(filename_data,'rb');  %rb - 二级制打开文件(1 Btye = 8 bits)
fseek(fid1,0,'eof');  %eof  打开文件末尾
fsize = ftell(fid1);  %获取文件长度

total_samplelen = fsize/4;  % Byte/8 * 2， IQ分离

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
    [data,~]=fread(fid_cali,[2,sample_1channel*chan_num_total],'int16','l');% 竖着从左至右填满数组2xX，IQ分解 int16 16字节数
    for idx_chan = 0:chan_num_total-1
        frame_start = sample_1channel*(0 + idx_chan)+1;
        frame_end = sample_1channel*(1 + idx_chan);
        chan_seq_cali(idx_chan+1,(loopIdx-1)*sample_1channel+1:loopIdx*sample_1channel) = ...%分三组写入数据
            data(1,frame_start:frame_end) + 1i*data(2,frame_start:frame_end);%雀食是IQ分解，把之前的IQ分解组合
    end
end
fclose(fid_cali);

%% 绘图检查接收信号是否正确
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
% sgtitle('实际接收到的信号');% 整体
% figure(2)
% 
% for idxChannel = 1:chan_num_total
%     subplot(chan_num_total,1,idxChannel)
%     pwelch(chan_seq_cali(idxChannel,StartSearch:EndSearch),[],[],[],f_s,'centered','psd');
% end
% sgtitle('psd');%整体
% 
%% 数据处理
data_count = 2184*1;
S_ref = chan_seq_cali(3,1:data_count);
S_tar = chan_seq_cali(1,1:data_count);
h_temp = xcorr(S_ref,S_tar);
[row_h,col_h] = size(h_temp);
col_h = col_h +1;
% N = floor(CIT*f_s)/100-1;
% N = 1000;
% col_max = find(max(abs(h_temp(col_h/2-N:col_h/2+N))) == abs(h_temp(col_h/2-N:col_h/2+N)));%模糊函数函数预处理一组
% array_sample_shift = col_max - N -1;

% figure(3)
% plot(abs(h_temp));

%提取信道
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


%滑动窗口
array_start_time = 0:T_slide:total_duration-CIT;
A_TD = zeros(length(array_start_time),length(array_Doppler_frequency));
ref_array=data_cor1(1,:); % ref_array里面装着参考对齐后的数据
tar_array=data_cor1(2,:); % tar_array里面装着监视1对齐后的数据

A_TD2 = zeros(length(array_start_time),length(array_Doppler_frequency));
ref_array2=data_cor2(1,:); % ref_array2里面装着参考对齐后的数据,和ref_array里面的数据一模一样
tar_array2=data_cor2(2,:); % tar_array2里面装着监视2对齐后的数据

%数据处理
ref_integer=zeros(length(array_start_time),CIT_region);
tar_integer=zeros(length(array_start_time),CIT_region);%准备fft数组表

ref_integer2=zeros(length(array_start_time),CIT_region);
tar_integer2=zeros(length(array_start_time),CIT_region);%准备fft数组表

for i= 1:length(array_start_time)-2
    ref_integer(i,:)=ref_array((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer(i,:)=tar_array((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end                                                    %分类并填写fft数组表

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



%% 绘制CAF图
% 
% A_TD2(95,:)=0;
% A_TD(95,:)=0;

%%画ref和sur1的CAF
thres_A_TRD = -30;
fig1 = figure(1);
set(fig1,'position',[50,50,900,600]);
plot_A_DT = abs(A_TD');
plot_A_DT = mag2db(plot_A_DT/max(max(plot_A_DT)));
h1 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT);
xlim([array_start_time(1),array_start_time(end)]);
ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
% 反转y轴方向
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

% saveas(gcf, 'E:\0617\5s_'+string(p)+'-1fft.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5a.jpg', 'jpg')



%%画ref和sur2的CAF 
thres_A_TRD = -30;
fig2 = figure(2);
set(fig2,'position',[50,50,900,600]);
plot_A_DT2 = abs(A_TD2');
plot_A_DT2 = mag2db(plot_A_DT2/max(max(plot_A_DT2)));
h2 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT2);
xlim([array_start_time(1),array_start_time(end)]);
ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
% 反转y轴方向
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


%saveas(gcf, 'E:\0617\5s_'+string(p)+'-3fft.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5b.jpg', 'jpg')
% close all;

%% 提取CAF结果每一列的最大值并进行迭代

% 最后两列数据是-Inf，删除
A_TD = A_TD(1:end-2, :);
A_TD2 = A_TD2(1:end-2,:);
% 删除零频附近的值
A_TD(:, 101) = 0;
A_TD2(:, 101) = 0;

%%在画图对比一下

%%画ref和sur1的CAF
% thres_A_TRD = -30;
% fig3 = figure(3);
% set(fig3,'position',[50,50,900,600]);
% plot_A_DT = abs(A_TD');
% plot_A_DT = mag2db(plot_A_DT/max(max(plot_A_DT)));
% h1 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT);
% xlim([array_start_time(1),array_start_time(end)]);
% ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
% % 反转y轴方向
% set(gca, 'YDir', 'normal');
% set(gcf,'unit','centimeters','position',[5 3 30 15]);
% set(get(gca,'XLabel'),'FontSize',22);
% set(get(gca,'YLabel'),'FontSize',22);
% colorbar;
% title('ref and sur1')
% xlabel('Time (s)')
% ylabel('Doppler frequency (Hz)')
% colormap('jet');
% clim([thres_A_TRD,0]);

% saveas(gcf, 'E:\0617\5s_'+string(p)+'-1fft.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5a.jpg', 'jpg')



%%画ref和sur2的CAF 
% thres_A_TRD = -30;
% fig4 = figure(4);
% set(fig4,'position',[50,50,900,600]);
% plot_A_DT2 = abs(A_TD2');
% plot_A_DT2 = mag2db(plot_A_DT2/max(max(plot_A_DT2)));
% h2 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT2);
% xlim([array_start_time(1),array_start_time(end)]);
% ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
% % 反转y轴方向
% set(gca, 'YDir', 'normal');
% set(gcf,'unit','centimeters','position',[5 3 30 15]);
% set(get(gca,'XLabel'),'FontSize',22);
% set(get(gca,'YLabel'),'FontSize',22);
% colorbar;
% xlabel('Time (s)')
% ylabel('Doppler frequency (Hz)')
% title('ref and sur2')
% colormap('jet');
% clim([thres_A_TRD,0]);


%saveas(gcf, 'E:\0617\5s_'+string(p)+'-3fft.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5b.jpg', 'jpg')
% close all;

%% 提取两个CAF矩阵每一行的最大值所在的列
% 找到每一行中的最大值及其对应的列索引
[maxValues, columnIndices] = max(abs(A_TD), [], 2);
[maxValues2, columnIndices2] = max(abs(A_TD2), [], 2);
% columnIndices 就是每一行中最大值所在的列数,把列数置换成对应的多普勒频率
maxDF1 = (columnIndices - 101) * step_dop;
maxDF2 = (columnIndices2 - 101) * step_dop;
%% 初始位置确定
%确定物体和收发机的初始位置
% 定义收发机位置
%发射机位置 
xT = 0;
yT = 0;
%sur1的位置
xR1 = 3;
yR1 = 0;
%sur2的位置--假设sur2和tx呈45度
xR2 = sqrt(2);
yR2 = -sqrt(2);
%物体的初始位置
xTar = 2;
yTar = -sqrt(2);

%%角度初始化
fai_sur1 = zeros(1,length(array_start_time)-1);
fai_sur2 = zeros(1,length(array_start_time)-1);
fai_tx = zeros(1,length(array_start_time)-1);

fai_sur1(1) = atan((yTar)/ (xTar - xR1));
fai_sur2(1) = atan((yTar - yR2)/ (xTar - xR2));
fai_tx(1) = atan(yTar/xTar);

%存放每一个时刻目标的位置
xtar = zeros(1,length(array_start_time)-1);
ytar = zeros(1,length(array_start_time)-1);
xtar(1) = xTar;
ytar(1) = yTar;
%存放方程每一次迭代解的xy方向速度值
v_xy = zeros(length(array_start_time)-2,2);
%存放两条链路多普勒频率的临时变量
fd = zeros(1,2);
%解方程需要的F矩阵,参见ICCC论文
F = zeros(2,2);
%载波频率
fc = 60.48e9;
%光速
c = 3e8;
%开始解方程与迭代
for i = 1:1:length(array_start_time)-2
    fd = [maxDF1(i);maxDF2(i)];
    F = -2*fc/c*[cos((fai_sur1(i) - fai_tx(i))/2)*cos((fai_sur1(i) + fai_tx(i))/2) , ...
        cos((fai_sur1(i) - fai_tx(i))/2)*sin((fai_sur1(i) + fai_tx(i))/2); ...
        cos((fai_sur2(i) - fai_tx(i))/2)*cos((fai_sur2(i) + fai_tx(i))/2), ...
        cos((fai_sur2(i) - fai_tx(i))/2)*sin((fai_sur2(i) + fai_tx(i))/2)];
    v_xy(i,:) = F \ fd;
    %判断解是否正确，如果解出来离谱的速度值，要把它限制在一定的范围内
    % if v_xy(i,1) > 1 || v_xy(i,2) > 1
    %     v_xy(i,1) = 0.5;
    %     v_xy(i,2) = 0.5;
    end
    xtar(i+1) = xtar(i) + v_xy(i,1) * T_slide;
    ytar(i+1) = ytar(i) + v_xy(i,2) * T_slide;
    %更新下一时刻的AOA角度以便下一次进行迭代
    fai_sur1(i+1) = atan((ytar(i+1))/ (xtar(i+1) - xR1));
    fai_sur2(i+1) = atan((ytar(i+1) - yR2)/ (xtar(i+1) - xR2));
    fai_tx(i+1) = atan(ytar(i+1)/xtar(i+1));
end

%% 绘制运动物体的轨迹
% 绘制轨迹
fig3 = figure(3);
plot(xtar, ytar, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
% 设置图形标题和坐标轴标签
title('物体运动轨迹');
xlabel('X坐标');
ylabel('Y坐标');
% xlim([]);
% ylim([]);
% 显示网格
grid on;


