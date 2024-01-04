%% 初始参数配置
clear;
close all;
clc;
%addpath('fun');

chan_num_ref = 1;
chan_num_tar = 2;
total_duration= 4;

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

filename_data = "E:\Desktop\Project\Data\Right_data\3.bin";

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


%% 读取文件，获得三路信道的Raw Rata
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

%% 数据预处理:检查两行数据是否对齐，如果没对齐则裁掉错开的地方
% 假设S_ref是参考信号，对齐S_tar1和S_tar2到S_ref的时间轴
data_count = 2184*1;
S_ref = chan_seq_cali(3,1:data_count);
S_tar = chan_seq_cali(1,1:data_count);
h_temp = xcorr(S_ref,S_tar);
[~,col_h] = size(h_temp);
col_h = col_h +1;

% 提取三路信道的raw data
N = floor(CIT*f_s)/1000-1;
col_max = find(max(abs(h_temp(col_h/2-N:col_h/2+N))) == abs(h_temp(col_h/2-N:col_h/2+N)));
array_sample_shift = col_max -N -1;
if array_sample_shift>0
    S_ref = chan_seq_cali(1,1+array_sample_shift:end);
    S_tar1 = chan_seq_cali(2,1:end-array_sample_shift);
    S_tar2 = chan_seq_cali(3,1:end-array_sample_shift);
else
    S_ref = chan_seq_cali(1,1:end+array_sample_shift);
    S_tar1 = chan_seq_cali(2,1-array_sample_shift:end);
    S_tar2 = chan_seq_cali(3,1-array_sample_shift:end);
end
%% 数据处理:杂波消除+CAF 

%滑动窗口
array_start_time = 0:T_slide:total_duration-CIT;
A_TD = zeros(length(array_start_time),length(array_Doppler_frequency));
A_TD2 = zeros(length(array_start_time),length(array_Doppler_frequency));


%数据处理
ref_integer=zeros(length(array_start_time),CIT_region);
tar_integer=zeros(length(array_start_time),CIT_region);%准备fft数组表

ref_integer2=zeros(length(array_start_time),CIT_region);
tar_integer2=zeros(length(array_start_time),CIT_region);%准备fft数组表

for i= 1:length(array_start_time)-2
    ref_integer(i,:)=S_ref((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer(i,:)=S_tar1((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end                                                    %分类并填写fft数组表

for i= 1:length(array_start_time)-2
    ref_integer2(i,:)=S_ref((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer2(i,:)=S_tar2((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end  

%%%滑动窗口式执行杂波消除和CAF
for i= 1:length(array_start_time)-2
    %%%杂波消除
    tar_integer(i,:) = ClutterCancellation_Doppler(tar_integer(i,:),ref_integer(i,:));
    %%%CAF
    final=fftshift(fft(tar_integer(i,:).*conj(ref_integer(i,:)),CIT_region));
    A_TD(i,:) = final(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
    %%%杂波消除
    tar_integer2(i,:) = ClutterCancellation_Doppler(tar_integer2(i,:),ref_integer2(i,:));
    %%%CAF
    final2=fftshift(fft(tar_integer2(i,:).*conj(ref_integer2(i,:)),CIT_region));
    A_TD2(i,:) = final2(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
end


%% 数据绘图(CAF+CFAR)
% 
% A_TD2(95,:)=0;
% A_TD(95,:)=0;

%%% 画参考与监视1的CAF
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
title('1 and 2')
xlabel('Time (s)')
ylabel('Doppler frequency (Hz)')
colormap('jet');
caxis([thres_A_TRD,0]);
% saveas(gcf, 'E:\0617\5s_'+string(p)+'-1fft.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5a.jpg', 'jpg')



%%% CFAR算法
%提前消除0频最大多普勒分量
% plot_A_DT(101, :) = -1000;
cfar2D = phased.CFARDetector2D('GuardBandSize',5,'TrainingBandSize',5,...
  'ProbabilityFalseAlarm',0.7);
resp=plot_A_DT;
rngGrid=array_Doppler_frequency.';
dopGrid=array_start_time.';
rangeIndx(1)=60;
rangeIndx(2)=140;
dopplerIndx(1)= 11;
dopplerIndx(2)= 379;
[columnInds,rowInds] = meshgrid(dopplerIndx(1):dopplerIndx(2),...
  rangeIndx(1):rangeIndx(2));
CUTIdx = [rowInds(:) columnInds(:)]';
%CFAR检测结果
detections_1 = cfar2D(resp,CUTIdx);
helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections_1)
% %saveas(gcf, 'E:\0617\5s_'+string(p)+'-2cfar.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5c.jpg', 'jpg')
% 
% 
% 
% 
% 
% 
% 
%%% 画参考与监视2的CAF
thres_A_TRD = -30;
fig3 = figure(3);
set(fig3,'position',[50,50,900,600]);
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
title('1 and 3')
colormap('jet');
caxis([thres_A_TRD,0]);
%saveas(gcf, 'E:\0617\5s_'+string(p)+'-3fft.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5b.jpg', 'jpg')
% 
% 
% 
% 
%%% CFAR算法
%提前消除0频最大多普勒分量
% plot_A_DT2(101, :) = -1000;
cfar2D = phased.CFARDetector2D('GuardBandSize',5,'TrainingBandSize',5,...
  'ProbabilityFalseAlarm',0.7);
resp=plot_A_DT2;
rngGrid=array_Doppler_frequency.';
dopGrid=array_start_time.';
rangeIndx(1)= 60;
rangeIndx(2)= 140;
dopplerIndx(1)= 11;
dopplerIndx(2)= 379;
[columnInds,rowInds] = meshgrid(dopplerIndx(1):dopplerIndx(2),...
  rangeIndx(1):rangeIndx(2));
CUTIdx = [rowInds(:) columnInds(:)]';
%CFAR检测结果
detections_2 = cfar2D(resp,CUTIdx);
helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections_2)
% saveas(gcf, 'E:\0617\5s_'+string(p)+'-4cfar.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5f.jpg', 'jpg')
% close all;


%% 提取CFAR结果并将其映射到CAF结果中的每一列最大值的多普勒频率
 Map_1 = zeros(size(resp));
 Map_2 = zeros(size(resp));
 Map_1(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2)) = ...
    reshape(double(detections_1),rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);
  Map_2(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2)) = ...
    reshape(double(detections_2),rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);
  %将CFAR得到的二进制数组矩阵映射到之前CAF结果对应的行和列中
foundColumns_1 = [];
foundColumns_2 = [];
  %遍历每一列，找到元素不全为0的列，证明从该列开始检测到多普勒频率了，要记录该列的索引
for col = 1:size(Map_1, 2)
    % 检查当前列是否所有行的元素都不为0
    if any(Map_1(:, col) ~= 0) && any(Map_2(:, col) ~= 0)
        foundColumns_1 = [foundColumns_1, col];
        foundColumns_2 = [foundColumns_2, col];
    end
end

%%%The goal of the path matching algorithm is to increase the number of Doppler points calculated by CFAR to 369
%%path matching algorithm
endCol_first_1 = 0;
endCol_first_2 = 0;
resultDFList_1 = zeros(1,dopplerIndx(2)-dopplerIndx(1)+1);
resultDFList_2 = zeros(1,dopplerIndx(2)-dopplerIndx(1)+1);
% Extract and process 11st column
lastCol = Map_1(:,dopplerIndx(1));
lastCol_2 = Map_2(:,dopplerIndx(1));
nonZeroRow = find(lastCol ~= 0);
nonZeroRow_2 = find(lastCol_2 ~= 0);
%Weight parameter when there are multiple Doppler frequencies
alpha_1 = 10;
beta_1 = 1;
alpha_2 = 10;
beta_2 = 1;
%%Find the first interpolation Doppler of the two channels
if isempty(nonZeroRow) % If the first moment is 0, you have to look later
    for i = dopplerIndx(1)+1 : dopplerIndx(2)
            testCol_first = Map_1(:,i);
            testNonZeroRow_first = find(testCol_first ~= 0);
            if isempty(testNonZeroRow_first)
                continue
            else
                endCol_first_1 = i;
                endRow_first_1 = mean(Row2Df(testNonZeroRow_first));
                break
            end
    end
    resultDF = endRow_first_1; 
    resultDFList_1(endCol_first_1-dopplerIndx(1)+1) = resultDF;
    resultDFList_1(1:endCol_first_1-dopplerIndx(1)) = [];
else
    endCol_first_1 = dopplerIndx(1);
    resultDF = mean(Row2Df(nonZeroRow));
    resultDFList_1(1) = resultDF;
end

if isempty(nonZeroRow_2) % If the first moment is 0, you have to look later
    for i = dopplerIndx(1)+1 : dopplerIndx(2)
            testCol_first = Map_2(:,i);
            testNonZeroRow_first = find(testCol_first ~= 0);
            if isempty(testNonZeroRow_first)
                continue
            else
                endCol_first_2 = i;
                endRow_first_2 = mean(Row2Df(testNonZeroRow_first));
                break
            end
    end
    resultDF = endRow_first_2; %%Insert 0 Doppler
    resultDFList_2(endCol_first_2-dopplerIndx(1)+1) = resultDF;
    resultDFList_2(1:endCol_first_2-dopplerIndx(1)) = [];
else
    endCol_first_2 = dopplerIndx(1);
    resultDF = mean(Row2Df(nonZeroRow_2));
    resultDFList_2(1) = resultDF;
end

% Extract and process endCol_first+1 st to 379th columns for channel 1
endCol_1 = 0;
for i = endCol_first_1+1 : dopplerIndx(2) %12 to 379
    if i < endCol_1+1
        continue
    end

    thisCol = Map_1(:,i);
    nonZeroRow = find(thisCol ~= 0);
    if isempty(nonZeroRow)
        startCol = i;
        % Find forwad until there exists a unempty column
        for j = i+1 : dopplerIndx(2)
            testCol = Map_1(:,j);
            testNonZeroRow = find(testCol ~= 0);
            if isempty(testNonZeroRow)
                continue
            else
                endCol_1 = j-1;
                endRow = mean(Row2Df(testNonZeroRow));
                break
            end
        end
        % Assign a value for sequence, 
        % with index from startCol to endCol, 
        % and value from resultRow to endRow
        slope = (endRow-resultDF)/((endCol_1-startCol+2)*T_slide);
        for k = startCol:endCol_1
            resultDFList_1(k-endCol_first_1+1) =  resultDF + slope*((k-startCol+1)*T_slide);
        end

    elseif length(nonZeroRow) == 1  %Just one
        % Simply take it
        resultDF = Row2Df(nonZeroRow);
        resultDFList_1(i-endCol_first_1+1) = resultDF;
    else  %There's multiple Dopplers
        % Find the nearest index
        distance = abs(Row2Df(nonZeroRow) - resultDF); % distance
        tempolate = alpha_1 * exp(-distance) + beta_1 * abs(Row2Df(nonZeroRow));
        [miniRow, miniCol] = find(tempolate==max(tempolate));
        resultDF = Row2Df(nonZeroRow(miniRow(1)));
        resultDFList_1(i-endCol_first_1+1) = resultDF;
    end
end

% Extract and process endCol_first+1 st to 379th columns for channel 2
endCol_2 = 0;
for i = endCol_first_2+1 : dopplerIndx(2) %12 to 379
    if i < endCol_2+1
        continue
    end

    thisCol = Map_2(:,i);
    nonZeroRow = find(thisCol ~= 0);
    if isempty(nonZeroRow)
        startCol = i;
        % Find forwad until there exists a unempty column
        for j = i+1 : dopplerIndx(2)
            testCol = Map_2(:,j);
            testNonZeroRow = find(testCol ~= 0);
            if isempty(testNonZeroRow)
                continue
            else
                endCol_2 = j-1;
                endRow = mean(Row2Df(testNonZeroRow));
                break
            end
        end
        % Assign a value for sequence, 
        % with index from startCol to endCol, 
        % and value from resultRow to endRow
        slope = (endRow-resultDF)/((endCol_2-startCol+2)*T_slide);
        for k = startCol:endCol_2
            resultDFList_2(k-endCol_first_2+1) =  resultDF + slope*((k-startCol+1)*T_slide);
        end

    elseif length(nonZeroRow) == 1  %just one
        % Simply take it
        resultDF = Row2Df(nonZeroRow);
        resultDFList_2(i-endCol_first_2+1) = resultDF;
    else  %There's multiple Dopplers
        % Find the nearest index
        distance = abs(Row2Df(nonZeroRow) - resultDF); % distance
        tempolate = alpha_2 * exp(-distance) + beta_2 * abs(Row2Df(nonZeroRow));
        [miniRow, miniCol] = find(tempolate==max(tempolate));
        resultDF = Row2Df(nonZeroRow(miniRow(1)));
        resultDFList_2(i-endCol_first_2+1) = resultDF;
    end
end

%如果两路信道的多普勒矩阵没有对齐，就要让他们对齐
[aligned_DF1, aligned_DF2] = align_matrices(resultDFList_1, resultDFList_2);

%%加卡尔曼滤波平滑
aligned_DF1 = KalmanSmoother(aligned_DF1);
aligned_DF2 = KalmanSmoother(aligned_DF2);

fig5 = figure(5);
plot((0.01:0.01:round(numel(aligned_DF1)/100,2))-0.01,aligned_DF1, '-','Color', [1, 0.5, 0],'LineWidth', 3, 'DisplayName', 'Sur1');
hold on;
plot((0.01:0.01:round(numel(aligned_DF1)/100,2))-0.01,aligned_DF2, '-','Color', [0, 0.5, 1],'LineWidth', 3, 'DisplayName', 'Sur2');
title('路径匹配算法与卡尔曼滤波后的时间-多普勒图');
xlabel('时间(s)');
ylabel('多普勒频移(Hz)');
legend('show','Location', 'southwest');
% 
% %%记录这些不全为0的列的行索引
% nonZeroRowIndices_1 = cell(numel(foundColumns_1),2);
% nonZeroRowIndices_2 = cell(numel(foundColumns_1),2);
% % 遍历每一个找到的列
% for i = 1:numel(foundColumns_1)
%     col_1 = foundColumns_1(i);
%     col_2 = foundColumns_2(i);
%     % 找到当前列中不为0的元素的行索引
%     nonZeroRows_1 = find(Map_1(:, col_1) ~= 0);
%     nonZeroRows_2 = find(Map_2(:, col_2) ~= 0);
%     % 存储到单元格数组中
%     nonZeroRowIndices_1{i, 1} = col_1;
%     nonZeroRowIndices_1{i, 2} = nonZeroRows_1;
%     nonZeroRowIndices_2{i, 1} = col_2;
%     nonZeroRowIndices_2{i, 2} = nonZeroRows_2;
% end
% % 将这些不为0的列对应的行索引映射到之前画CAF的矩阵中，并比较他们的最大值，取最大值的行索引对应当前列的多普勒频移
% % 处理前先抹掉0频处的值
% plot_A_DT(101, :) = -1000;
% plot_A_DT2(101, :) = -1000;
% % 遍历每一个不为0的列
% maxRowIndices_1 = zeros(numel(foundColumns_1),2);
% maxRowIndices_2 = zeros(numel(foundColumns_2),2);
% for i = 1:size(nonZeroRowIndices_1, 1)
%     % 获取当前不为0的列的列索引和对应列中不为0的行的行索引
%     colIndex_1 = nonZeroRowIndices_1{i, 1};
%     rowIndex_1 = nonZeroRowIndices_1{i, 2};
%     colIndex_2 = nonZeroRowIndices_2{i, 1};
%     rowIndex_2 = nonZeroRowIndices_2{i, 2};
%     % 从 data1 中获取对应列中这些行的值
%     values_1 = plot_A_DT(rowIndex_1, colIndex_1);
%     values_2 = plot_A_DT2(rowIndex_2, colIndex_2);
%     % 找到最大值对应的行索引
%     [~, maxIndex_1] = max(values_1);
%     [~, maxIndex_2] = max(values_2);
%     % 存储最大值对应的行索引和对应列索引
%     maxRowIndices_1(i,1) = colIndex_1;
%     maxRowIndices_1(i,2) = rowIndex_1(maxIndex_1);
%     maxRowIndices_2(i,1) = colIndex_2;
%     maxRowIndices_2(i,2) = rowIndex_2(maxIndex_2);
% end
% % 最后两列数据是-Inf，删除
% plot_A_DT = plot_A_DT(:, 1:end-2);
% plot_A_DT2 = plot_A_DT2(:,1:end-2);
% % 删除零频附近的值
% % A_TD(:, 101) = 0;
% % A_TD2(:, 101) = 0;
% %%重要，这个删去无关信息的起始时间delete_index需要自适应地动态调整
% delete_index = 40;
% plot_A_DT = plot_A_DT(:,delete_index:end);
% plot_A_DT2 = plot_A_DT2(:,delete_index:end);

%% 初始位置确定
%%两个链路多普勒频率矩阵的确定
maxDF1 = aligned_DF1;
maxDF2 = aligned_DF2;
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
xTar = 2.5;
yTar = -1.5;

sensing_time = array_start_time(1:numel(aligned_DF1));
%%角度初始化
fai_sur1 = zeros(1,length(sensing_time)-1);
fai_sur2 = zeros(1,length(sensing_time)-1);
fai_tx = zeros(1,length(sensing_time)-1);

fai_sur1(1) = atan((yTar)/ (xTar - xR1));
fai_sur2(1) = atan((yTar - yR2)/ (xTar - xR2));
fai_tx(1) = atan(yTar/xTar);

%存放每一个时刻目标的位置
xtar = zeros(1,length(sensing_time)-1);
ytar = zeros(1,length(sensing_time)-1);
xtar(1) = xTar;
ytar(1) = yTar;
%存放方程每一次迭代解的xy方向速度值
v_xy = zeros(length(sensing_time)-2,2);
%存放两条链路多普勒频率的临时变量
fd = zeros(1,2);
%解方程需要的F矩阵
F = zeros(2,2);
%载波频率
fc = 60.48e9;
%光速
c = 3e8;
%开始解方程与迭代
for i = 1:1:length(sensing_time)-2
    fd = [maxDF1(i);maxDF2(i)];
    F = -2*fc/c*[cos((fai_sur1(i) - fai_tx(i))/2)*cos((fai_sur1(i) + fai_tx(i))/2) , ...
        cos((fai_sur1(i) - fai_tx(i))/2)*sin((fai_sur1(i) + fai_tx(i))/2); ...
        cos((fai_sur2(i) - fai_tx(i))/2)*cos((fai_sur2(i) + fai_tx(i))/2), ...
        cos((fai_sur2(i) - fai_tx(i))/2)*sin((fai_sur2(i) + fai_tx(i))/2)];
    v_xy(i,:) = F \ fd;
    %判断解是否正确，如果解出来离谱的速度值，就要用之前的值进行加权平均
    % if abs(v_xy(i,1)) > 0.6 
    %     v_xy(i,1) = (v_xy(i-1,1) + v_xy(i-2,1) + v_xy(i-3,1))/3;
    % elseif abs(v_xy(i,2)) > 0.6
    %     v_xy(i,2) = (v_xy(i-1,2) + v_xy(i-2,2) + v_xy(i-3,2))/3;
    % end
    xtar(i+1) = xtar(i) + v_xy(i,1) * T_slide;
    ytar(i+1) = ytar(i) + v_xy(i,2) * T_slide;
    %更新下一时刻的AOA角度以便下一次进行迭代
    fai_sur1(i+1) = atan((ytar(i+1))/ (xtar(i+1) - xR1));
    fai_sur2(i+1) = atan((ytar(i+1) - yR2)/ (xtar(i+1) - xR2));
    fai_tx(i+1) = atan(ytar(i+1)/xtar(i+1));
end

%% 绘制运动物体的轨迹
% 绘制轨迹
fig6 = figure(6);
plot(xtar, ytar, '-','Color', [1, 0.5, 0],'LineWidth', 3,'DisplayName', 'Estimated Trajectory');
% 设置图形标题和坐标轴标签
hold on;
plot(xtar(1), ytar(1), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'Initial Point');
plot(xtar(end), ytar(end), 'x', 'MarkerSize', 8, 'MarkerEdgeColor', 'b',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'End Point');
title('物体运动轨迹');
xlabel('X坐标(m)');
ylabel('Y坐标(m)');
hold off;
legend('show','Location', 'southwest');
% xlim([1.5 2]);
% ylim([-2 -1]);
% 显示网格
grid on;

%% Rotate and Mirror
coordinates = [xtar; ytar];
center = mean(coordinates, 2);
translated = coordinates - center;

theta_rotate = 7/12 * pi;
R = [cos(theta_rotate), -sin(theta_rotate); sin(theta_rotate), cos(theta_rotate)];
rotated = R * translated;

% Mirror horizontally
rotated = [-rotated(1,:) ; rotated(2,:)];
rotated = rotated + center;
fig7 = figure(7);

plot(rotated(1, :), rotated(2, :), '-','Color', [1, 0.5, 0],'LineWidth', 3,'DisplayName', 'Estimated Trajectory');
hold on;
plot(rotated(1, 1), rotated(2, 1), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'Initial Point');
plot(rotated(1, end), rotated(2, end), 'x', 'MarkerSize', 8, 'MarkerEdgeColor', 'b',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'End Point');
hold off;

legend('show');
title('物体运动轨迹');
xlabel('X坐标(m)');
ylabel('Y坐标(m)');
% axis equal;
grid on;

% %% save the position data
% temp_save(1,:) = xtar;
% temp_save(2,:) = ytar;
% save('./results_data/2.mat', 'temp_save');

% %% Get Ground Truth
% % Read picture and convert it to binary image
% img = imread("E:\Desktop\Project\imgs\eight\eight_hand.png");  
% 
% Ih = rgb2gray(img); 
% 
% thresh2 = graythresh(Ih); 
% imgGray = im2bw(Ih,thresh2);
% 
% 
% % Cut the figure to drop the useless white
% 
% [row, col] = find(imgGray == 0);
% left = min(col);
% right = max(col);
% top = min(row);
% bottom = max(row);
% 
% croppedImg = imcrop(imgGray, [left, top, right - left, bottom - top]);
% 
% 
% %% Adjust the data to expected range
% [Y, X] = find(croppedImg==0);
% Y = -Y;
% 
% lowerLeft = [min(rotated(1, :)), min(rotated(2, :))];
% upperRight = [max(rotated(1, :)), max(rotated(2, :))];
% 
% targetX = [lowerLeft(1), upperRight(1)];
% 
% normalizedX = (X - min(X)) / (max(X) - min(X));
% 
% mappedX = normalizedX * (targetX(2) - targetX(1)) + targetX(1);
% 
% 
% targetY = [lowerLeft(2), upperRight(2)];
% 
% normalizedY = (Y - min(Y)) / (max(Y) - min(Y));
% 
% mappedY = normalizedY * (targetY(2) - targetY(1)) + targetY(1);
% 
% 
% 
% %% Trajectory to scatter
% matrix = ones(size(croppedImg,1), size(croppedImg,2));
% for i = 1:1:length(rotated(1,:))
%     point = [rotated(1,i), rotated(2,i)];
%     pointMatrix = [0,0];
% 
%     pointMatrix(1) = (point(1)-lowerLeft(1))/(upperRight(1)-lowerLeft(1))*size(croppedImg,2);
%     pointMatrix(1) = ceil(pointMatrix(1)+1);
%     if pointMatrix(1) > size(croppedImg,2)
%         pointMatrix(1) = size(croppedImg,2);
%     end
% 
%     pointMatrix(2) = (upperRight(2)-point(2))/(upperRight(2)-lowerLeft(2))*size(croppedImg,1);
%     pointMatrix(2) = ceil(pointMatrix(2)+1);
%     if pointMatrix(2) > size(croppedImg,1)
%         pointMatrix(2) = size(croppedImg,1);
%     end
% 
%     matrix(pointMatrix(2), pointMatrix(1)) = 0;
% end
% 
% [matrixY, matrixX] = find(matrix==0);
% matrixY = -matrixY;
% 
% 
% % Under sample estimated to match the length
% 
% undersampleFactor = floor(length(X) / length(matrixX));
% 
% undersampledX = X(1:undersampleFactor:end);
% undersampledY = Y(1:undersampleFactor:end);
% undersampledX = undersampledX(1:length(matrixX),:);
% undersampledY = undersampledY(1:length(matrixY),:);
% 
% 
% %% Draw trajectories
% fig8 = figure(8);
% % scatter(mappedX, mappedY, ...
% %     'Marker', '.', ...
% %     'MarkerEdgeColor', [0.5, 1, 1], ...
% %     'DisplayName', 'Handwritten Trajectory');
% % hold on;
% 
% % plot(undersampledX, undersampledY, ...
% %     '-', ...
% %     'Color', [1, 0.5, 0], ...
% %     'LineWidth', 0.01, ...
% %     'DisplayName', 'Estimated Trajectory');
% % hold on;
% 
% scatter(undersampledX, undersampledY, ...
%     'Marker', '.', ...
%     'MarkerEdgeColor', [0.5, 1, 1], ...
%     'DisplayName', 'Handwritten Trajectory');
% hold on;
% 
% plot(rotated(1, :), rotated(2, :), ...
%     '-', ...
%     'Color', [1, 0.5, 0], ...
%     'LineWidth', 3, ...
%     'DisplayName', 'Estimated Trajectory');
% hold on;
% 
% plot(rotated(1, 1), rotated(2, 1), ...
%     'x', ...
%     'MarkerSize', 8, ...
%     'MarkerEdgeColor', 'r',...
%     'MarkerFaceColor', 'none', ...
%     'LineWidth', 2.5, ...
%     'DisplayName', 'Initial Point');
% hold on;
% 
% plot(rotated(1, end), rotated(2, end), ...
%     'x', ...
%     'MarkerSize', 8, ...
%     'MarkerEdgeColor', 'b',...
%     'MarkerFaceColor', 'none', ...
%     'LineWidth', 2.5, ...
%     'DisplayName', 'End Point');
% hold off;
% 
% legend('show');
% grid on;
% 
% % axis equal; 
% 
% 
% %% Draw scatters
% fig9 = figure(9);
% scatter(matrixX,size(croppedImg,1)+matrixY, ...
%     'Marker', '.', ...
%     'MarkerEdgeColor', [1, 0.5, 0], ...
%     'DisplayName', 'Estimated Trajectory');
% hold on;
% 
% scatter(undersampledX, size(croppedImg,1)+undersampledY, ...
%     'Marker', '.', ...
%     'MarkerEdgeColor', [0, 0, 1], ...
%     'DisplayName', 'Handwritten Trajectory');
% hold off;
% 
% xlim([0, size(croppedImg,2)]);
% ylim([0, size(croppedImg,1)]);
% legend('show');
% grid on;
% 
% 
% %% Evaluate similarity
% 
% 
% cdf_estimated = cumsum(histcounts2(matrixX, matrixY, 'Normalization', 'cdf'));
% cdf_handwritten = cumsum(histcounts2(undersampledX, undersampledY, 'Normalization', 'cdf'));
% 
% 
% fig10 = figure(10);
% plot(linspace(0, 1, length(cdf_estimated)), cdf_estimated, 'LineWidth', 2, 'DisplayName', 'Estimated Trajectory');
% hold on;
% plot(linspace(0, 1, length(cdf_handwritten)), cdf_handwritten, 'LineWidth', 2, 'DisplayName', 'Handwritten Trajectory');
% hold off;
% 
% 
% legend('Estimated Trajectory', 'Handwritten Trajectory');
% xlabel('Similarity');
% ylabel('Cumulative Probability');
% title('Cumulative Distribution Function (CDF) of Trajectories');

