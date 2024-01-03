%% 
clear;
close all;
clc;
%addpath('fun');

chan_num_ref = 1;
chan_num_tar = 2;
total_duration= 4;

f_c = 0.5e9;
f_s = 1e6;

CIT = 0.1; %
CIT_region = CIT*f_s;

N_slide = 10;
T_slide = CIT / N_slide;% 0.1s N_slide 

max_dop = 1000;
step_dop = 1/CIT;
array_Doppler_frequency = -max_dop:step_dop:max_dop;


chan_num_total = chan_num_tar+chan_num_ref;

filename_data = "D:\Github\Passive-Handwriting-Tracking\data\0413\8_small.bin";

fid1=fopen(filename_data,'rb');  %rb - (1 Btye = 8 bits)
fseek(fid1,0,'eof');  %eof  
fsize = ftell(fid1);  %

total_samplelen = fsize/4;  % /8 * 2 IQ

sample_1channel = 1996; %unit   unitIQ unitchannel byte

samplelen = sample_1channel * (chan_num_total); %block = channel*
%sampleBytes = samplelen*4;
loopNum = floor(total_samplelen/samplelen); % = (/block)

samples_per_chan = sample_1channel * loopNum ;%channel:*()
%DurationTruth = samples_per_chan/f_s;
fclose(fid1);


%% Raw Rata
fid_cali=fopen(filename_data,'rb');%byte
fseek(fid_cali,0,'bof');%bof 
chan_seq_cali = zeros(chan_num_total,samples_per_chan);% 3xX 
for loopIdx = 1:loopNum
    [data,~]=fread(fid_cali,[2,sample_1channel*chan_num_total],'int16','l');% 2xXIQ?  int16 16
    for idx_chan = 0:chan_num_total-1
        frame_start = sample_1channel*(0 + idx_chan)+1;
        frame_end = sample_1channel*(1 + idx_chan);
        chan_seq_cali(idx_chan+1,(loopIdx-1)*sample_1channel+1:loopIdx*sample_1channel) = ...%
            data(1,frame_start:frame_end) + 1i*data(2,frame_start:frame_end);%IQ
    end
end
fclose(fid_cali);

%% :
% S_refS_tar1S_tar2S_ref
data_count = 2184*1;
S_ref = chan_seq_cali(3,1:data_count);
S_tar = chan_seq_cali(1,1:data_count);
h_temp = xcorr(S_ref,S_tar);
[~,col_h] = size(h_temp);
col_h = col_h +1;

% raw data
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
%% :+CAF 

%
array_start_time = 0:T_slide:total_duration-CIT;
A_TD = zeros(length(array_start_time),length(array_Doppler_frequency));
A_TD2 = zeros(length(array_start_time),length(array_Doppler_frequency));


%
ref_integer=zeros(length(array_start_time),CIT_region);
tar_integer=zeros(length(array_start_time),CIT_region);%fft

ref_integer2=zeros(length(array_start_time),CIT_region);
tar_integer2=zeros(length(array_start_time),CIT_region);%fft

for i= 1:length(array_start_time)-2
    ref_integer(i,:)=S_ref((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer(i,:)=S_tar1((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end                                                    %fft

for i= 1:length(array_start_time)-2
    ref_integer2(i,:)=S_ref((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
    tar_integer2(i,:)=S_tar2((round(array_start_time(i)*f_s+1)):(round(array_start_time(i)*f_s)+round(CIT*f_s)));
end  

%%%CAF
for i= 1:length(array_start_time)-2
    %%%
    tar_integer(i,:) = ClutterCancellation_Doppler(tar_integer(i,:),ref_integer(i,:));
    %%%CAF
    final=fftshift(fft(tar_integer(i,:).*conj(ref_integer(i,:)),CIT_region));
    A_TD(i,:) = final(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
    %%%
    tar_integer2(i,:) = ClutterCancellation_Doppler(tar_integer2(i,:),ref_integer2(i,:));
    %%%CAF
    final2=fftshift(fft(tar_integer2(i,:).*conj(ref_integer2(i,:)),CIT_region));
    A_TD2(i,:) = final2(CIT_region/2+1-max_dop/step_dop:CIT_region/2+1+max_dop/step_dop);
end


%% (CAF+CFAR)
% 
% A_TD2(95,:)=0;
% A_TD(95,:)=0;

%%% 1CAF
thres_A_TRD = -30;
fig1 = figure(1);
set(fig1,'position',[50,50,900,600]);
plot_A_DT = abs(A_TD');
plot_A_DT = mag2db(plot_A_DT/max(max(plot_A_DT)));
h1 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT);
xlim([array_start_time(1),array_start_time(end)]);
ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
% y
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



%%% CFAR
%0
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
%CFAR
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
%%% 2CAF
thres_A_TRD = -30;
fig3 = figure(3);
set(fig3,'position',[50,50,900,600]);
plot_A_DT2 = abs(A_TD2');
plot_A_DT2 = mag2db(plot_A_DT2/max(max(plot_A_DT2)));
h2 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT2);
xlim([array_start_time(1),array_start_time(end)]);
ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
% y
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
%%% CFAR
%0
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
%CFAR
detections_2 = cfar2D(resp,CUTIdx);
helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections_2)
% saveas(gcf, 'E:\0617\5s_'+string(p)+'-4cfar.jpg', 'jpg')
% saveas(gcf, 'C:\Users\Wu\Desktop\labsource\code\5f.jpg', 'jpg')
% close all;


%% CFARCAF
 Map_1 = zeros(size(resp));
 Map_2 = zeros(size(resp));
 Map_1(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2)) = ...
    reshape(double(detections_1),rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);
  Map_2(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2)) = ...
    reshape(double(detections_2),rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);
  %CFARCAF
foundColumns_1 = [];
foundColumns_2 = [];
  %0
for col = 1:size(Map_1, 2)
    % 0
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

%
[aligned_DF1, aligned_DF2] = align_matrices(resultDFList_1, resultDFList_2);

%% 
aligned_DF1 = KalmanSmoother(aligned_DF1);
aligned_DF2 = KalmanSmoother(aligned_DF2);

%% 
fig5 = figure(5);
plot(aligned_DF1, '-','Color', [1, 0.5, 0],'LineWidth', 3);
hold on;
plot(aligned_DF2, '-','Color', [0, 0.5, 1],'LineWidth', 3);

% 
% %%0
% nonZeroRowIndices_1 = cell(numel(foundColumns_1),2);
% nonZeroRowIndices_2 = cell(numel(foundColumns_1),2);
% % 
% for i = 1:numel(foundColumns_1)
%     col_1 = foundColumns_1(i);
%     col_2 = foundColumns_2(i);
%     % 0
%     nonZeroRows_1 = find(Map_1(:, col_1) ~= 0);
%     nonZeroRows_2 = find(Map_2(:, col_2) ~= 0);
%     % 
%     nonZeroRowIndices_1{i, 1} = col_1;
%     nonZeroRowIndices_1{i, 2} = nonZeroRows_1;
%     nonZeroRowIndices_2{i, 1} = col_2;
%     nonZeroRowIndices_2{i, 2} = nonZeroRows_2;
% end
% % 0CAF
% % 0
% plot_A_DT(101, :) = -1000;
% plot_A_DT2(101, :) = -1000;
% % 0
% maxRowIndices_1 = zeros(numel(foundColumns_1),2);
% maxRowIndices_2 = zeros(numel(foundColumns_2),2);
% for i = 1:size(nonZeroRowIndices_1, 1)
%     % 00
%     colIndex_1 = nonZeroRowIndices_1{i, 1};
%     rowIndex_1 = nonZeroRowIndices_1{i, 2};
%     colIndex_2 = nonZeroRowIndices_2{i, 1};
%     rowIndex_2 = nonZeroRowIndices_2{i, 2};
%     %  data1 
%     values_1 = plot_A_DT(rowIndex_1, colIndex_1);
%     values_2 = plot_A_DT2(rowIndex_2, colIndex_2);
%     % 
%     [~, maxIndex_1] = max(values_1);
%     [~, maxIndex_2] = max(values_2);
%     % 
%     maxRowIndices_1(i,1) = colIndex_1;
%     maxRowIndices_1(i,2) = rowIndex_1(maxIndex_1);
%     maxRowIndices_2(i,1) = colIndex_2;
%     maxRowIndices_2(i,2) = rowIndex_2(maxIndex_2);
% end
% % -Inf
% plot_A_DT = plot_A_DT(:, 1:end-2);
% plot_A_DT2 = plot_A_DT2(:,1:end-2);
% % 
% % A_TD(:, 101) = 0;
% % A_TD2(:, 101) = 0;
% %%delete_index
% delete_index = 40;
% plot_A_DT = plot_A_DT(:,delete_index:end);
% plot_A_DT2 = plot_A_DT2(:,delete_index:end);

%% 
%%
maxDF1 = aligned_DF1;
maxDF2 = aligned_DF2;
%
% 
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
xTar = 2.5;
yTar = -1.5;

sensing_time = array_start_time(1:numel(aligned_DF1));
%%
fai_sur1 = zeros(1,length(sensing_time)-1);
fai_sur2 = zeros(1,length(sensing_time)-1);
fai_tx = zeros(1,length(sensing_time)-1);

fai_sur1(1) = atan((yTar)/ (xTar - xR1));
fai_sur2(1) = atan((yTar - yR2)/ (xTar - xR2));
fai_tx(1) = atan(yTar/xTar);

%
xtar = zeros(1,length(sensing_time)-1);
ytar = zeros(1,length(sensing_time)-1);
xtar(1) = xTar;
ytar(1) = yTar;
%xy
v_xy = zeros(length(sensing_time)-2,2);
%
fd = zeros(1,2);
%F
F = zeros(2,2);
%
fc = 60.48e9;
%
c = 3e8;
%
for i = 1:1:length(sensing_time)-2
    fd = [maxDF1(i);maxDF2(i)];
    F = -2*fc/c*[cos((fai_sur1(i) - fai_tx(i))/2)*cos((fai_sur1(i) + fai_tx(i))/2) , ...
        cos((fai_sur1(i) - fai_tx(i))/2)*sin((fai_sur1(i) + fai_tx(i))/2); ...
        cos((fai_sur2(i) - fai_tx(i))/2)*cos((fai_sur2(i) + fai_tx(i))/2), ...
        cos((fai_sur2(i) - fai_tx(i))/2)*sin((fai_sur2(i) + fai_tx(i))/2)];
    v_xy(i,:) = F \ fd;
    %
    % if abs(v_xy(i,1)) > 0.6 
    %     v_xy(i,1) = (v_xy(i-1,1) + v_xy(i-2,1) + v_xy(i-3,1))/3;
    % elseif abs(v_xy(i,2)) > 0.6
    %     v_xy(i,2) = (v_xy(i-1,2) + v_xy(i-2,2) + v_xy(i-3,2))/3;
    % end
    xtar(i+1) = xtar(i) + v_xy(i,1) * T_slide;
    ytar(i+1) = ytar(i) + v_xy(i,2) * T_slide;
    %AOA
    fai_sur1(i+1) = atan((ytar(i+1))/ (xtar(i+1) - xR1));
    fai_sur2(i+1) = atan((ytar(i+1) - yR2)/ (xtar(i+1) - xR2));
    fai_tx(i+1) = atan(ytar(i+1)/xtar(i+1));
end

%% 
fig6 = figure(6);
plot(xtar, ytar, '-','Color', [1, 0.5, 0],'LineWidth', 3,'DisplayName', 'Estimated Trajectory');
% 
hold on;
plot(xtar(1), ytar(1), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'Initial Point');
plot(xtar(end), ytar(end), 'x', 'MarkerSize', 8, 'MarkerEdgeColor', 'b',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'End Point');
title('');
xlabel('X(m)');
ylabel('Y(m)');
hold off;
legend('show','Location', 'southwest');
% xlim([1.5 2]);
% ylim([-2 -1]);
% 
grid on;

%% Rotate and Mirror
% Created by Eason Hua
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
title('');
xlabel('X(m)');
ylabel('Y(m)');
% axis equal;
grid on;


%% Get Ground Truth
% Created by Eason Hua
close all;


% Read picture and convert it to binary image
img = imread('D:\Github\Passive-Handwriting-Tracking\imgs\eight\eight_hand.png');  
% RGB
Ih = rgb2gray(img); 
% 
thresh2 = graythresh(Ih); 
imgGray = im2bw(Ih,thresh2);


% Cut the figure to drop the useless white
% 
[row, col] = find(imgGray == 0);
left = min(col);
right = max(col);
top = min(row);
bottom = max(row);
% 
croppedImg = imcrop(imgGray, [left, top, right - left, bottom - top]);


%% Adjust the data to expected range
% Created by Eason Hua
[Y, X] = find(croppedImg==0);
Y = -Y;

lowerLeft = [min(rotated(1, :)), min(rotated(2, :))];
upperRight = [max(rotated(1, :)), max(rotated(2, :))];
% 
targetX = [lowerLeft(1), upperRight(1)];
%  [0, 1] 
normalizedX = (X - min(X)) / (max(X) - min(X));
% 
mappedX = normalizedX * (targetX(2) - targetX(1)) + targetX(1);

% 
targetY = [lowerLeft(2), upperRight(2)];
%  [0, 1] 
normalizedY = (Y - min(Y)) / (max(Y) - min(Y));
% 
mappedY = normalizedY * (targetY(2) - targetY(1)) + targetY(1);


% Under sample to match the length
% 计算欠采样因子
undersampleFactor = floor(length(X) / length(rotated(1,:)));
% 对序列进行欠采样
undersampledX = X(1:undersampleFactor:end);
undersampledY = Y(1:undersampleFactor:end);


%% Trajectory to scatter
% Created by Eason Hua
matrix = ones(size(croppedImg,1), size(croppedImg,2));
for i = 1:1:length(rotated(1,:))
    point = [rotated(1,i), rotated(2,i)];
    pointMatrix = [0,0];

    pointMatrix(1) = (point(1)-lowerLeft(1))/(upperRight(1)-lowerLeft(1))*size(croppedImg,2);
    pointMatrix(1) = ceil(pointMatrix(1)+1);
    if pointMatrix(1) > size(croppedImg,2)
        pointMatrix(1) = size(croppedImg,2);
    end

    pointMatrix(2) = (upperRight(2)-point(2))/(upperRight(2)-lowerLeft(2))*size(croppedImg,1);
    pointMatrix(2) = ceil(pointMatrix(2)+1);
    if pointMatrix(2) > size(croppedImg,1)
        pointMatrix(2) = size(croppedImg,1);
    end

    matrix(pointMatrix(2), pointMatrix(1)) = 0;
end

[matrixY, matrixX] = find(matrix==0);
matrixY = -matrixY;


%% Draw trajectories
figure;
% scatter(mappedX, mappedY, ...
%     'Marker', '.', ...
%     'MarkerEdgeColor', [0.5, 1, 1], ...
%     'DisplayName', 'Handwritten Trajectory');
% hold on;

% plot(undersampledX, undersampledY, ...
%     '-', ...
%     'Color', [1, 0.5, 0], ...
%     'LineWidth', 0.01, ...
%     'DisplayName', 'Estimated Trajectory');
% hold on;

scatter(undersampledX, undersampledY, ...
    'Marker', '.', ...
    'MarkerEdgeColor', [0.5, 1, 1], ...
    'DisplayName', 'Handwritten Trajectory');
hold on;

plot(rotated(1, :), rotated(2, :), ...
    '-', ...
    'Color', [1, 0.5, 0], ...
    'LineWidth', 3, ...
    'DisplayName', 'Estimated Trajectory');
hold on;

plot(rotated(1, 1), rotated(2, 1), ...
    'x', ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'r',...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2.5, ...
    'DisplayName', 'Initial Point');
hold on;

plot(rotated(1, end), rotated(2, end), ...
    'x', ...
    'MarkerSize', 8, ...
    'MarkerEdgeColor', 'b',...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 2.5, ...
    'DisplayName', 'End Point');
hold off;

legend('show');
grid on;
% 
% axis equal; 


%% Draw scatters
figure;
scatter(matrixX,size(croppedImg,1)+matrixY, ...
    'Marker', '.', ...
    'MarkerEdgeColor', [1, 0.5, 0], ...
    'DisplayName', 'Estimated Trajectory');
hold on;

scatter(undersampledX, size(croppedImg,1)+undersampledY, ...
    'Marker', '.', ...
    'MarkerEdgeColor', [0, 0, 1], ...
    'DisplayName', 'Handwritten Trajectory');
hold off;

xlim([0, size(croppedImg,2)]);
ylim([0, size(croppedImg,1)]);
legend('show');
grid on;


%% Evaluate similarity

% 计算散点图的累积分布函数
cdf_estimated = cumsum(histcounts2(matrixX, matrixY, 'Normalization', 'cdf'));
cdf_handwritten = cumsum(histcounts2(undersampledX, undersampledY, 'Normalization', 'cdf'));

% 绘制CDF曲线
figure;
plot(linspace(0, 1, length(cdf_estimated)), cdf_estimated, 'LineWidth', 2, 'DisplayName', 'Estimated Trajectory');
hold on;
plot(linspace(0, 1, length(cdf_handwritten)), cdf_handwritten, 'LineWidth', 2, 'DisplayName', 'Handwritten Trajectory');
hold off;

% 设置图例等其他图形属性
legend('Estimated Trajectory', 'Handwritten Trajectory');
xlabel('Similarity');
ylabel('Cumulative Probability');
title('Cumulative Distribution Function (CDF) of Trajectories');

