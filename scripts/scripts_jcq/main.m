%% 
clear;
close all;
clc;
%addpath('fun');

%% 
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

filename_data = "D:\Github\Passive-Handwriting-Tracking\data\0413\seven_small.bin";

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


%% CAF
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
plot_A_DT = abs(A_TD');
plot_A_DT = mag2db(plot_A_DT/max(max(plot_A_DT)));

% fig1 = figure(1);
% set(fig1,'position',[50,50,900,600]);

% h1 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT);
% xlim([array_start_time(1),array_start_time(end)]);
% % ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
% ylim([-300,300]);
% % y
% set(gca, 'YDir', 'normal');
% set(gcf,'unit','centimeters','position',[5 3 30 15]);
% set(get(gca,'XLabel'),'FontSize',22);
% set(get(gca,'YLabel'),'FontSize',22);
% colorbar;
% title('1 and 2')
% xlabel('Time (s)')
% ylabel('Doppler frequency (Hz)')
% colormap('jet');
% clim([-30,0]);


%%% CFAR
% Modifided by Eason Hua
ProbabilityFalseAlarm = 0.6;
GuardBandSize = 5;
TrainingBandSize = 5;


% plot_A_DT(101, :) = -1000;
cfar2D = phased.CFARDetector2D('GuardBandSize',GuardBandSize, ...
                                'TrainingBandSize',TrainingBandSize, ...
                                'ProbabilityFalseAlarm',ProbabilityFalseAlarm);

resp=plot_A_DT;
rngGrid=array_Doppler_frequency.';
dopGrid=array_start_time.';
rangeIndx(1)=70;
rangeIndx(2)=130;
dopplerIndx(1)= 11;
dopplerIndx(2)= 379;
[columnInds,rowInds] = meshgrid(dopplerIndx(1):dopplerIndx(2),...
  rangeIndx(1):rangeIndx(2));
CUTIdx = [rowInds(:) columnInds(:)]';
%CFAR
detections_1 = cfar2D(resp,CUTIdx);
helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections_1)



%%% 2CAF
plot_A_DT2 = abs(A_TD2');
plot_A_DT2 = mag2db(plot_A_DT2/max(max(plot_A_DT2)));

% fig3 = figure(3);
% set(fig3,'position',[50,50,900,600]);

% h2 = imagesc(array_start_time,array_Doppler_frequency,plot_A_DT2);
% xlim([array_start_time(1),array_start_time(end)]);
% % ylim([array_Doppler_frequency(1),array_Doppler_frequency(end)]);
% ylim([-300,300]);
% % y
% set(gca, 'YDir', 'normal');
% set(gcf,'unit','centimeters','position',[5 3 30 15]);
% set(get(gca,'XLabel'),'FontSize',22);
% set(get(gca,'YLabel'),'FontSize',22);
% colorbar;
% xlabel('Time (s)')
% ylabel('Doppler frequency (Hz)')
% title('1 and 3')
% colormap('jet');
% clim([-30,0]);



%%% CFAR
% plot_A_DT2(101, :) = -1000;
% Modifided by Eason Hua
cfar2D = phased.CFARDetector2D('GuardBandSize',GuardBandSize, ...
                                'TrainingBandSize',TrainingBandSize,...
                              'ProbabilityFalseAlarm',ProbabilityFalseAlarm);

resp=plot_A_DT2;
rngGrid=array_Doppler_frequency.';
dopGrid=array_start_time.';
rangeIndx(1)= 70;
rangeIndx(2)= 130;
dopplerIndx(1)= 11;
dopplerIndx(2)= 379;
[columnInds,rowInds] = meshgrid(dopplerIndx(1):dopplerIndx(2),...
  rangeIndx(1):rangeIndx(2));
CUTIdx = [rowInds(:) columnInds(:)]';
%CFAR
detections_2 = cfar2D(resp,CUTIdx);
helperDetectionsMap(resp,rngGrid,dopGrid,rangeIndx,dopplerIndx,detections_2)


%% CFARCAF
% dection is a column vector, use reshape to convert into a matrix
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


%% Path Match
% Created by Eason Hua
resultList1 = [];
% Extract and process 11st column
lastCol1 = Map_1(:,dopplerIndx(1));
nonZeroRow1 = find(lastCol1 ~= 0);
if isempty(nonZeroRow1)
    resultRow1 = [101];
else
    resultRow1 = mean(nonZeroRow1);
end
resultList1 = [resultRow1];

% Extract and process 12nd to 379th columns
startCol1 = 0;
endCol1 = 0;
for i = dopplerIndx(1)+1 : dopplerIndx(2)
    if i < endCol1+1
        continue
    end

    thisCol1 = Map_1(:,i);
    nonZeroRow1 = find(thisCol1 ~= 0);
    if isempty(nonZeroRow1)
        startCol1 = i;
        % Find forwad until there exists a unempty column
        for j = i+1 : dopplerIndx(2)
            testCol1 = Map_1(:,j);
            testNonZeroRow1 = find(testCol1 ~= 0);
            if isempty(testNonZeroRow1)
                continue
            else
                endCol1 = j-1;
                endRow1 = mean(testNonZeroRow1);
                break
            end
        end
        % Assign a value for sequence, 
        % with index from startCol to endCol, 
        % and value from resultRow to endRow
        slope1 = (endRow1-resultRow1)/(endCol1-startCol1);
        for k = startCol1:endCol1
            resultList1(k-dopplerIndx(1)) =  resultRow1 + slope1*(k-startCol1);
        end
        i = endCol1+1;

    elseif length(nonZeroRow1) == 1
        % Simply take it
        resultRow1 = nonZeroRow1;

    else
        % Find the nearest index
        tempolate1 = abs(nonZeroRow1 - resultRow1);
        [miniRow1, miniCol1] = find(tempolate1==min(tempolate1));
        resultRow1 = nonZeroRow1(miniRow1(1));
    end
    resultList1 = [resultList1, resultRow1];
end



resultList2 = [];
% Extract and process 11st column
lastCol2 = Map_2(:,dopplerIndx(1));
nonZeroRow2 = find(lastCol2 ~= 0);
if isempty(nonZeroRow2)
    resultRow2 = [101];
else
    resultRow2 = mean(nonZeroRow2);
end
resultList2 = [resultRow2];

% Extract and process 12nd to 379th columns
startCol2 = 0;
endCol2 = 0;
for i = dopplerIndx(1)+1 : dopplerIndx(2)
    if i < endCol2+1
        continue
    end

    thisCol2 = Map_2(:,i);
    nonZeroRow2 = find(thisCol2 ~= 0);
    if isempty(nonZeroRow2)
        startCol2 = i;
        % Find forwad until there exists a unempty column
        for j = i+1 : dopplerIndx(2)
            testCol2 = Map_2(:,j);
            testNonZeroRow2 = find(testCol2 ~= 0);
            if isempty(testNonZeroRow2)
                continue
            else
                endCol2 = j-1;
                endRow2 = mean(testNonZeroRow2);
                break
            end
        end
        % Assign a value for sequence, 
        % with index from startCol to endCol, 
        % and value from resultRow to endRow
        slope2 = (endRow2-resultRow2)/(endCol2-startCol2);
        for k = startCol2 : endCol2
            resultList2(k-dopplerIndx(1)) =  resultRow2 + slope2*(k-startCol2);
        end
        i = endCol2+1;

    elseif length(nonZeroRow2) == 1
        % Simply take it
        resultRow2 = nonZeroRow2;

    else
        % Find the nearest index
        tempolate2 = abs(nonZeroRow2 - resultRow2);
        [miniRow2, miniCol2] = find(tempolate2==min(tempolate2));
        resultRow2 = nonZeroRow2(miniRow2(1));
    end
    resultList2 = [resultList2, resultRow2];
end


figure;
plot(resultList1, '-','Color', [1, 0.5, 0],'LineWidth', 3);
hold on;
plot(resultList2, '-', 'Color', [0.5, 1, 0], 'LineWidth',3);
hold off;

legend('show');
grid on;

%% 0
% Commented by Eason Hua
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


% 0CAF
% 0
plot_A_DT(101, :) = -1000;
plot_A_DT2(101, :) = -1000;
% % 0
% Commented by Eason Hua
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



% -Inf
% Modifided by Eason Hua
% 
% A_TD(:, 101) = 0;
% A_TD2(:, 101) = 0;
%%delete_index
delete_index = 40;
plot_A_DT = plot_A_DT(:,delete_index:end-2);
plot_A_DT2 = plot_A_DT2(:,delete_index:end-2);

%% CAF
% columnIndices ,
% Modifided by Eason Hua
maxDF1 = (resultList1 - 101) * step_dop;
maxDF2 = (resultList2 - 101) * step_dop;


%% KalmanSmoother
% Created by Eason Hua
maxDF1 = KalmanSmoother(maxDF1);
maxDF2 = KalmanSmoother(maxDF2);

%% 
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
xTar = 2;
yTar = -sqrt(2);

sensing_time = array_start_time(delete_index+10:end-10);
sensing_time = sensing_time - (delete_index -1 ) * 0.01;
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
% Modifided by Eason Hua
for i = 1:1:length(maxDF1)
    fd = [maxDF1(i);maxDF2(i)];
    F = -2*fc/c*[cos((fai_sur1(i) - fai_tx(i))/2)*cos((fai_sur1(i) + fai_tx(i))/2) , ...
        cos((fai_sur1(i) - fai_tx(i))/2)*sin((fai_sur1(i) + fai_tx(i))/2); ...
        cos((fai_sur2(i) - fai_tx(i))/2)*cos((fai_sur2(i) + fai_tx(i))/2), ...
        cos((fai_sur2(i) - fai_tx(i))/2)*sin((fai_sur2(i) + fai_tx(i))/2)];
    v_xy(i,:) = F \ fd;

    %
    % Modifided by Eason Hua
    if abs(v_xy(i,1)) > 0.4
        v_xy(i,1) = (v_xy(i-1,1) + v_xy(i-2,1) + v_xy(i-3,1))/3;
    elseif abs(v_xy(i,2)) > 0.4
        v_xy(i,2) = (v_xy(i-1,2) + v_xy(i-2,2) + v_xy(i-3,2))/3;
    end

    xtar(i+1) = xtar(i) + v_xy(i,1) * T_slide;
    ytar(i+1) = ytar(i) + v_xy(i,2) * T_slide;
    %AOA
    fai_sur1(i+1) = atan((ytar(i+1))/ (xtar(i+1) - xR1));
    fai_sur2(i+1) = atan((ytar(i+1) - yR2)/ (xtar(i+1) - xR2));
    fai_tx(i+1) = atan(ytar(i+1)/xtar(i+1));
end

%% 
% Modifided by Eason Hua
xtar = xtar(:,1:length(maxDF1));
ytar = ytar(:,1:length(maxDF2));

figure(5);
subplot(1,2,1);
plot(xtar, ytar, '-','Color', [1, 0.5, 0],'LineWidth', 3,'DisplayName', 'Estimated Trajectory');
% 
hold on;
plot(xtar(1), ytar(1), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'Initial Point');
plot(xtar(end), ytar(end), 'x', 'MarkerSize', 8, 'MarkerEdgeColor', 'b',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'End Point');
hold off;

legend('show');
title('KalmanSmoother');
grid on;

%% Rotate and Mirror
% Created by Eason Hua
coordinates = [xtar; ytar];
% 
center = mean(coordinates, 2);

% 
translated = coordinates - center;

% pi/2
theta_rotate = 7/12 * pi;
R = [cos(theta_rotate), -sin(theta_rotate); sin(theta_rotate), cos(theta_rotate)];

% 
rotated = R * translated;


% % 
% rotated = rotated + center;
% 
% % 
% subplot(1,3,2);
% plot(rotated(1, :), rotated(2, :), '-','Color', [1, 0.5, 0],'LineWidth', 3,'DisplayName', 'Rotated Trajectory');
% hold on;
% plot(rotated(1, 1), rotated(2, 1), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r',...
%     'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'Initial Point');
% plot(rotated(1, end), rotated(2, end), 'x', 'MarkerSize', 8, 'MarkerEdgeColor', 'b',...
%     'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'End Point');
% hold off;
% 
% legend('show');
% title('Rotated clockwise 7/12*pi');
% % axis equal;
% grid on;
% 
% 
% 
% % 
% rotated = rotated - center;

% Mirror horizontally
rotated = [-rotated(1,:) ; rotated(2,:)];
% 
rotated = rotated + center;


% 
subplot(1,2,2);
plot(rotated(1, :), rotated(2, :), '-','Color', [1, 0.5, 0],'LineWidth', 3,'DisplayName', 'Rotated Trajectory');
hold on;
plot(rotated(1, 1), rotated(2, 1), 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'Initial Point');
plot(rotated(1, end), rotated(2, end), 'x', 'MarkerSize', 8, 'MarkerEdgeColor', 'b',...
    'MarkerFaceColor', 'none', 'LineWidth', 2.5,'DisplayName', 'End Point');
hold off;

legend('show');
title('Rotated and Mirrored');
% axis equal;
grid on;


