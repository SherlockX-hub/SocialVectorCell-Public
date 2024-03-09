%% PCDetect.m
% This function is used to detect place cells in one session.
%
% Input:
%       PC_name: a 1*2 string, the type name of place cells in another form;
%                 e.g. ["PlaceCells", "PlaceFields"];
%       DeconvSignals: a m*n matrix, n represents the number of cells;
%       calcium_time: a m*1 vector, the timestamp of calcium imaging;
%       pos: a p*3 matrix, [time x y];
%       pos_modified: a p*3 matrix, [time x y], under speed filter;
%       SFP: a width*height*cellNum matrix, spatial footprint;
%       color_PC: a 3*1 vector, [R G B], R,G,B < 1;
%       save_folder: a string, the path that save files;
%       fig_fmt: a n*1 string, the saving format of figures;
%       p: a struct, containing all necessory parameters;
%          p.datatime, p.binWidth, p.smooth, p.binMinTime, p.maxGap,
%          p.limits, p.blank, p.blank_fig, p.draw_fig (bool, draw fig or
%          not)
%
% Output:
%

% v2: 2022 Sept.
% Make some improvements;
% Use smoothed ratemap to calculate information content;
% Add half-stability threshold of 0.3.

% Created by Xiang Zhang, 2021.

function PCDetect(PC_name, DeconvSignals, calcium_time, pos, pos_modified, ...
        SFP, color_PC, save_folder, fig_fmt, p)
    %% initialization;
    if isempty(PC_name), PC_name = ["PlaceCells", "PlaceFields"];
    elseif length(PC_name) == 1, PC_name(2) = "PlaceFields"; end
    
    try draw_fig = p.draw_fig; catch, draw_fig = 1; end
    
    disp(['Start detecting ', PC_name{1}, '...']);
    
    % make folders;
    if ~exist(save_folder,'dir'), mkdir(save_folder); end
    mkdir(strcat(save_folder, '\maps'));
    mkdir(strcat(save_folder, '\shuffles'));
    mkdir(strcat(save_folder, '\stability'));
    if draw_fig
        mkdir(strcat(save_folder, '\cells'));
        mkdir(strcat(save_folder, '\ratemaps'));
        mkdir(strcat(save_folder, '\cellpaths'));
        mkdir(strcat(save_folder, '\merge'));
    end
    
    disp('Folds are made.');
    
    % information;
    numNeurons = size(DeconvSignals,2);
    information_all = struct('rate',cell(numNeurons,1), 'content',cell(numNeurons,1));
    threshold_all = nan(numNeurons,1);
    PC = nan(numNeurons,1);
    PF = struct('PC',cell(numNeurons,1), 'num',cell(numNeurons,1), ...
        'size',cell(numNeurons,1), 'spike_num', cell(numNeurons,1), ...
        'peak',cell(numNeurons,1), 'peakXY',cell(numNeurons,1), 'meanRate',cell(numNeurons,1));
    
    % stability;
    % [splitDataArray, ~, ~, ~, ~] = splitData(pos);
    [splitDataArray_modified, ~, ~, ~, ~] = splitData(pos_modified);
    
    %% main function;
    parfor i = 1:numNeurons
        k = find(DeconvSignals(:,i) > 3 * std(DeconvSignals(:,i)));
        
        if ~isempty(k) && length(k) > 10 % identify they are all 0;
            spike_pos = calcium_time(k); %#ok<PFBNS> % the time of spike
            spike_pos = spikePos(spike_pos, pos(:,1), pos_modified); %#ok<PFBNS>
            if isempty(spike_pos)
                information_all(i).rate = nan;
                information_all(i).content = nan;
                continue;
            end
            
            % information;
            information_all(i) = infoCal(i, pos_modified, spike_pos, p, save_folder);
            % p.limits =  [0 behav.trackLength 0 behav.trackLength];
            
            % shuffle;
            threshold_all(i) = cellShuffle(i, information_all(i).content, calcium_time, k, ...
                pos_modified, p, save_folder, fig_fmt);
            
            % stability;
            map_stability = cellStability(i, spike_pos, ...
                pos_modified, splitDataArray_modified, p, save_folder);
            
            % place cell;
            if information_all(i).content > threshold_all(i) && map_stability.half_correlation > 0.3
                    % || map_stability.time_correlation > 0.3)
                [PC(i), PF(i)] = isPCProcess(i, pos, pos_modified, p, information_all(i).content, ...
                    save_folder, draw_fig, SFP, fig_fmt, color_PC); %#ok<PFOUS>
            end
        else
            information_all(i).rate = nan;
            information_all(i).content = nan;
        end
    end
    disp('Information contents are calculated and place cells are found.');
    
    %% save files and draw figures;
    PC(isnan(PC)) = [];
    information_all = changeStruct(information_all);
    saveFiles(save_folder, information_all, threshold_all);
    eval(strcat(PC_name(1), '= PC;'));
    save_path = strcat(save_folder, filesep, PC_name(1), '.mat');
    save(save_path, PC_name(1));
    eval(strcat(PC_name(2), '= PF;'));
    save_path = strcat(save_folder, filesep, PC_name(2), '.mat');
    save(save_path, PC_name(2));
    
    if draw_fig
        PCPlot(PC, PC_name{1}, SFP, save_folder, fig_fmt, color_PC);
        close();
    end
    
    %% show;
    % fprintf('session %s has been finished calculating.\n',dir_name{1});
    disp(strcat('Finished detecting ', 32, PC_name{1}, 32, 'and files are saved.'));
    
end

%% functions;
% split data in two ways;
function [splitDataArray, pos_first, pos_second, pos_odd, pos_even] = splitData(pos)
    % spilt data in four ways;
    % 1 Row:    First half of data (Half and half type split)
    % 2 Row:    Second half of data (Half and half type split)
    % 3 Row:    First half of data (Binned type split) (1 min and 1 min)
    % 4 Row:    Second half of data (Binned type split) (1 min and 1 min)
    % 1 Col:    Post
    % 2 Col:    Posx
    % 3 Col:    Posy
    splitDataArray = cell(4,3);
    
    duration = pos(end,1); % msec; - pos(1,1);
    half_ind = find(pos(:,1) <= duration/2); % + pos(1,1));
    splitDataArray{1,1} = pos(half_ind,1);
    splitDataArray{1,2} = pos(half_ind,2);
    splitDataArray{1,3} = pos(half_ind,3);
    pos_first = pos(half_ind,:);
    % pos_first = [spiltDataArray{1,1} spiltDataArray{1,2}
    % spiltDataArray{1,3}]; % = cell2mat(spiltDataArray(1,:));
    
    half_ind = find(pos(:,1) > duration/2); % + pos(1,1));
    splitDataArray{2,1} = pos(half_ind,1);
    splitDataArray{2,2} = pos(half_ind,2);
    splitDataArray{2,3} = pos(half_ind,3);
    pos_second = pos(half_ind,:);
    % pos_second = [spiltDataArray{2,1} spiltDataArray{2,2}
    % spiltDataArray{2,3}]; % = cell2mat(spiltDataArray(2,:));
    
    splitDataArray(3:4, :) = splitOddEven(pos);
    pos_odd = [splitDataArray{3,1} splitDataArray{3,2} splitDataArray{3,3}]; % = cell2mat(splitDataArray(3,:));
    pos_even = [splitDataArray{4,1} splitDataArray{4,2} splitDataArray{4,3}]; % = cell2mat(splitDataArray(4,:));
end

% split data in odd and even times;
function splitDataArray = splitOddEven(pos)
    duration = pos(end,1); %  - pos(1,1);
    [x3,x4,y3,y4,t3,t4] = deal(zeros(size(pos,1),1));
    [totSamp3,totSamp4] = deal(0);
    
    nBins = ceil(duration / 60000);
    t_start = 0; % pos(1,1);
    t_stop = t_start + 60000;
    for spilt_i = 1:nBins
        ind = find(pos(:,1) >= t_start & pos(:,1) < t_stop);
        samps = length(ind);
        if samps > 0
            if mod(spilt_i,2) % odd;
                x3(totSamp3+1:totSamp3+samps) = pos(ind,2);
                y3(totSamp3+1:totSamp3+samps) = pos(ind,3);
                t3(totSamp3+1:totSamp3+samps) = pos(ind,1);
                totSamp3 = totSamp3 + samps;
            else % even;
                x4(totSamp4+1:totSamp4+samps) = pos(ind,2);
                y4(totSamp4+1:totSamp4+samps) = pos(ind,3);
                t4(totSamp4+1:totSamp4+samps) = pos(ind,1);
                totSamp4 = totSamp4 + samps;
            end
        end
        t_start = t_stop;
        t_stop = t_start + 60000;
    end
    
    splitDataArray{1,2} = x3(1:totSamp3);
    splitDataArray{1,3} = y3(1:totSamp3);
    splitDataArray{1,1} = t3(1:totSamp3);
    splitDataArray{2,2} = x4(1:totSamp4);
    splitDataArray{2,3} = y4(1:totSamp4);
    splitDataArray{2,1} = t4(1:totSamp4);
end

function information = infoCal(cell_i, pos, spike_pos, p, save_folder)
    % calculate rate map;
    map = map_new(pos, spike_pos(:,1), 'datatime', p.datatime, 'binWidth', p.binWidth, 'smooth', p.smooth, ...
        'minTime', p.binMinTime, 'blanks', p.blank, 'maxGap', p.maxGap, 'limits', p.limits);
    
    % spike position;
    map.spike_pos = spike_pos;
    save(strcat(save_folder, '\maps\map',num2str(cell_i),'.mat'),'map');
    
    % information about map
    [information, ~, ~] = analyses.mapStatsPDF(map);
end

function threshold = cellShuffle(cell_i, information_content, calcium_time, k, pos_modified, p, save_folder, fig_fmt)
    information_random = nan(1001,1);
    information_random(1) = information_content;
    if strcmpi(p.datatime, 'msec')
        time_random = randi(round(calcium_time(end) - calcium_time(1)),1000,1);
    else
        time_random = randi(round(1000*(calcium_time(end) - calcium_time(1))),1000,1) / 1000;
    end
    
    for m = 1:1000
        spike_pos_random = [];
        spike_pos_random(:,1) = calcium_time(1) + mod(calcium_time(k) - calcium_time(1) + time_random(m), ...
            calcium_time(end) - calcium_time(1));
        % calculate rate map;
        map = map_new(pos_modified, spike_pos_random(:,1), 'datatime', p.datatime, 'binWidth', p.binWidth, 'smooth', p.smooth, ...
            'minTime', p.binMinTime, 'blanks', p.blank, 'maxGap', p.maxGap, 'limits', p.limits);
        [information, ~, ~] = analyses.mapStatsPDF(map);
        information_random(m) = information.content;
    end
    save([save_folder, '\shuffles\information_random', num2str(cell_i), '.mat'], 'information_random');
    
    threshold = prctile(information_random, 95);
    drawShufflePlot(information_content, cell_i, information_random, threshold, save_folder, fig_fmt);
end

function drawShufflePlot(information, cell_i, information_random, threshold, save_folder, fig_fmt)
    figure;
    histogram(information_random);
    hold on;
    plot([threshold,threshold],get(gca,'Ylim'),'m--');
    scatter(information,0); % ,[],'r','filled');
    hold off;
    for i = 2%1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '\shuffles\shuffle', num2str(cell_i), fig_fmt(i)));
    end
    close();
end

function map_stability = cellStability(cell_i, spike_pos, ...
        pos_modified, splitDataArray_modified, p, save_folder)
    map_stability = struct;
    map_stability = cellStability_half(map_stability, spike_pos, ...
        pos_modified(:,1), pos_modified, splitDataArray_modified, p);
    map_stability = cellStability_time(map_stability, spike_pos, splitDataArray_modified, p);
    
    save([save_folder, '\stability\map_stability', num2str(cell_i), '.mat'], 'map_stability');
end

% half stability;
function map_stability = cellStability_half(map_stability, spike_pos, ...
        behav_time, pos_modified, splitDataArray_modified, p)
    spike_time = spike_pos(:,1);
    spike_pos_first = spike_time(spike_time <= (behav_time(end) - behav_time(1)) / 2 + behav_time(1));
    spike_pos_second = spike_time(spike_time > (behav_time(end) - behav_time(1)) / 2 + behav_time(1));
    
    spike_pos_first = spikePos(spike_pos_first, behav_time, pos_modified);
    spike_pos_second = spikePos(spike_pos_second, behav_time, pos_modified);
    
    pos_first_modified = cell2mat(splitDataArray_modified(1,:));
    pos_second_modified = cell2mat(splitDataArray_modified(2,:));
    
    if ~isempty(spike_pos_first) && ~isempty(spike_pos_second)
        map_first = map_new(pos_first_modified, spike_pos_first(:,1), ...
            'datatime', p.datatime, 'binWidth', p.binWidth, 'smooth', p.smooth, ...
            'minTime', p.binMinTime, 'maxGap', p.maxGap, 'limits', p.limits);
        map_first.spike_pos = spike_pos_first;
        map_second = map_new(pos_second_modified, spike_pos_second(:,1), ...
            'datatime', p.datatime, 'binWidth', p.binWidth, 'smooth', p.smooth, ...
            'minTime', p.binMinTime, 'maxGap', p.maxGap, 'limits', p.limits);
        map_second.spike_pos = spike_pos_second;
        
        % correlation;
        r_half = zeroLagCorrelation(map_first.z,map_second.z);
    else
        map_first = [];
        map_second = [];
        r_half = nan;
    end
    map_stability.map_first = map_first;
    map_stability.map_second = map_second;
    map_stability.half_correlation = r_half;
end

% time binned stability;
function map_stability = cellStability_time(map_stability, spike_pos, splitDataArray_modified, p)
    
    splitDataArray_spike = splitOddEven(spike_pos);
    
    spike_pos_odd = cell2mat(splitDataArray_spike(1,:));
    spike_pos_even = cell2mat(splitDataArray_spike(2,:));
    
    if ~isempty(spike_pos_odd) && ~isempty(spike_pos_even)
        pos_odd_modified = cell2mat(splitDataArray_modified(3,:));
        pos_even_modified = cell2mat(splitDataArray_modified(4,:));
        
        map_odd = map_new(pos_odd_modified, spike_pos_odd(:,1), 'datatime', p.datatime, ...
            'binWidth', p.binWidth, 'smooth', p.smooth, 'minTime', p.binMinTime, ...
            'maxGap', p.maxGap, 'limits', p.limits);
        map_odd.spike_pos = spike_pos_odd;
        map_even = map_new(pos_even_modified, spike_pos_even(:,1), 'datatime', p.datatime, ...
            'binWidth', p.binWidth, 'smooth', p.smooth, 'minTime', p.binMinTime, ...
            'maxGap', p.maxGap, 'limits', p.limits);
        map_even.spike_pos = spike_pos_even;
        
        % correlation;
        r_time = zeroLagCorrelation(map_odd.z,map_even.z);
    else
        map_odd = [];
        map_even = [];
        r_time = nan;
    end
    map_stability.map_odd = map_odd;
    map_stability.map_even = map_even;
    map_stability.time_correlation = r_time;
end

function [PC, PF] = isPCProcess(cell_i, pos, pos_modified, p, info_content, save_folder, draw_fig, SFP, fig_fmt, color_PC)
    PC = cell_i;
    
    load(strcat(save_folder, '\maps\map',num2str(cell_i),'.mat'), 'map');
    % place field;
    PF = PFall(cell_i, map, p, pos_modified);
    
    if draw_fig
        % ratemap;
        ratemapPlot(map, cell_i, save_folder, fig_fmt);
        
        % cell path;
        cellpathPlot(cell_i, pos, p.limits, map.spike_pos, save_folder, fig_fmt);
        
        % merge plots;
        cellMerge(cell_i, map, info_content, pos, p.limits, save_folder, fig_fmt);
        
        % plot cells;
        cellPlot(SFP, cell_i, save_folder, fig_fmt, color_PC);
        close all;
    end
end

function spike_pos = spikePos(spike_pos, time, pos)
    if isempty(spike_pos)
        return;
    end
    for j = 1:size(spike_pos,1)
        if time(end) > spike_pos(j,1)
            msFrameTomsBehav = find(time > spike_pos(j,1));
            spike_pos(j,2) = pos(msFrameTomsBehav(1),2); % figure out the first position X after spike time;
            spike_pos(j,3) = pos(msFrameTomsBehav(1),3); % figure out the first position Y after spike time;
        else
            spike_pos(j,2) = nan;
            spike_pos(j,3) = nan;
        end
    end
    spike_pos(isnan(spike_pos(:,2)),:) = [];
end

function ratemapPlot(map_figure, cell_i, save_folder, fig_fmt)
    figure;
    plot.colorMap(map_figure.z);
    axis equal;
    axis off;
    title(sprintf('%s%3.2f','Rate Map. Peak = ',nanmax(nanmax(map_figure.z))));
    for i = 1%:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '\ratemaps\cell', num2str(cell_i), 'ratemap', fig_fmt(i)));
    end
end

function cellpathPlot(cell_i, pos, behavLimit, spike_pos, save_folder, fig_fmt)
    figure;
    plot(pos(:,2),pos(:,3), 'color',[0.5,0.5,0.5],'LineWidth',1);
    axis equal;
    xlimit = get(gca,'xlim');
    ylimit = get(gca,'ylim');
    if xlimit(2)-xlimit(1) < behavLimit(2) - behavLimit(1)
        xlim(behavLimit(1:2));
    end
    if ylimit(2)-ylimit(1) < behavLimit(2) - behavLimit(1)
        ylim(behavLimit(1:2));
    end
    hold on;
    sz = 20;
    scatter(spike_pos(:,2),spike_pos(:,3),sz,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
    hold off;
    for i = 2%:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '\cellpaths\cell', num2str(cell_i), 'path', fig_fmt(i)));
    end
    
    axis off;
    for i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '\cellpaths\cell', num2str(cell_i), 'path-noAxis', fig_fmt(i)));
    end
end

function cellPlot(SFP, cell_i, save_folder, fig_fmt, color_PC)
    % figure;
    % imagesc(permute(SFP(cell_i,:,:),[2 3 1]));
    % for i = 1:length(fig_fmt)
    %     saveas(gcf,strcat(save_folder, '\cells\cell', num2str(cell_i), fig_fmt(i)));
    % end
    
    markPartSFP_color_v3(SFP, SFP(cell_i,:,:), color_PC);
    for i = 1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '\cells\cell', num2str(cell_i), 'inAll', fig_fmt(i)));
    end
end

function cellMerge(cell_i, map, info_content, pos, behavLimit, save_folder, fig_fmt)
    figure('Position', [500 300 1000 400]);
    
    subplot(1,2,1);
    plot.colorMap(map.z, 'bar','off');
    axis equal;
    axis off;
    title(sprintf('Peak: %3.2f, Info: %3.2f', nanmax(nanmax(map.z)), info_content));
    
    subplot(1,2,2);
    plot(pos(:,2),pos(:,3), 'color',[0.5,0.5,0.5],'LineWidth',1);
    axis equal;
    xlimit = get(gca,'xlim');
    ylimit = get(gca,'ylim');
    if xlimit(2)-xlimit(1) < behavLimit(2) - behavLimit(1)
        xlim(behavLimit(1:2));
    end
    if ylimit(2)-ylimit(1) < behavLimit(2) - behavLimit(1)
        ylim(behavLimit(1:2));
    end
    hold on;
    sz = 20;
    scatter(map.spike_pos(:,2),map.spike_pos(:,3),sz,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
    hold off;
    axis off;
    
    for i = 2%1:length(fig_fmt)
        saveas(gcf,strcat(save_folder, '\merge\cell', num2str(cell_i), fig_fmt(i)));
    end
end

function PF = PFall(cell_i, map_figure, p, pos)
    [~, fields] = analyses.placefield(map_figure, 'binWidth', p.binWidth, 'pos', pos, 'minPeak', 0.5, 'threshold', 0.6);
    
    field_size = size(length(fields),1);
    field_peak = size(length(fields),1);
    field_peakCoordinate = size(length(fields),2);
    field_meanRate = size(length(fields),1);
    field_spike_num = size(length(fields),1);
    for field_i = 1:length(fields)
        field_size(field_i,1) = fields(field_i).size;
        field_peak(field_i,1) = fields(field_i).peak;
        field_peakCoordinate(field_i,1) = fields(field_i).x * p.binWidth;
        field_peakCoordinate(field_i,2) = fields(field_i).y * p.binWidth;
        field_meanRate(field_i,1) = fields(field_i).meanRate;
        field_spike_num(field_i,1) = nansum(map_figure.Nspikes(fields(field_i).PixelIdxList));
    end
    
    PF.PC = cell_i;
    PF.num = length(fields);
    PF.size = field_size;
    PF.peak = field_peak;
    PF.peakXY = field_peakCoordinate;
    PF.meanRate = field_meanRate;
    PF.spike_num = field_spike_num;
    if PF.num == 0
        PF.size = 0;
        PF.peak = 0;
        PF.peakXY = [0,0];
        PF.meanRate = 0;
        PF.spike_num = 0;
    end
end

function PCPlot(PC, PC_name, SFP, save_folder, fig_fmt, color_PC)
    if ~isempty(PC)
        SFP_PC = SFP(PC,:,:);
        eval(['SFP_' PC_name ' = SFP_PC;']);
        figure;
        imagesc(max(permute(SFP_PC,[2 3 1]),[],3));
        save(strcat(save_folder, '\SFP_', PC_name, '.mat'),['SFP_' PC_name]);
        %     if strcmpi(PC_name, 'selfPC')
        %         PC_name2 = 'placecells';
        %         saveas(gcf,strcat(save_folder, '\', PC_name2, '.fig'));
        %     else
        %         saveas(gcf,strcat(save_folder, '\', PC_name, '.fig'));
        %     end
        
        markPartSFP_color_v3(SFP, SFP_PC, color_PC);
        for i = 1:length(fig_fmt)
            saveas(gcf,strcat(save_folder, '\', PC_name, 'inAll', fig_fmt(i)));
        end
        close all;
    end
end

function information_all = changeStruct(information_all)
    information_all2.content = cat(1,information_all.content);
    information_all2.rate = cat(1,information_all.rate);
    information_all = information_all2;
end

function saveFiles(save_folder, information_all, threshold_all) % , PlaceCells, PlaceFields)
    save(strcat(save_folder, '\information_all.mat'),'information_all');
    save(strcat(save_folder, '\threshold_all.mat'),'threshold_all');
    % save(strcat(save_folder, '\PlaceCells.mat'),'PlaceCells');
    % save(strcat(save_folder, '\PlaceFields.mat'),'PlaceFields');
end