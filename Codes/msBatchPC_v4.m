%% msBatchPC_v4_ZX.m
% This code is used to find four kinds of place cells recorded by miniscope.
% v3: Use data from behav file.
% Compute multiple cells in the same time instead of shuffle.
% Considering that there are more than one calcium imaging recording in one
% session.
%
% v4: intergrate PC analysis into a function.
% Egocentric place cells will be analyzed in a Cartesian coordinates;
% Change ratemap p.smooth to 2;
% Change place field minPeak to 0.5, threshold to 0.6;
% Add half-stability threshold of 0.3.

% Created by Xiang Zhang, 2021.

% warning off all;
% clc;
clear;
tic;

%% paths of recordings;
info = readtable('E:\Social\Dir.xls.xlsx');
file_list = info.path;
% file_list = info.path(strcmp(info.sessionType, 'chase'));

%% coding path;
addpath(genpath('G:\ZX\Codes\BNT-master'));

%% selective;
doSelfPCProcessInit = 1;
doOtherPCProcessInit = 1;
doAlloPCProcessInit = 1;
doEgoPCProcessInit = 1;

%% general parameters;
% ratemap;
p.datatime = 'msec';
p.binWidth = 2; % cm
p.binWidthInit = 2; % cm
p.smooth = 2; % bins
p.binMinTime = 0;
p.maxGap = 0.300;
% p.type = 'circular';1
% p.minCoverage = 80;
% p.minNumSpikes = 30; % if coverage >= p.minCoverage && numSpikesSpeedFiltered >= p.minNumSpikes
p.lowSpeedThreshold = 0.0025; % cm / ms;
p.blank = 'on';

% cell color;
if doSelfPCProcessInit, color_selfPC = [213,129,119] / 255; end
if doOtherPCProcessInit, color_otherPC = [240,214,116] / 255; end
if doAlloPCProcessInit, color_alloPC = [162,197,157] / 255; end
if doEgoPCProcessInit, color_egoPC = [170,174,211] / 255; end

% figure;
figFmt = [".fig",".png"];

%% main function;
for folder_i = 1:length(file_list)%length(file_list):-1:1
    dir_name = file_list{folder_i};
    cd(dir_name);
    if strcmpi(dir_name(1), 'E') || strcmpi(dir_name(1), 'F')
        sInd = strfind(dir_name,'\');
        mouseID = dir_name(sInd(2)+1:sInd(3)-1);
        session = dir_name(sInd(3)+1:end);
        sInd2 = strfind(session,'\');
        if ~isempty(sInd2)
            session(sInd2) = '_';
        end
    elseif strcmpi(dir_name(1), 'D')
        sInd = strfind(dir_name,'\');
        mouseID_session = dir_name(sInd(end)+1:end);
        sInd2 = strfind(mouseID_session,'_');
        mouseID = mouseID_session(1:sInd2(1)-1);
        session = mouseID_session(sInd2(1)+1:end);
    end
    disp(['Start calculating session: ', mouseID, 32, session]);
    
    ms_type = 1; % v4;
    miniscopeRecording = dir('My_V4_Miniscope*');
    if isempty(miniscopeRecording)
        miniscopeRecording = dir('ms.mat');
        ms_type = 2; % v3;
        if isempty(miniscopeRecording), continue; end
    end
    
    for ms_i = 1:length(miniscopeRecording)
        ms_dir = miniscopeRecording(ms_i).name;
        
        if ms_type == 1
            save_folder = ['PC', erase(ms_dir, 'My_V4_Miniscope')];
        elseif ms_type == 2
            sInd = strfind(ms_dir, '.');
            save_folder = ['PC', erase(ms_dir(1:sInd(1)-1), 'ms')];
        end
        
        %% verify place cells version;
        v_file = strcat(save_folder, filesep, 'version.txt');
        skipAnalyze = 0;
        if exist(v_file, 'file')
            fid = fopen(v_file, 'r');
            while ~feof(fid)
                fline = fgetl(fid);
                ver_PC = strfind(fline, 'Find place cells.');
                if ver_PC, skipAnalyze = 1; end
            end
            fclose(fid);
        end
        
        if skipAnalyze, continue; end
        
        %% load information files;
        if ms_type == 1, load([ms_dir, '\ms.mat'], 'ms');
        elseif ms_type == 2, load('ms.mat', 'ms'); end
        
        if ~isfield(ms, 'numNeurons'), continue; end
        
        if ms_type == 1, load([ms_dir, '\SFP.mat'], 'SFP');
        elseif ms_type == 2, load('SFP.mat', 'SFP'); end
        
        try load('My_WebCam\behav.mat', 'behav');
        catch
            try load('behav.mat'); catch, continue; end
        end
        
        if ~isfield(behav,'position'), continue; end
        
        if isfield(behav, 'wirelessStart') && strcmpi(ms_dir,'My_V4_Miniscope_2') && ~isfield(ms,'timeOriginal')
            ms.timeOriginal = ms.time;
            ms.time = ms.time + behav.time(behav.wirelessStart);
            save([ms_dir, '\ms.mat'], 'ms');
        end
        
        doSelfPCProcess = doSelfPCProcessInit;
        doOtherPCProcess = doOtherPCProcessInit;
        doAlloPCProcess = doAlloPCProcessInit;
        doEgoPCProcess = doEgoPCProcessInit;
        if length(behav.position) == 1
            doOtherPCProcess = 0;
            doAlloPCProcess = 0;
            doEgoPCProcess = 0;
        end
        if ~isfield(behav, 'hdDir')
            doEgoPCProcess = 0;
        end
        
        % save folders;
        if ~exist(save_folder, 'dir'), mkdir(save_folder); end
        
        save_folder_selfPC = [save_folder, filesep, 'selfPC'];
        save_folder_otherPC = [save_folder, filesep, 'otherPC'];
        save_folder_alloPC = [save_folder, filesep, 'alloPC'];
        save_folder_egoPC = [save_folder, filesep, 'egoPC'];
        
        if doSelfPCProcess, if exist([save_folder_selfPC filesep 'selfPC.mat'], 'file'), doSelfPCProcess = 0; end, end
        if doOtherPCProcess, if exist([save_folder_otherPC filesep 'otherPC.mat'], 'file'), doOtherPCProcess = 0; end, end
        if doAlloPCProcess, if exist([save_folder_alloPC filesep 'alloPC.mat'], 'file'), doAlloPCProcess = 0; end, end
        if doEgoPCProcess, if exist([save_folder_egoPC filesep 'egoPC.mat'], 'file'), doEgoPCProcess = 0; end, end
        
        if ~doSelfPCProcess && ~doOtherPCProcess && ~doAlloPCProcess && ~doEgoPCProcess, continue; end
        
        %% parameters;
        % position;
        pos_selfPC = [];
        pos_selfPC(:,1) = behav.time;
        pos_selfPC(:,2:3) = behav.position{ms_i}; % time; position X; position Y
        
        % speed;
        if p.lowSpeedThreshold > 0
            speed = speed2D(pos_selfPC(:,2),pos_selfPC(:,3),pos_selfPC(:,1));
            lowSpeed_index = speed < p.lowSpeedThreshold;
        end
        
        disp('Initialization finished.');
        
        %% self's place cells;
        % parameters;
        if doSelfPCProcess
            mkdir(save_folder_selfPC);
            % positions;
            save([save_folder_selfPC, '\pos_selfPC.mat'], 'pos_selfPC');
            
            pos_selfPC_modified = pos_selfPC;
            if p.lowSpeedThreshold > 0
                pos_selfPC_modified(lowSpeed_index, 2:3) = NaN;
            end
            save([save_folder_selfPC, '\pos_selfPC_modified.mat'], 'pos_selfPC_modified');
            
            if sum(isnan(pos_selfPC_modified(:,2))) == size(pos_selfPC_modified,1)
                doSelfPCProcess = 0;
                doAlloPCProcess = 0;
                doEgoPCProcess = 0;
            end
        end
        
        if doSelfPCProcess
            % ratemap;
            p.limits = [0 behav.trackLength 0 behav.trackLength];
            p.binWidth = p.binWidthInit;
            p.blank = 'on';
            if ~isfield(behav,'shape')
                p.blank = 'off';
            else
                if behav.shape == 1, p.blank = 'off'; end
            end
            
            PCDetect(["selfPC", "selfPF"], ...
                ms.DeconvSignals, ms.time, pos_selfPC, pos_selfPC_modified, SFP, ...
                color_selfPC, save_folder_selfPC, figFmt, p);
        end
        
        %% others' place cells;
        if doOtherPCProcess
            mkdir(save_folder_otherPC);
            % positions;
            pos_otherPC = [];
            pos_otherPC(:,1) = behav.time;
            pos_otherPC(:,2:3) = behav.position{3 - ms_i};
            if max(pos_otherPC(:,2)) - min(pos_otherPC(:,2)) < 10 && ...
                    max(pos_otherPC(:,3)) - min(pos_otherPC(:,3)) < 10
                continue;
            end
            save([save_folder_otherPC, '\pos_otherPC.mat'], 'pos_otherPC');
            
            pos_otherPC_modified = pos_otherPC;
            if p.lowSpeedThreshold > 0
                pos_otherPC_modified(lowSpeed_index, 2:3) = NaN;
            end
            save([save_folder_otherPC, '\pos_otherPC_modified.mat'], 'pos_otherPC_modified');
            
            if sum(isnan(pos_otherPC_modified(:,2))) == size(pos_otherPC_modified,1)
                doOtherPCProcess = 0;
                doAlloPCProcess = 0;
                doEgoPCProcess = 0;
            end
        end
        
        if doOtherPCProcess
            % ratemap;
            p.limits = [0 behav.trackLength 0 behav.trackLength];
            p.binWidth = p.binWidthInit;
            p.blank = 'on';
            if ~isfield(behav,'shape')
                p.blank = 'off';
            else
                if behav.shape == 1, p.blank = 'off'; end
            end
            
            PCDetect(["otherPC", "otherPF"], ...
                ms.DeconvSignals, ms.time, pos_otherPC, pos_otherPC, SFP, ...
                color_otherPC, save_folder_otherPC, figFmt, p);
        end
        
        %% allocentric place cells; self-vector place cell;
        if doAlloPCProcess
            mkdir(save_folder_alloPC);
            % positions;
            pos_alloPC = [];
            pos_alloPC(:,1) = behav.time;
            pos_alloPC(:,2:3) = behav.position{3 - ms_i} - behav.position{ms_i};
            save([save_folder_alloPC, '\pos_alloPC.mat'], 'pos_alloPC');
            
            pos_alloPC_modified = pos_alloPC;
            if p.lowSpeedThreshold > 0
                pos_alloPC_modified(lowSpeed_index, 2:3) = NaN;
            end
            save([save_folder_alloPC, '\pos_alloPC_modified.mat'], 'pos_alloPC_modified');
            
            % ratemap;
            p.binWidth = 2 * p.binWidthInit; % cm
            p.limits = [-behav.trackLength behav.trackLength -behav.trackLength behav.trackLength];
            p.blank = 'on';
            
            % PC;
            PCDetect(["alloPC", "alloPF"], ...
                ms.DeconvSignals, ms.time, pos_alloPC, pos_alloPC, SFP, ...
                color_alloPC, save_folder_alloPC, figFmt, p);
            
            % recover;
            % p.binWidth = p.binWidthInit; % cm
        end
        
        %% egocentric place cells;
        if doEgoPCProcess
            mkdir(save_folder_egoPC);
            % positions;
            pos_egoPC = [];
            pos_egoPC(:,1) = behav.time;
            vectorTwo = behav.position{3 - ms_i} - behav.position{ms_i};
            pos_egoPC(:,2) = mod(atan2(vectorTwo(:,2), vectorTwo(:,1)) - behav.hdDir{ms_i} * pi / 180, 2*pi); % direction;
            pos_egoPC(:,3) = vecnorm(vectorTwo,2,2); % distance;
            pos_egoPC_polar = pos_egoPC;
            [pos_egoPC(:,2), pos_egoPC(:,3)] = pol2cart(pos_egoPC(:,2), pos_egoPC(:,3));
            save([save_folder_egoPC, '\pos_egoPC.mat'], 'pos_egoPC');
            save([save_folder_egoPC, '\pos_egoPC_polar.mat'], 'pos_egoPC_polar');
            
            pos_egoPC_modified = pos_egoPC;
            if p.lowSpeedThreshold > 0
                pos_egoPC_modified(lowSpeed_index, 2:3) = NaN;
            end
            save([save_folder_egoPC, '\pos_egoPC_modified.mat'], 'pos_egoPC_modified');
            
            % ratemap;
            if isfield(behav,'shape')
                if behav.shape == 2 || behav.shape == 3
                    p.limits = [-behav.trackLength behav.trackLength -behav.trackLength behav.trackLength];
                else
                    p.limits = [-sqrt(2)*behav.trackLength sqrt(2)*behav.trackLength ...
                        -sqrt(2)*behav.trackLength sqrt(2)*behav.trackLength];
                end
            else
                p.limits = [-sqrt(2)*behav.trackLength sqrt(2)*behav.trackLength ...
                    -sqrt(2)*behav.trackLength sqrt(2)*behav.trackLength];
            end
            
            % ratemap;
            p.binWidth = 2 * p.binWidthInit; % cm
            p.blank = 'on';
            
            % PC;
            PCDetect(["egoPC", "egoPF"], ...
                ms.DeconvSignals, ms.time, pos_egoPC, pos_egoPC, SFP, ...
                color_egoPC, save_folder_egoPC, figFmt, p);
            
            % recover;
            % p.binWidth = p.binWidthInit; % cm
        end
        
        %% cell plots;
        markAllPC_v4('PC', save_folder_selfPC, save_folder_otherPC, save_folder_alloPC, save_folder_egoPC, SFP, 0);
        for i = 1:length(figFmt)
            saveas(gcf,strcat(save_folder, filesep, 'PCinAll', figFmt(i)));
        end
        close();
        
        %% create version file;
        v_file = strcat(save_folder, filesep, 'version.txt');
        fid = fopen(v_file, 'a');
        fprintf(fid, '\n%s\n', string(datetime));
        fprintf(fid, '%s\n', 'Find place cells.');
        fclose(fid);
        
    end % end of one miniscope recording;
    
    %% show;
    % fprintf('session %s has been finished calculating.\n',dir_name{1});
    disp(strcat('session ', 32, dir_name, ' has been finished calculating and files are saved.'));
end

%% end and clear
rmpath(genpath('G:\ZX\Codes\BNT-master'));
% fclose('all');
close all;
fprintf('%s%s\n','Finished: ', datestr(now));
toc;

% clear;
% clc;

%% functions;