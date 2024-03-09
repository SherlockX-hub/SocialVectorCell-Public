%% review 2.2;
% Simulation of cells.
% Created by Xiang Zhang, Jan., 2024.

clear;
close all;

%% parameters;
file_lists = {
    'E:\Social\15028\20210405\15028-10071\PC';
    'E:\Social\15029\20210420\15029-10068\PC';
    'E:\Social\15030\20210420\15030-10074\PC';
    'E:\Social\15032\20210407\15032-15037\PC';
    'E:\Social\15033\20210406\15033-15048\PC';
    'E:\Social\15035\20210422\15035-10075\PC';
    'E:\Social\15028\20210405\15028-10068\PC_2';
    'E:\Social\15028\20210413\15028-10071\PC_2';
    'E:\Social\15030\20210407\15030-10074\PC_2';
    'E:\Social\15029\20210411\15029-10075\PC_2';
    'E:\Social\10114\20211117\2021111703\PC';
    'E:\Social\10117\20211125\2021112503\PC';
    'E:\Social\10118\20211122\2021112203\PC';
    'E:\Social\10120\20211122\2021112203\PC';
    'E:\Social\11094\20221126\2022112601\PC';
    'E:\Social\11094\20221202\2022120201\PC';
    'E:\Social\11122\20221128\2022112802\PC';
    'E:\Social\11122\20221204\2022120402\PC';
    };

bin_size = 10;
bin_num = ceil(70 / bin_size);

% save_folder0 = 'D:\Social\lm\session';

%% simulate cells;
for file_i = 1:length(file_lists)
    file_temp = file_lists{file_i};
    sInd = strfind(file_temp, '\');
    mouse_ID = file_temp(sInd(2)+1:sInd(3)-1);
    session_ID = file_temp(sInd(4)+1:sInd(5)-1);

    %% load files;
    if strcmp(file_temp(end), '2'), load([file_temp, '\..\My_V4_Miniscope_2\ms.mat'], 'ms');
        session_ID = [session_ID, '_2']; %#ok<AGROW>
    else, load([file_temp, '\..\My_V4_Miniscope\ms.mat'], 'ms');
    end
    load([file_temp, '\..\My_WebCam\behav.mat'], 'behav');
    
    load([file_temp, '\selfPC\pos_selfPC.mat'], 'pos_selfPC');
    load([file_temp, '\otherPC\pos_otherPC.mat'], 'pos_otherPC');
    load([file_temp, '\alloPC\pos_alloPC.mat'], 'pos_alloPC');
    load([file_temp, '\egoPC\pos_egoPC.mat'], 'pos_egoPC');
    load([file_temp, '\egoPC\pos_egoPC_polar.mat'], 'pos_egoPC_polar');
    
    %% Use miniscope timestamps
% %     pos_selfPC=[behav.time behav.position{1}];
% %     pos_alloPC=[behav.time (behav.position{2}-behav.position{1}+70)/2];
    events = ms.DeconvSignals > 3 * std(ms.DeconvSignals);
    pos_ms = interp1(behav.time, pos_selfPC, ms.time);
    pos_ms_allo=interp1(behav.time, pos_alloPC, ms.time);
    
    %% Use behavioral timestamps
% %     pos_ms=[behav.time behav.position{1}];    
% %     pos_ms_allo=[behav.time (behav.position{2}-behav.position{1}+70)/2];
%     deconv = interp1(ms.time,ms.DeconvSignals,behav.time);
%     events = deconv> 3 * nanstd(deconv);
%     pos_ms = pos_selfPC;

    %% all maps;
    pos_ms(:, 4) = ceil(pos_ms(:,3) / bin_size) + (ceil(pos_ms(:,2) / bin_size) - 1) * bin_num; % pos index;
    pos_map = arrayfun(@(x) sum(pos_ms(:,4) == x), 1:bin_num^2)';
    pos_map = reshape(pos_map, bin_num,bin_num);
    
    [~,smap,map] = map_helper(pos_ms(:,2),pos_ms(:,3),events,bin_size);
    msmap = cell2mat(cellfun(@(z) z(:),smap,'UniformOutput',0));
    
    %% simulate cells;
    events_sim = zeros(size(events,1), size(events,2)*10); %zeros(size(events)); % 
    for cell_ii = 1:ms.numNeurons*10
        cell_i = cell_ii - ceil(cell_ii / ms.numNeurons - 1) * ms.numNeurons;

        k = events(:,cell_i);
        k_sim = rand(size(k));
        
        for pos_i = 1:max(pos_ms(:, 4))
            pos_temp = find(pos_ms(:, 4) == pos_i);
            if isempty(pos_temp), continue; end
            
            %% former method 4 - method five: keep event order;
            if numel(pos_temp) <= 300, k_shift = round(numel(pos_temp)/2); % 5 sec;
            else, k_shift = 150 + randi(numel(pos_temp) - 300); end
            k_sim(pos_temp) = circshift(k(pos_temp), k_shift);

        end
        k_sim(k_sim~=1) = 0;
        events_sim(:, cell_ii) = k_sim;
    end
    
    save_folder = ['D:\Social\5_10_10\session', num2str(file_i)];
    if ~exist(save_folder, 'dir'), mkdir(save_folder); end
    save([save_folder, '\events_sim.mat'], 'events_sim');
end