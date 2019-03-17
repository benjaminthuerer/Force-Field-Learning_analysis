% Load c3d files for each subject and rearrange the trial order and compute
% the relevant statistics in a table
%%
clear

s_rate = double(1000);
s_time = double(1)/s_rate;

x_pos = [0, 0.0866, 0.0866, 0, -0.0866, -0.0866, 0];
y_pos = [0.298, 0.248, 0.148, 0.098, 0.148, 0.248, 0.198];
% set study folder and add to path
mainStudyPath = uigetdir('','choose the study folder, e.g. C:\...\DAVOS');
scriptPath = [mainStudyPath '\scripts\matlab'];
cd(scriptPath);
addpath(scriptPath)

% load csv data with all relevant subject and trial information
subjDistribFile = [scriptPath '\' 'subj_distribution.csv'];
fid = fopen(subjDistribFile,'r');
C = textscan(fid,repmat('%s',1,6),'delimiter',';','CollectOutput',1);
subjDistribution = C{1};
clear C subjDistribFile
fclose(fid);

% load blocked trial order
blockedTrialOrderFile = [scriptPath '\' 'blocked_trial_order.csv'];
fid = fopen(blockedTrialOrderFile,'r');
C = textscan(fid,repmat('%s',1,18),'delimiter',';','CollectOutput',1);
blockedTrialOrder = C{1};
clear C blockedTrialOrderFile
fclose(fid);
blockedHeader = blockedTrialOrder(1,:);
blockedTrialOrder = f_sort_trials(blockedTrialOrder);

%load random trial order
randomTrialOrderFile = [scriptPath '\' 'random_trial_order.csv'];
fid = fopen(randomTrialOrderFile,'r');
C = textscan(fid,repmat('%s',1,18),'delimiter',';','CollectOutput',1);
randomTrialOrder = C{1};
clear C randomTrialOrderFile
fclose(fid);
randomHeader = randomTrialOrder(1,:);
randomTrialOrder = f_sort_trials(randomTrialOrder);

% set results folder (folder with .m files from 'read_c3d_files.m')
resultsPath = uigetdir('','choose the results folder, e.g. C:\...\read_c3d_davos\results');
outputFolder = uigetdir('','choose the output folder, e.g. C:\...\statistics');

% set relevant szenarios for this study
szenarios = {'practice','savings','transfer'};


for s = 1:length(szenarios)
    num_trials = find(blockedTrialOrder(:,s)==0);
    if isempty(num_trials)
        num_trials = length(blockedTrialOrder);
    else
        num_trials = num_trials(1)-1;
    end
    
    szenario = szenarios{s};
    trialDuration = zeros(length(subjDistribution)-1,num_trials);
    trialPD = zeros(length(subjDistribution)-1,num_trials,550);
    FFC_factor = zeros(length(subjDistribution)-1,num_trials);
    trialPD_raw = zeros(length(subjDistribution)-1,num_trials,1199);
    trialPD_change = zeros(length(subjDistribution)-1,num_trials,549);
    trialPD_change_raw = zeros(length(subjDistribution)-1,num_trials,1200);
    trialVEL = zeros(length(subjDistribution)-1,num_trials,551);
    trialACC = zeros(length(subjDistribution)-1,num_trials,550);
    statistics_header = cell(length(subjDistribution)-1,5);
    for i = 2:length(subjDistribution) %i=1 is the header of the table
        targetOrder = cell2mat(subjDistribution(i,4));
        condition = cell2mat(subjDistribution(i,2));
        load([resultsPath '\' 'output_x_x_' cell2mat(subjDistribution(i,5)) '_' szenarios{s} '.mat']);
    
        % different trial orders and rearrangements for each szenario
        data = cell(1,length(folder_struct_selected));
        
        if strcmp(szenarios{s},'practice')
            dataset = eval([condition 'TrialOrder']);
            vector_trials = dataset(:,1+(str2double(targetOrder(2))-1)*3);
            sOrder =  f_sort_trials_again(folder_struct_selected,vector_trials);
            data(1,sOrder) = folder_struct_selected;
            folder_struct_selected = data;
            clear dataset vector_trials sOrder data
            cellLength = 155;
        elseif strcmp(szenarios{s},'savings')
            dataset = eval([condition 'TrialOrder']);
            stop_idx = find(dataset(:,2)==0);
            vector_trials = dataset(1:stop_idx(1)-1,2+(str2double(targetOrder(2))-1)*3);
            sOrder =  f_sort_trials_again(folder_struct_selected,vector_trials);
            data(1,sOrder) = folder_struct_selected;
            folder_struct_selected = data;  
            clear dataset vector_trials sOrder data
            cellLength = 51;
        elseif strcmp(szenarios{s},'transfer')
            dataset = eval([condition 'TrialOrder']);
            if strcmp(cell2mat(subjDistribution(i,5)),'DAVOS1078')
                stop_idx = 15;
            elseif strcmp(cell2mat(subjDistribution(i,5)),'DAVOS1080')
                stop_idx = 28;
            else
                stop_idx = find(dataset(:,2)==0);
            end 
            vector_trials = dataset(1:stop_idx(1)-1,3+(str2double(targetOrder(2))-1)*3);
            sOrder =  f_sort_trials_again(folder_struct_selected,vector_trials);
            data(1,sOrder) = folder_struct_selected;
            folder_struct_selected = data; 
            clear dataset vector_trials sOrder data
            cellLength = 51;
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % from here parameter or statistics output
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % create header with all subject, group, target information
        statistics_header{i-1,1} = cell2mat(subjDistribution(i,1));
        statistics_header{i-1,2} = cell2mat(subjDistribution(i,2));
        statistics_header{i-1,3} = cell2mat(subjDistribution(i,3));
        statistics_header{i-1,4} = cell2mat(subjDistribution(i,4));
        statistics_header{i-1,5} = cell2mat(subjDistribution(i,5));
        
        %% first parameter: trial duration
        for trial = 1:length(folder_struct_selected)
            trialDuration(i-1,trial) = folder_struct_selected{trial}.trial_stop_seconds - folder_struct_selected{trial}.trial_start_seconds;
        end
        
        %% second parameter: perpendicular displacement
        meanTrlDuration = 550; % mean trial duration over all participants. If you change this --> change preallocation of all variables at the top!!!!
        
        % Butterworth filter for position and force data
        d_filter = designfilt('lowpassiir','DesignMethod','butter', ...
                'FilterOrder',2,'HalfPowerFrequency',12, 'SampleRate', 1000); 
        f_filter = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 24, 'SampleRate', 1000, 'DesignMethod', 'butter');
        
        for trial = 1:length(folder_struct_selected)
             
            x_path = folder_struct_selected{trial}.raw_pos_x_euclid_trial;
            x_path = filtfilt(d_filter,x_path);
            velocityX = diff(x_path)./s_time; 
            
            t_vector = length(x_path);
            n_vector = linspace(1,t_vector,meanTrlDuration);
            v_vector = 1:t_vector;  
            x_path_raw = x_path';
            
            x_path = interp1(v_vector,x_path,n_vector); % from here x_path is interpolated
                 
            % add start and end target position
            dataset = eval([condition 'TrialOrder']);
            trialIdx = num2str(dataset(trial,s+((str2double(targetOrder(2))-1)*3))); % find target number to add coordinates at the end of x-path
            trialIdx = str2double(trialIdx(2));
            x_path_raw = [x_pos(7), x_path_raw, x_pos(trialIdx)];
            x_path = [x_pos(7), x_path, x_pos(trialIdx)];
  
            y_path = folder_struct_selected{trial}.raw_pos_y_euclid_trial;
            y_path = filtfilt(d_filter,y_path);
            velocityY = diff(y_path)./s_time;
            
            y_path_raw = y_path';
            y_path = interp1(v_vector,y_path,n_vector);  % from here y_path is interpolated
            y_path_raw = [y_pos(7), y_path_raw, y_pos(trialIdx)];
            y_path = [y_pos(7), y_path, y_pos(trialIdx)];
            
            %compute PD
            y_v1 = y_path(end)-y_path(1);
            x_v1 = x_path(end)-x_path(1);
            vec1 = [x_v1;y_v1];
            PD_progress = zeros(1,length(x_path));
            for v = 1:length(x_path)
                x_v2 = x_path(v)-x_path(1);
                y_v2 = y_path(v)-y_path(1);
                vec2 = [x_v2;y_v2];
                vecAngle = atan2(det([vec1,vec2]),dot(vec1,vec2));
                PD_progress(v) = norm(vec2)*sin(vecAngle);
            end
            
            %compute raw PD (not interpolated)
            y_v1_raw = y_path_raw(end)-y_path_raw(1);
            x_v1_raw = x_path_raw(end)-x_path_raw(1);
            vec1_raw = [x_v1_raw;y_v1_raw];
            PD_progress_raw = zeros(1,length(x_path_raw));
            length_x_path = length(x_path_raw);
            for v = 1:length(x_path_raw)
                x_v2_raw = x_path_raw(v)-x_path_raw(1);
                y_v2_raw = y_path_raw(v)-y_path_raw(1);
                vec2_raw = [x_v2_raw;y_v2_raw];
                vecAngle_raw = atan2(det([vec1_raw,vec2_raw]),dot(vec1_raw,vec2_raw));
                PD_progress_raw(v) = norm(vec2_raw)*sin(vecAngle_raw);
            end
            
            if length(PD_progress_raw) > 1200 %trials longer than 1.2 seconds --> outliers
                PD_progress_raw = NaN(1,1200);
                length_x_path = 1200;
            else
                PD_progress_raw(end+1:1200) = NaN;
            end
            
            % subject x trial x sample cell array
            trialPD_raw(i-1,trial,:) = PD_progress_raw(2:end);
            trlPD_chng_smth = diff(PD_progress_raw(2:length_x_path-1));
            
            trlPD_chng_smth(end+1:1200) = NaN;
            trialPD_change_raw(i-1,trial,:) = trlPD_chng_smth;
            
            % subject x trial x sample cell array
            trialPD(i-1,trial,:) = PD_progress(2:end-1);
            trialPD_change(i-1,trial,:) = diff(PD_progress(2:end-1)); % filtfilt butterworth! instead of moving-average        

            %% computing error clamp force
            vec2_x_target = [-1, -1, 12, 1, 1, -12];
            x_path = folder_struct_selected{trial}.raw_pos_x_euclid_trial;
            x_path_raw = x_path';
                   
            % add start and end target position
            dataset = eval([condition 'TrialOrder']);
            trialIdx = num2str(dataset(trial,s+((str2double(targetOrder(2))-1)*3))); % find target number to add coordinates at the end of x-path
            trialIdx = str2double(trialIdx(2));
  
            y_path = folder_struct_selected{trial}.raw_pos_y_euclid_trial;
            y_path_raw = y_path';
            
            y_v1_raw = y_pos(trialIdx)-y_pos(7);
            x_v1_raw = x_pos(trialIdx)-x_pos(7);
            vec1_raw = [x_v1_raw;y_v1_raw]; % target direction vector
            
            start_sample = folder_struct_selected{trial}.trial_start_sample;
            stop_sample = folder_struct_selected{trial}.trial_stop_sample;
            x_forces = folder_struct_selected{trial}.Right_FS_ForceX(start_sample:stop_sample);
            y_forces = folder_struct_selected{trial}.Right_FS_ForceY(start_sample:stop_sample);
            
            x_forces = filtfilt(f_filter,x_forces);
            y_forces = filtfilt(f_filter,y_forces);
            
            length_x_path = length(x_path_raw);
            
            FFMatrix = [0 -1; 1 0]; %CCW
            vec2_raw = FFMatrix * vec1_raw; % orthogonal to target direction vector
            
            % find orthogonal vector of vec1_raw for each sample and
            % compute the projection of the force fector
            proj_force = [];
            VecAngle = [];
            for v = 1:length(x_path_raw)
                vecForce_raw = [x_forces(v);y_forces(v)];
                
                % force projection on vec2
                proj_force(v) = norm((dot(vecForce_raw,vec2_raw)/dot(vec2_raw,vec2_raw))*vec2_raw);
                VecAngle(v) = rad2deg(atan2(norm(cross([vecForce_raw' 0],[vec2_raw' 0])),dot(vecForce_raw, vec2_raw)));
            end
            VecAngle = VecAngle > 90;
            proj_force(VecAngle) = proj_force(VecAngle)*(-1); %sign force: if positive, subject is pushing into couter-clockwise direction
            
            % interpolate force (for plotting)
            t_vector = length(proj_force);
            n_vector = linspace(1,t_vector,100);
            v_vector = 1:t_vector;
            interp_proj_force = interp1(v_vector,proj_force,n_vector);
            interpForce(i-1,trial,:) = interp_proj_force;
            
            % compute raw force field compensation factor
            velocityX = filtfilt(f_filter,velocityX);
            velocityY = filtfilt(f_filter,velocityY);
            velocity = [velocityX'; velocityY'];
            
            % find the right FFMatrix according to subjects condition
            if strcmp(cell2mat(subjDistribution(i,2)),'blocked') && strcmp(cell2mat(subjDistribution(i,4)),'T1')
                FFMatrix = [0 20;-20 0];
            elseif strcmp(cell2mat(subjDistribution(i,2)),'blocked') && strcmp(cell2mat(subjDistribution(i,4)),'T2')
                FFMatrix = [0 10;-10 0];
            elseif strcmp(cell2mat(subjDistribution(i,2)),'blocked') && strcmp(cell2mat(subjDistribution(i,4)),'T3')
                FFMatrix = [0 15;-15 0];
            elseif strcmp(cell2mat(subjDistribution(i,2)),'blocked') && strcmp(cell2mat(subjDistribution(i,4)),'T4')
                FFMatrix = [0 15;-15 0];
            elseif strcmp(cell2mat(subjDistribution(i,2)),'blocked') && strcmp(cell2mat(subjDistribution(i,4)),'T5')
                FFMatrix = [0 20;-20 0];
            elseif strcmp(cell2mat(subjDistribution(i,2)),'blocked') && strcmp(cell2mat(subjDistribution(i,4)),'T6')
                FFMatrix = [0 10;-10 0];
            elseif strcmp(cell2mat(subjDistribution(i,2)),'random')
                FFMatrix = [0 15;-15 0];
            end
            
            % change FFMatrix after the first 6 EC-trials
            if strcmp(szenarios(s),'savings')
                if trial > 6
                    FFMatrix = [0 15;-15 0];
                end
            elseif strcmp(szenarios(s),'transfer')
                FFMatrix = [0 15;-15 0];
            end
            
            forceIdealArray = FFMatrix * velocity;

            % Compute force values for Data Points from force Array
            idealForce = [];
            for fct=1:length(forceIdealArray)
                idealForce(fct)= sqrt(forceIdealArray(1, fct)^2 + forceIdealArray(2, fct)^2);
            end

            % Linear regression fit force ideal and force measured
            p = polyfit(idealForce, proj_force(1:end-1), 1);

            FFC_factor(i-1,trial) = p(1)*100;   % The ForceFieldCompensationFactor
            
        end

        fprintf([szenario ': subject nr ' num2str(i-1) '\n']);
        
    end

    %% take the mean for each group of the PD 
    SB = [];
    SR = [];
    WB = [];
    WR = [];
    for sbj = 1:length(subjDistribution)-1
        if strcmp(cell2mat(subjDistribution(sbj+1,3)),'SB')
            SB = [SB; sbj];
        elseif  strcmp(cell2mat(subjDistribution(sbj+1,3)),'SR')
             SR = [SR; sbj];
        elseif  strcmp(cell2mat(subjDistribution(sbj+1,3)),'WB')
            WB = [WB; sbj];
        elseif  strcmp(cell2mat(subjDistribution(sbj+1,3)),'WR')
             WR = [WR; sbj];
        end
    end
    
    % be careful! there are still error clamps!!! practice: 145:150;
    % savings: 1:6, 13, 20, 27, 34, 41:46
    % transfer: 1:6, 13, 20, 27, 34, 41:46

%     trials_indices = [7:12, 14:19]; %first 6 trials
%     trials_indices = 14:19; %second 6 trials
%     trials_indices = [7:12, 14:19, 21:26, 28:33, 35:40]; % first 12 trials in savings and transfer
    
%     trials_indices = 31:36; %first 48 trials in practice
%     trials_indices = 139:144; %last 6 trials in practice
%     trials_indices = 133:144; %last 12 trials in practice

%     meanTrial_SB = squeeze(nanmean(nanmean(trialPD(SB,trials_indices,:),2),1));
%     meanTrial_SR = squeeze(nanmean(nanmean(trialPD(SR,trials_indices,:),2),1));
%     meanTrial_WB = squeeze(nanmean(nanmean(trialPD(WB,trials_indices,:),2),1));
%     meanTrial_WR = squeeze(nanmean(nanmean(trialPD(WR,trials_indices,:),2),1));
%     
%     semTrial_SB = squeeze(nanstd(nanmean(trialPD(SB,trials_indices,:),2),1))/sqrt(numel(SB));
%     semTrial_SR = squeeze(nanstd(nanmean(trialPD(SR,trials_indices,:),2),1))/sqrt(numel(SR));
%     semTrial_WB = squeeze(nanstd(nanmean(trialPD(WB,trials_indices,:),2),1))/sqrt(numel(WB));
%     semTrial_WR = squeeze(nanstd(nanmean(trialPD(WR,trials_indices,:),2),1))/sqrt(numel(WR));
%     
%     meanTrial_change_SB = squeeze(nanmean(nanmean(trialPD_change(SB,trials_indices,:),2),1));
%     meanTrial_change_SR = squeeze(nanmean(nanmean(trialPD_change(SR,trials_indices,:),2),1));
%     meanTrial_change_WB = squeeze(nanmean(nanmean(trialPD_change(WB,trials_indices,:),2),1));
%     meanTrial_change_WR = squeeze(nanmean(nanmean(trialPD_change(WR,trials_indices,:),2),1));
%     
%     semTrial_change_SB = squeeze(nanstd(nanmean(trialPD_change(SB,trials_indices,:),2),1))/sqrt(numel(SB));
%     semTrial_change_SR = squeeze(nanstd(nanmean(trialPD_change(SR,trials_indices,:),2),1))/sqrt(numel(SR));
%     semTrial_change_WB = squeeze(nanstd(nanmean(trialPD_change(WB,trials_indices,:),2),1))/sqrt(numel(WB));
%     semTrial_change_WR = squeeze(nanstd(nanmean(trialPD_change(WR,trials_indices,:),2),1))/sqrt(numel(WR));
%     
%     meanTrial_Vel_SB = squeeze(nanmean(nanmean(trialVEL(SB,trials_indices,:),2),1));
%     meanTrial_Vel_SR = squeeze(nanmean(nanmean(trialVEL(SR,trials_indices,:),2),1));
%     meanTrial_Vel_WB = squeeze(nanmean(nanmean(trialVEL(WB,trials_indices,:),2),1));
%     meanTrial_Vel_WR = squeeze(nanmean(nanmean(trialVEL(WR,trials_indices,:),2),1));
%     
%     meanTrial_Acc_SB = squeeze(nanmean(nanmean(trialACC(SB,trials_indices,:),2),1));
%     meanTrial_Acc_SR = squeeze(nanmean(nanmean(trialACC(SR,trials_indices,:),2),1));
%     meanTrial_Acc_WB = squeeze(nanmean(nanmean(trialACC(WB,trials_indices,:),2),1));
%     meanTrial_Acc_WR = squeeze(nanmean(nanmean(trialACC(WR,trials_indices,:),2),1));
%     
%     plot_pd
        
    %% compute statistic tables for export
      
    above_thrsh = trialPD>0.0045;
    below_thrsh = trialPD<-0.0045;
    first_binarized_error = bsxfun(@plus,above_thrsh,below_thrsh);
    
    if s == 1
        ff_trls = [1:144];
        errClamp_trls = [145:150];
    else
        ff_trls = [7:12, 14:19, 21:26, 28:33, 35:40];
        errClamp_trls = 1:46;
        errClamp_trls(ff_trls) = [];
    end
    
    stat_FFCF = FFC_factor(:,errClamp_trls);
    
    stat_area_curve = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_pred_PD120 = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_pred_PD150 = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_pred_PD180 = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_max_PD = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_mean_PD = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_binary_error = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_direction_change = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_first_reaction = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_direction_change_raw = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_first_reaction_raw = zeros(length(subjDistribution)-1,length(ff_trls));
    stat_binary_time_elapse = zeros(length(subjDistribution)-1,length(ff_trls));
    for sbj = 1:size(trialPD,1)
        for trl = 1:size(trialPD,2)
%             stat_pred_PD120(sbj,trl) = max(trialPD(sbj,trl,2:120),[],3); % prediction in the first 120 ms
%             stat_pred_PD150(sbj,trl) = max(trialPD(sbj,trl,2:150),[],3); % prediction in the first 150 ms
%             stat_pred_PD180(sbj,trl) = max(trialPD(sbj,trl,2:180),[],3); % prediction in the first 180 ms
            stat_max_PD(sbj,trl) = max(abs(trialPD(sbj,trl,2:end-1)),[],3); % maximum PD over the whole trial
            stat_mean_PD(sbj,trl) = squeeze(mean(trialPD(sbj,trl,2:end-1),3)); % mean PD over the whole trial
            rectify_area_curve = squeeze(abs(trialPD(sbj,trl,2:end-1)));
            idx_channel = rectify_area_curve <= 0.0045;
            rectify_area_curve(idx_channel) = 0;
            stat_area_curve(sbj,trl) = sum(rectify_area_curve)/meanTrlDuration; % area under the curve
            
            % find "0" in trial --> where gets the PD negative or positive
            zeroIDX = [];
            for i = 150:300
                if (trialPD(sbj,trl,i)>0 && trialPD(sbj,trl,i+1)<0) || (trialPD(sbj,trl,i)<0 && trialPD(sbj,trl,i+1)>0) 
                    zeroIDX(i) = 1;
                end
            end
            zeroIDX = find(zeroIDX);
            
            % if "0" turnpoint was before 150 ms, take 150ms as idx
            % (exclude prediction)
            if isempty(zeroIDX)
                neg_idx = 150;
            else
                neg_idx = zeroIDX(1);
            end
            neg_idx = neg_idx+10;
            copy_neg_idx = neg_idx;
            
            % compute for interpolated PD
            if ~isempty(isnan(squeeze(trialPD(sbj,trl,:)))==0) %if usual trial and not outlier (outliers are NaN)
                
                % find "0" turnpoint for the derivation (velocity)
                zeroIDX = [];
                for i = neg_idx:length(trialPD_change(sbj,trl,:))-1
                    if (trialPD_change(sbj,trl,i)>0 && trialPD_change(sbj,trl,i+1)<0) || (trialPD_change(sbj,trl,i)<0 && trialPD_change(sbj,trl,i+1)>0) 
                        zeroIDX(i) = 1;
                    end
                end
                zeroIDX = find(zeroIDX);
            
                if isempty(zeroIDX)
                    stat_direction_change(sbj,trl) = NaN;
                    stat_first_reaction(sbj,trl) = NaN;
                else
                    stat_direction_change(sbj,trl) = zeroIDX(1); % time point at velocity change (peak of PD)

                    if (zeroIDX(1) - neg_idx) < 10 %in case that zeroIDX and neg_idx are to narrow
                        neg_idx = 100;
                    end
                    [~,min_idx] = min(abs(diff(trialPD_change(sbj,trl,neg_idx+1:zeroIDX(1)-1))),[],3);
                    stat_first_reaction(sbj,trl) = min_idx(end) + neg_idx+1; % time point at acceleration change near to velocity change
                end
            else
                stat_direction_change(sbj,trl) = NaN;
                stat_first_reaction(sbj,trl) = NaN;
            end
            
            neg_idx = copy_neg_idx;
            
            % compute for non-interpolated (raw) PD
            if ~isempty(isnan(squeeze(trialPD_raw(sbj,trl,1:200))) == 0) %if usual trial and not outlier (outliers are NaN)
                stat_pred_PD120(sbj,trl) = max(trialPD_raw(sbj,trl,2:120),[],3); % prediction in the first 120 ms
                stat_pred_PD150(sbj,trl) = max(trialPD_raw(sbj,trl,2:150),[],3); % prediction in the first 150 ms
                stat_pred_PD180(sbj,trl) = max(trialPD_raw(sbj,trl,2:180),[],3); % prediction in the first 180 ms
                stat_max_PD200(sbj,trl) = max(abs(trialPD_raw(sbj,trl,200:end-1)),[],3); % maximum PD over the whole trial

                % find "0" turnpoint for the derivation (velocity)
                zeroIDX = [];
                for i = neg_idx:length(trialPD_change_raw(sbj,trl,:))-1
                    if (trialPD_change_raw(sbj,trl,i)>0 && trialPD_change_raw(sbj,trl,i+1)<0) || (trialPD_change_raw(sbj,trl,i)<0 && trialPD_change_raw(sbj,trl,i+1)>0) 
                        zeroIDX(i) = 1;
                    end
                end
                zeroIDX = find(zeroIDX);
            
                if isempty(zeroIDX)
                    stat_direction_change_raw(sbj,trl) = NaN;
                    stat_first_reaction_raw(sbj,trl) = NaN;
                else
                    stat_direction_change_raw(sbj,trl) = zeroIDX(1); % time point at velocity change (peak of PD)

                    if (zeroIDX(1) - neg_idx) < 10 %in case that zeroIDX and neg_idx are to narrow
                        neg_idx = 100;
                    end
                    [~,min_idx] = min(abs(diff(trialPD_change_raw(sbj,trl,neg_idx+1:zeroIDX(1)-1))),[],3);
                    stat_first_reaction_raw(sbj,trl) = min_idx(end) + neg_idx+1; % time point at acceleration change near to velocity change
                end
            else
                stat_direction_change_raw(sbj,trl) = NaN;
                stat_first_reaction_raw(sbj,trl) = NaN;
            end
            
            % binarized data (find PD bigger than target radius)
            idx_first_binary = find(first_binarized_error(sbj,trl,:));
            if isempty(idx_first_binary)
                stat_binary_error(sbj,trl) = NaN;
                stat_binary_time_elapse(sbj,trl) = NaN;
            else
                stat_binary_error(sbj,trl) = idx_first_binary(1); % first error bigger than target radius
                stat_binary_time_elapse(sbj,trl) = sum(idx_first_binary); % number of time bins bigger than target radius
            end
        end
    end

    % binary statistics
    idx_binary = isnan(stat_binary_error);
    
%     stat_binary_predPD = stat_pred_PD;
%     stat_binary_predPD(idx_binary) = 0;
    stat_binary_pdmax = stat_max_PD;
    stat_binary_pdmax(idx_binary) = 0; % max PD only for trials which have a binarized error
    stat_binary_pdmean = stat_mean_PD;
    stat_binary_pdmean(idx_binary) = 0;
    stat_binary_direction_change = stat_direction_change;
    stat_binary_direction_change(idx_binary) = NaN;
    stat_binary_first_reaction = stat_first_reaction;
    stat_binary_first_reaction(idx_binary) = NaN;
    
    %remove error clamp trials
    parameters = {'trialDuration','stat_pred_PD120','stat_pred_PD150','stat_pred_PD180','stat_max_PD','stat_max_PD200','stat_mean_PD','stat_direction_change','stat_first_reaction','stat_binary_pdmax', ...
        'stat_binary_pdmean','stat_binary_direction_change','stat_binary_first_reaction','stat_binary_time_elapse','stat_area_curve','stat_direction_change_raw','stat_first_reaction_raw'};
    for stats = 1:length(parameters)
        var = eval(parameters{stats});
        var(:,errClamp_trls) = [];
        assignin('base',parameters{stats},var);
    end

    trialPD(:,errClamp_trls,:) = [];
    
    %% save .mat and .csv for statistics
   
    
    
    % write csv
    parameters(end+1) = {'stat_FFCF'};
    parameters(end+1) = {'statistics_header'};
    parameters(end+1) = {'trialPD'};
    
    for stats = 1:length(parameters)
        save([outputFolder '\' parameters{stats} '_' szenario '.mat'],parameters{stats}, '-v7.3');
        
        fid = fopen([outputFolder '\' 'Results_' parameters{stats} '_' szenario '.csv'],'w');
        for z = 1:size(eval(parameters{stats}),1)
            for si = 1:size(eval(parameters{stats}),2)            
                var = eval(parameters{stats});
                var = var(z,si);

                if size(var,1) == 0
                    var = '';
                end

                if isnumeric(var) == 1
                    var = num2str(var);
                end

                if iscell(var) ==1
                    var = var{1};
                end
                
                fprintf(fid,var);

                if si ~= size(trialDuration,2)
                    fprintf(fid,[';']);
                end
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        fprintf(['saved ' parameters{stats} '\n']);
    end
    
    
end



    
