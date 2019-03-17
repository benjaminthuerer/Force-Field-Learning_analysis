clear
dbstop if error
%mainStudyPath = 'C:\Users\cleaner\Desktop\lab\Thürer\study\DAVOS';
mainStudyPath = uigetdir('','choose the study folder, e.g. C:\...\DAVOS');

scriptPath = [mainStudyPath '\scripts\matlab'];
toolbox_path = [scriptPath '\c3d\MOtoNMS-master\MOtoNMS-master\src\C3D2MAT_btk'];


%folderPathWithC3Dfiles = ['C:\Users\cleaner\Desktop\lab\Thürer\study\DAVOS\behavioral_data\c3d_data\test_c3d'];
folderPathWithC3Dfiles = uigetdir(mainStudyPath,'choose the input folder, e.g. C:\...\StudyIII_Sleep');
%pathOutputFolder = ['C:\Users\cleaner\Desktop\lab\Thürer\study\DAVOS\behavioral_data\c3d_data\test_c3d_output'];
pathOutputFolder = uigetdir(mainStudyPath,'choose the output folder, e.g. C:\...\results');

if ~isdir(pathOutputFolder)
    mkdir(pathOutputFolder);
end


%% run statistics separated for every c3d-file
% first check folder for all subject names
list = dir(folderPathWithC3Dfiles);
logicList = cell2mat({list.isdir});
intList = find(logicList);
intList(1:2) = [];

subjNames = {list.name};

% create a loop which does statistics for every subject
for i = 1:length(intList)
    subjFolder = [folderPathWithC3Dfiles '\' cell2mat(subjNames(intList(i))) ];
    subjNameToSave = cell2mat(subjNames(intList(i)));
    
    % check subjec folder for all datasets
    subj_list = dir(subjFolder);
    subj_logicList = cell2mat({subj_list.isdir});
    subj_intList = find(subj_logicList);
    subj_intList(1:2) = [];
    
    fileNames = {subj_list.name};
    for ii = 1:length(subj_intList) %1=training; 2=savings; 3=transfer
        fileFolder = [subjFolder '\' cell2mat(fileNames(subj_intList(ii))) ];
        if ii == 1
            szenarioToSave = 'practice';
        elseif ii == 2
            szenarioToSave = 'savings';
        elseif ii == 3
            szenarioToSave = 'transfer';
        else
            error(['more than 3 szenarios in subject folder ' subjNameToSave]);
        end

        % start c3d read (scripted by Freddy)
        
        % documentation on how to use at http://rehabenggroup.github.io/MOtoNMS/manual/C3DtoMAT.html
        addpath(scriptPath)
        addpath([scriptPath '\c3d\btk'])
        addpath(toolbox_path)
        cd(toolbox_path)

        all_folder_structs = {};

        folder_content = C3D2MAT_fw_folder_2(fileFolder);
        %saveMat(folder_content)



        %all_folder_structs.folderPathWithC3Dfiles
        folder_struct_selected = {};
        for iC3Dfile = 1:size(folder_content.all,2)

        temp_struct = {};
        temp_struct.file_name = folder_content.c3dFileName{1, iC3Dfile};
        temp_struct.file_path_and_name = folder_content.c3dFilePathAndName{1, iC3Dfile};

        ['working on file: ' temp_struct.file_path_and_name]

        temp_trial_chars = strsplit(temp_struct.file_name,'_');

        temp_struct.subject = str2num(temp_trial_chars{1});
        temp_struct.trial = str2num(temp_trial_chars{2});
        temp_struct.repeat = str2num(temp_trial_chars{3});



        temp_struct.fsample = folder_content.Markers{1, iC3Dfile}.Rate;
        temp_struct.unit = folder_content.Markers{1, iC3Dfile}.Units;

        temp_struct.raw_pos_x_euclid = folder_content.Markers{1, iC3Dfile}.RawData(:,1);
        temp_struct.raw_pos_y_euclid = folder_content.Markers{1, iC3Dfile}.RawData(:,2);

        temp_struct.marker_labels = folder_content.all{1, iC3Dfile}.children.EVENTS.children.LABELS.info.values;
        temp_struct.marker_seconds = folder_content.all{1, iC3Dfile}.children.EVENTS.children.TIMES.info.values;
        
          %Adjust the marker!
        temp_struct.trial_start_seconds = temp_struct.marker_seconds(2);
        temp_struct.trial_stop_seconds = temp_struct.marker_seconds(3);
%         temp_struct.trial_start_seconds = temp_struct.marker_seconds(1);
%         temp_struct.trial_stop_seconds = temp_struct.marker_seconds(2);
        temp_struct.trial_start_sample = max(1,round(temp_struct.trial_start_seconds*temp_struct.fsample));
        temp_struct.trial_stop_sample = min(round(temp_struct.trial_stop_seconds*temp_struct.fsample),numel(temp_struct.raw_pos_x_euclid));

        temp_struct.raw_pos_x_euclid_trial = temp_struct.raw_pos_x_euclid(temp_struct.trial_start_sample:temp_struct.trial_stop_sample);
        temp_struct.raw_pos_y_euclid_trial = temp_struct.raw_pos_y_euclid(temp_struct.trial_start_sample:temp_struct.trial_stop_sample);

        temp_struct.Right_L1Ang = folder_content.analogs{1, iC3Dfile}.Right_L1Ang_______;
        temp_struct.Right_L2Ang = folder_content.analogs{1, iC3Dfile}.Right_L2Ang_______;
        temp_struct.Right_L1Vel = folder_content.analogs{1, iC3Dfile}.Right_L1Vel_______;
        temp_struct.Right_L2Vel = folder_content.analogs{1, iC3Dfile}.Right_L2Vel_______;
        temp_struct.Right_L1Acc = folder_content.analogs{1, iC3Dfile}.Right_L1Acc_______;
        temp_struct.Right_L2Acc = folder_content.analogs{1, iC3Dfile}.Right_L2Acc_______;
        temp_struct.FS_TimeStamp = folder_content.analogs{1, iC3Dfile}.Right_FS_TimeStamp;
        temp_struct.Right_M1TorCMD = folder_content.analogs{1, iC3Dfile}.Right_M1TorCMD____;
        temp_struct.Right_M2TorCMD = folder_content.analogs{1, iC3Dfile}.Right_M2TorCMD____;
        temp_struct.Right_FS_ForceX = folder_content.analogs{1, iC3Dfile}.Right_FS_ForceX___;
        temp_struct.Right_FS_ForceY = folder_content.analogs{1, iC3Dfile}.Right_FS_ForceY___;
        temp_struct.Right_FS_ForceZ = folder_content.analogs{1, iC3Dfile}.Right_FS_ForceZ___;
        temp_struct.Right_FS_TorqueX = folder_content.analogs{1, iC3Dfile}.Right_FS_TorqueX__;
        temp_struct.Right_FS_TorqueY = folder_content.analogs{1, iC3Dfile}.Right_FS_TorqueY__;
        temp_struct.Right_FS_TorqueZ = folder_content.analogs{1, iC3Dfile}.Right_FS_TorqueZ__;


        folder_struct_selected{iC3Dfile} = temp_struct(:);

        end
        clear folder_content

        save([pathOutputFolder filesep 'output_' subjNameToSave '_' szenarioToSave '.mat'],'folder_struct_selected', '-v7.3');

        cd(pathOutputFolder)
    end
end

%% loading section

% clear
% 
% load(uigetfile('','choose the the file to import folder, e.g. C:\...\output_save.mat'));
% 
% figure;
% hold on
% for iTrial = 1:numel(folder_struct_selected)
%     duration(iTrial) = folder_struct_selected{iTrial}.trial_stop_seconds - folder_struct_selected{iTrial}.trial_start_seconds;
%     x = folder_struct_selected{iTrial}.raw_pos_x_euclid_trial';
%     y = folder_struct_selected{iTrial}.raw_pos_y_euclid_trial';
%     
%     plot(x,y,'LineWidth',0.5,...
%         'MarkerEdgeColor','k')
%     
% end
% hold off
% 
% 
% min_vel = 0;
% max_vel = 0;
% vel = {};
% for iTrial = 1:numel(folder_struct_selected)
%     vel_x = folder_struct_selected{iTrial}.Right_L1Vel(folder_struct_selected{iTrial}.trial_start_sample:folder_struct_selected{iTrial}.trial_stop_sample);
%     vel_y = folder_struct_selected{iTrial}.Right_L2Vel(folder_struct_selected{iTrial}.trial_start_sample:folder_struct_selected{iTrial}.trial_stop_sample);
% 
%     vel{iTrial} = sqrt(vel_x.^2 + vel_x.^2);
%     
%     min_vel = min(min_vel, min(vel{iTrial}));
%     max_vel = max(max_vel, max(vel{iTrial}));
% end
% 
% figure;
% hold on
% for iTrial = 1:numel(folder_struct_selected)
% 
%     x = folder_struct_selected{iTrial}.raw_pos_x_euclid_trial';
%     y = folder_struct_selected{iTrial}.raw_pos_y_euclid_trial';
%     z = zeros(size(x));
%     
%     %col = linspace(0,1,numel(x));  % This is the color
%     col = linspace(min_vel,max_vel,numel(x));  % This is the color
%     surface([x;x],[y;y],[z;z],[col;col], ...
%         'facecol','no', ...
%         'edgecol','interp', ...
%         'linew',2);
%     
%     plot(x,y,'LineWidth',0.5,...
%         'MarkerEdgeColor','k')
%     
% end
% hold off
% 
% colormap(jet(256))
% colorbar


%% old 
% 
% baseDataPath = [mainStudyPath 'behavioral_data\c3d_data'];
% outDataPath = [baseDataPath '\ElaboratedData\DAVOS1001\sessionData\';
%     
% currentTrialName = '01_01_01';
% load([outDataPath currentTrialName '\all.mat'])
% load([outDataPath currentTrialName '\analogs.mat'])
% load([outDataPath currentTrialName '\Markers.mat'])
% 
% trial_markers_seconds = all.children.EVENTS.children.TIMES.info.values;
% 
% figure;
% plot((1:numel(Markers.RawData(:,1)))/2000,Markers.RawData(:,1))
% figure;
% plot((1:numel(Markers.RawData(:,2)))/2000,Markers.RawData(:,2))
% 
% figure;
% plot(Markers.RawData(:,1),Markers.RawData(:,2))
% 
% figure;
% plot(analogs.Right_L1Ang_______,analogs.Right_L2Ang_______)
% 
% figure;
% plot(analogs.Right_L1Acc_______,analogs.Right_L2Acc_______)


