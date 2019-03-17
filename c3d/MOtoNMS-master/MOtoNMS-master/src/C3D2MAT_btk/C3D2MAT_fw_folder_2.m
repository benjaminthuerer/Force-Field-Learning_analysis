%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               MOtoNMS                                   %
%                MATLAB MOTION DATA ELABORATION TOOLBOX                   %
%                 FOR NEUROMUSCULOSKELETAL APPLICATIONS                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read a c3d file and store main information in mat files
%
% USAGE: Input files folder path MUST included a folder named 'InputData'
%       Output folders may be managed in mkOutputPath.m
%       Cancel the created sessionData folder if you modify the code before
%       rerunning it: data already saved are not overwritten!

% The file is part of matlab MOtion data elaboration TOolbox for
% NeuroMusculoSkeletal applications (MOtoNMS). 
% Copyright (C) 2012-2014 Alice Mantoan, Monica Reggiani
%
% MOtoNMS is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% Matlab MOtion data elaboration TOolbox for NeuroMusculoSkeletal applications
% is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along 
% with MOtoNMS.  If not, see <http://www.gnu.org/licenses/>.
%
% Alice Mantoan, Monica Reggiani
% <ali.mantoan@gmail.com>, <monica.reggiani@gmail.com>

%%

function folder_content = C3D2MAT_fw_folder2(folderPath)

addSharedPath()
runTerminalNote()

%% Selection of input data 
pathName = folderPath;%uigetdir('Select your input data folder');
%c3dFiles = dir ([pathName filesep '*.c3d']);
c3dFiles = rdir([pathName filesep '\**\*.c3d']);

w = waitbar(0,'Elaborating data...Please wait!');

folder_content = {};
for k=1:length(c3dFiles)
    
    %correction of the name --> after uniformation it should not be necessary
    %trialsName{k} = regexprep(regexprep((regexprep(c3dFiles(k).name, ' ' , '')), '-',''), '.c3d', '');
    
    %Folders and paths creation
    c3dFilePathAndName = c3dFiles(k).name;
    [fpathstr,fname,fext] = fileparts(c3dFilePathAndName);
       
    if k == 1
       folder_content.all = {};
       folder_content.analogs = {};
       folder_content.Markers = {};
       folder_content.c3dFilePathAndName = {};
       folder_content.c3dFileName = {};
       folder_content.c3dFileExt = {};
       folder_content.c3dFileFolderPath = {};


    end
    
    %Data Reading    
    [all, analogs, Markers, AnalogData, FPdata, Events, ForcePlatformInfo, Rates] = getInfoFromC3D(c3dFilePathAndName);
    folder_content.all{k} = all;
    folder_content.analogs{k} = analogs;
    folder_content.Markers{k} = Markers;
    folder_content.c3dFilePathAndName{k} = c3dFilePathAndName;
    folder_content.c3dFileName{k} = [fname];
    folder_content.c3dFileExt{k} = [fext];
    folder_content.c3dFileFolderPath{k} = [fpathstr];
    
    
    
%     trialMatFolder=mkOutputPath(pathName,trialsName{k});
%     
%     sessionFolder=regexprep(trialMatFolder, [trialsName{k} filesep], '');
%     
%     %Consistency check and Storing: 
%     %only if the trial is not a static because it may have different data
%     %(in that case their are saved anyway in the static trial folder)
%     
%     if isempty(strfind(upper(trialsName{k}),'STATIC'))
%         %Common Session Info (excluding static trials)
%         if isempty(Markers)
%             dMLabels=[];
%         else
%             dMLabels=Markers.Labels;
%         end
%         if isempty(AnalogData)
%             AnalogDataLabels=[];
%         else
%             AnalogDataLabels=AnalogData.Labels;
%         end
%         
%         checkAndSaveSessionInfo(ForcePlatformInfo, Rates, dMLabels, AnalogDataLabels, sessionFolder);
%     end
%     %Data for each trials

    waitbar(k/length(c3dFiles));    
end
close(w)

%Saving trialsName list of read c3d file at the end
%save([sessionFolder 'trialsName.mat'],'trialsName')

%save_to_base(1)
% save_to_base() copies all variables in the calling function to the base
% workspace. This makes it possible to examine this function internal
% variables from the Matlab command prompt after the calling function
% terminates. Uncomment the following command if you want to activate it
     