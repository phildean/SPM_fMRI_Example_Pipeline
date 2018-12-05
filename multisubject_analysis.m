function [mem] = multisubject_analysis(shift)
% Updated Dec 2018
% Based on script by Adam McNamara, edited by Philip Dean

% ALL ANALYSIS SCRIPTS RUN FROM HERE
% Runs all the specific programs for each subject in turn:
%   preprocess
%   firstlevel_analysis
%       including behavioural_data_onset_duration
%   secondlevel_analysis
%
% Can add any other analysis into here if wanted. 

% subject independent information such as study paths
% the scripts assume this setup for your data. 
sh.studypath = 'E:\MRI\BECi_Study\Data';
sh.imagepath = 'E:\MRI\BECi_Study\Raw_Data';
sh.behavpath = 'E:\MRI\BECi_Study\Behavioural_Data';
% subj.path is the folder name for your participant
% subj.log_sess1 is the file name for your behavioural logfile
% subj.response_button is the response button used for the experiment (if counterbalanced)

mypwd = pwd;
addpath(mypwd);
global ev
ev = {};

%% subject-specific information such as logfile names and response buttons

%for jj = 1:11; subj(jj).path = ['Subject_' num2str(jj)]; end;  % this is one way of specifying subject path or could do the below for more specific/random subject names

subj(1).path = 'Subject_1'; subj(1).log_sess1='P01_S1_Targ1_BECi_nBack_SCE.log'; subj(1).response_button = 1;       % subj(1).log_sess2='P01_S2_Targ1_BECi_nBack_SCE_REDONE.log'; subj(1).log_sess3='PO1_S3_Targ1_BECi_nBack_SCE.log'; 
subj(2).path = 'Subject_2'; subj(2).log_sess1='P02_S1_Targ2_BECi_nBack_SCE.log'; subj(2).response_button = 2;       % subj(2).log_sess2='P02_S2_Targ2_BECi_nBack_SCE_REDONE.log'; subj(2).log_sess3='PO2_S3_Targ2_BECi_nBack_SCE.log'; 
subj(3).path = 'Subject_3'; subj(3).log_sess1='P03_S1_Targ2_BECi_nBack_SCE.log'; subj(3).response_button = 2;       % subj(3).log_sess2='P03_S2_Targ2_BECi_nBack_SCE.log'; subj(3).log_sess3='P03_S3_Targ2_BECi_nBack_SCE.log'; 
subj(4).path = 'Subject_4'; subj(4).log_sess1='P04_S1_Targ1_BECi_nBack_SCE.log'; subj(4).response_button = 1;       % subj(4).log_sess2='P04_S2_Targ1_BECi_nBack_SCE.log'; subj(4).log_sess3='P04_S3_Targ1_BECi_nBack_SCE.log';
subj(5).path = 'Subject_5'; subj(5).log_sess1='P05_S1_Targ1_BECi_nBack_SCE.log'; subj(5).response_button = 1;       % subj(5).log_sess2='P05_S2_Targ1_BECi_nBack_SCE.log'; subj(5).log_sess3='P05_S3_Targ1_BECi_nBack_SCE.log';
subj(6).path = 'Subject_6'; subj(6).log_sess1='P06_S1_Targ2_BECi_nBack_SCE.log'; subj(6).response_button = 2;       % subj(6).log_sess2='P06_S2_Targ2_BECi_nBack_SCE.log'; subj(6).log_sess3='P06_S3_Targ2_BECi_nBack_SCE.log'; 
subj(7).path = 'Subject_7'; subj(7).log_sess1='P07_S1_Targ2_BECi_nBack_SCE.log'; subj(7).response_button = 2;       % subj(7).log_sess2='P07_S2_Targ2_BECi_nBack_SCE.log'; subj(7).log_sess3='P07_S3_Targ2_BECi_nBack_SCE.log'; 
subj(8).path = 'Subject_8'; subj(8).log_sess1='P08_S1_Targ1_BECi_nBack_SCE.log'; subj(8).response_button = 1;       % subj(8).log_sess2='P08_S2_Targ1_BECi_nBack_SCE.log'; subj(8).log_sess3='P08_S3_Targ1_BECi_nBack_SCE.log'; 
subj(9).path = 'Subject_9'; subj(9).log_sess1='P09_S1_Targ2_BECi_nBack_SCE.log'; subj(9).response_button = 2;       % subj(9).log_sess2='P09_S2_Targ2_BECi_nBack_SCE.log'; subj(9).log_sess3='P09_S3_Targ2_BECi_nBack_SCE.log'; 
subj(10).path = 'Subject_10'; subj(10).log_sess1='P10_S1_Targ1_BECi_nBack_SCE.log'; subj(10).response_button = 1;   % subj(10).log_sess2='P10_S2_Targ1_BECi_nBack_SCE.log'; subj(10).log_sess3='P10_S3_Targ1_BECi_nBack_SCE.log'; 
subj(11).path = 'Subject_11'; subj(11).log_sess1='P11_S1_Targ1_BECi_nBack_SCE.log'; subj(11).response_button = 1;   % subj(11).log_sess2='P11_S2_Targ1_BECi_nBack_SCE.log'; subj(11).log_sess3='P11_S3_Targ1_BECi_nBack_SCE.log'; 


for i = 1:length(subj)

spm('Defaults', 'FMRI');        % Reset SPM defaults for fMRI (not sure necessary - safety catch?)

%first define directory name
   fprintf(subj(i).path);
   fprintf('\n');
  
%% Stage 1: do preprocessing ('drsbcnog')
    fprintf(' --- Preprocessing fMRI --- \n');
    preprocess('drsbcnog', fullfile(sh.studypath,subj(i).path));
    
%% Stage 2: do first level analysis ('me')
    %fprintf('  --- First level analysis: GLM --- \n');
    %firstlevel_analysis('mec', fullfile(sh.studypath,subj(i).path), sh, subj(i)); 
end;

%% Stage 3: do second level analysis
    %fprintf(' --- Second level analysis: Paired-Test GLM & Parametric --- \n');
    %secondlevel_analysis('o',sh.studypath);