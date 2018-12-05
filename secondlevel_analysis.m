function secondlevel_analysis(todo,D)
% Dec 2018 redone for spm12
% Based on script by Adam McNamara, edited by Philip Dean

% MAKE AND ESTIMATE SECOND LEVEL MODEL
% Make model using first-level contrast images
% Estimate model to create SPM file for second level analysis
% Create second level contrasts

% INPUT ARGUMENTS:

% 'todo'
% o = make one-sample t-test
% t = make two-sample t-test                            %not done yet
% p = make paired sample t-test
% f = make Factorial model

% 'D'
% This is the Directory, e.g.  'E:\MRI\BECi_Study\Data'

% So could call script as:
% secondlevel_analysis('o','E:\MRI\BECi_Study\Data')


%% Global Variables
spm('Defaults', 'FMRI');        % Reset SPM defaults for fMRI (not sure necessary - safety catch?)
global defaults;                % Reset Global defaults (not sure why needed?)

% Groups to be analysed
group(1).subj=[1 2 3 4 5 6 7 8 9 10 11];group(1).name='All_11SJ';
%group(2).subj=[1 2 3 4 5 7 8 9 10 11];group(2).name='10SJ_No6';


% Analysis to be looked at
analysis_type={'GLM' 'Parametric' 'Bayesian'};

way='E:\MRI\BECi_Study\scripts\batch_files';    % Path to the "jobs"/batch files needed

% Contrast and Subject "precursors"
contrast_begin = 'con_00';          %this is added to with 01; 02.....10; 11 in script
subj_begin = 'Subject_';            % this is the folder name for your subjects and is added to with 01;..10 etc in script

% if folder for 2nd level within groups stats doesnt exist, create the folder
if ~exist(fullfile(D,'Second_Level_Stats\Within_Groups'),'dir'); mkdir(fullfile(D,'Second_Level_Stats\Within_Groups')); end;
withingroups_stats_directory = fullfile(D,'Second_Level_Stats\Within_Groups');

% if folder for 2nd level between groups stats doesnt exist, create the folder
%if ~exist(fullfile(D,'Second_Level_Stats\Between_Groups'),'dir'); mkdir(fullfile(D,'Second_Level_Stats\Between_Groups')); end;
%betweengroups_stats_directory = fullfile(D,'Second_Level_Stats\Between_Groups');

tic 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make model: 2nd level One-sample t-test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'o')
    
    if ~exist(fullfile(withingroups_stats_directory,'OneSampleTTest'),'dir'); mkdir(fullfile(withingroups_stats_directory,'OneSampleTTest')); end;
    onesampleTtest_directory = fullfile(withingroups_stats_directory,'OneSampleTTest');
    
    for aa = 1%:length(analysis_type)            %change it if want just GLM, just Parametric, Just Bayesian etc
    
        template_SPM = load(fullfile(D,'Subject_01\Stats', analysis_type{aa}, 'SPM.mat'));
        contrast_num = length(template_SPM.SPM.xCon);
        for nn = 1:length(template_SPM.SPM.xCon)
            contrast_name{nn} = template_SPM.SPM.xCon(nn).name;
            contrast_name{nn}(contrast_name{nn}==' ') = '_';
            contrast_name{nn}(contrast_name{nn}=='-') = '';
        end;  
                
        if ~exist(fullfile(onesampleTtest_directory,analysis_type{aa}),'dir'); 
            mkdir(fullfile(onesampleTtest_directory,analysis_type{aa}));    %Make GLM/Parametric etc directory within onesamplettest
        end
            
    for gg=1:length(group)   
        
        if ~exist(fullfile(onesampleTtest_directory,analysis_type{aa}, group(gg).name),'dir');
        mkdir(fullfile(onesampleTtest_directory,analysis_type{aa}, group(gg).name));    %Make directory for each set of subjects (e.g. 'All_11SJ')
        end;
        
    for cc=1:contrast_num
        
        if ~exist(fullfile(onesampleTtest_directory,analysis_type{aa}, group(gg).name, contrast_name{cc}),'dir');
        mkdir(fullfile(onesampleTtest_directory,analysis_type{aa}, group(gg).name, contrast_name{cc}));
        end;
        
        load(fullfile(way,'second_level_onesampleT_spm12.mat'));
        
        output_directory = fullfile(onesampleTtest_directory,analysis_type{aa}, group(gg).name, contrast_name{cc});
        
        matlabbatch{1}.spm.stats.factorial_design.dir = {output_directory};
    
    for ss=1:length(group(gg).subj)   
        
        contrast_end=['0' num2str(cc)]; if cc > 9; contrast_end=contrast_end(2:end); end;
        subj_end = ['0' num2str(group(gg).subj(ss))]; if length(subj_end) > 2; subj_end=subj_end(2:end); end;
   
        insert_contrast = [contrast_begin contrast_end '.nii'];
        insert_subj = [subj_begin subj_end];
        
        subject_directory = fullfile(D, insert_subj, 'Stats', analysis_type{aa});
        
        P=cellstr(spm_select('FPList', subject_directory, insert_contrast));
        
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{ss} = P{1};
        
        clear P;

    end; %ss (inputting each subjects contrast)
    
    %Other Possible Inputs
    %matlabbatch{1}.spm.stats.factorial_design.cov.c;                    %Covariate vector
    %matlabbatch{1}.spm.stats.factorial_design.cov.cname;                %Covariate name
    %matlabbatch{1}.spm.stats.factorial_design.cov.iCFi;                 %Covariate interactions (default = none)
    %matlabbatch{1}.spm.stats.factorial_design.cov.iCC;                  %Covariate centering (default = overall mean)
    %matlabbatch{1}.spm.stats.factorial_design.multi_cov.files;          %Multiple covariate files
    %matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCFI;           %Multiple covariate interactions
    %matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCC;            %Multiple covariate centering

    %matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;   %Threshold masking (default = none)
    %matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;           %Implicit mask (default = yes [1])
    %matlabbatch{1}.spm.stats.factorial_design.masking.em = {};          %Explicit mask 

    %matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;           %Global calculation (default = omit; used for PET)
    %matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;   %Global normalisation (default = no grand mean scaling; used for PET)
    %matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;          %Global normalisation (default = none; used for PET)
    
 
    %keyboard                  %Used in bug-testing script   

    spm_jobman('run',matlabbatch);
    
    clear matlabbatch;
    
    load(fullfile(way,'estimate_spm12.mat'));
    stats_directory=fullfile(onesampleTtest_directory, analysis_type{aa}, group(gg).name, contrast_name{cc});
    statsfile=fullfile(stats_directory,'SPM.mat');
    matlabbatch{1}.spm.stats.fmri_est.spmmat={statsfile};
    
    spm_jobman('run',matlabbatch);
    
    clear matlabbatch;
    
    load(fullfile(way,'contrast_manager_spm12.mat'));
    
    stats_directory=fullfile(onesampleTtest_directory, analysis_type{aa}, group(gg).name, contrast_name{cc});
    statsfile=fullfile(stats_directory,'SPM.mat');
    matlabbatch{1}.spm.stats.con.spmmat={statsfile};
    
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Main_Effect';         % t Contrast (f contrast = fcon)
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1];
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Main_Effect_Minus';         % t Contrast (f contrast = fcon)
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1];
    
    spm_jobman('run',matlabbatch);
    
    end;    %cc (inputting each contrast (e.g. con_0001.img))
    
    end;    %gg (inputting each group type (e.g. ALL_11SJ))
    
    end;    %aa (inputting each analysis type (e.g. GLM))

end;    %end 'o' one sample t-test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make model: 2nd level paired-sample t-test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'p')
    
    if ~exist(fullfile(withingroups_stats_directory,'PairedTTest'),'dir'); mkdir(fullfile(withingroups_stats_directory,'PairedTTest')); end;
    pairsampleTtest_directory = fullfile(withingroups_stats_directory,'PairedTTest');
    
    for aa = 1%:length(analysis_type)            %change it if want just GLM, just Parametric, Just Bayesian etc
    
        template_SPM = load(fullfile(D,'Subject_01\Stats', analysis_type{aa}, 'SPM.mat'));
        contrast_num = length(template_SPM.SPM.xCon);
        for nn = 1:length(template_SPM.SPM.xCon)
            contrast_name{nn} = template_SPM.SPM.xCon(nn).name;
            contrast_name{nn}(contrast_name{nn}==' ') = '_';
            contrast_name{nn}(contrast_name{nn}=='-') = '';
        end;  
                
        if ~exist(fullfile(pairsampleTtest_directory,analysis_type{aa}),'dir'); 
            mkdir(fullfile(pairsampleTtest_directory,analysis_type{aa}));    %Make GLM/Parametric etc directory within twosamplettest
        end
    
        if contrast_num > 6
            pairs_array = [1 3; 2 4; 7 9; 8 10];
            pairs_name = {'PreVsPost_TaskgtRest' 'PreVsPost_RestgtTask' 'PreVsPost_ME_Left' 'PreVsPost_ME_Right'};
        else 
            pairs_array = [1 3; 2 4];
            pairs_name = {'PreVsPost_TaskgtRest' 'PreVsPost_RestgtTask'};
        end
        
        
    for gg=1%:length(group)
        
        if ~exist(fullfile(pairsampleTtest_directory,analysis_type{aa}, group(gg).name),'dir');
        mkdir(fullfile(pairsampleTtest_directory,analysis_type{aa}, group(gg).name));    %Make directory for each set of subjects (e.g. 'Control')
        end;
        
    for cc=4:size(pairs_array,1)
        if ~exist(fullfile(pairsampleTtest_directory,analysis_type{aa}, group(gg).name, pairs_name{cc}),'dir');
        mkdir(fullfile(pairsampleTtest_directory,analysis_type{aa}, group(gg).name, pairs_name{cc}));    %Make directory for each paired t-test (e.g. Pre vs post task gt rest_
        end;
    
        load(fullfile(way,'second_level_pairedT_spm12'));
        
        output_directory = fullfile(pairsampleTtest_directory,analysis_type{aa}, group(gg).name, pairs_name{cc});
        
        matlabbatch{1}.spm.stats.factorial_design.dir = {output_directory};
        
        first_contrast_pair = ['0' num2str(pairs_array(cc,1))]; if pairs_array(cc,1) > 9; first_contrast_pair = first_contrast_pair(2:end); end;
        second_contrast_pair = ['0' num2str(pairs_array(cc,2))]; if pairs_array(cc,2) > 9; second_contrast_pair = second_contrast_pair(2:end); end;
        
    for ss=1:length(group(gg).subj)   
        
        insert_contrast1 = [contrast_begin first_contrast_pair '.nii'];
        insert_contrast2 = [contrast_begin second_contrast_pair '.nii'];
        
        insert_subj = [subj(group(gg).subj(ss)).path];
        
        subject_directory = fullfile(D, insert_subj, 'Stats', analysis_type{aa});
        
        P=cellstr(spm_select('FPList', subject_directory, insert_contrast1));
        matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(ss).scans{1,:} = P{1}; % new paired one
        clear P;
        
        P=cellstr(spm_select('FPList', subject_directory, insert_contrast2));
        matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(ss).scans{2,:} = P{1}; % new paired one
        clear P;

    end; %ss (inputting each subjects contrast)
    
    %Other Possible Inputs
    %matlabbatch{1}.spm.stats.factorial_design.cov.c;                    %Covariate vector
    %matlabbatch{1}.spm.stats.factorial_design.cov.cname;                %Covariate name
    %matlabbatch{1}.spm.stats.factorial_design.cov.iCFi;                 %Covariate interactions (default = none)
    %matlabbatch{1}.spm.stats.factorial_design.cov.iCC;                  %Covariate centering (default = overall mean)
    %matlabbatch{1}.spm.stats.factorial_design.multi_cov.files;          %Multiple covariate files
    %matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCFI;           %Multiple covariate interactions
    %matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCC;            %Multiple covariate centering

    %matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;   %Threshold masking (default = none)
    %matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;           %Implicit mask (default = yes [1])
    %matlabbatch{1}.spm.stats.factorial_design.masking.em = {};          %Explicit mask 

    %matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;           %Global calculation (default = omit; used for PET)
    %matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;   %Global normalisation (default = no grand mean scaling; used for PET)
    %matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;          %Global normalisation (default = none; used for PET)
    
 
    %keyboard                  %Used in bug-testing script   

    spm_jobman('run',matlabbatch);
    
    clear matlabbatch;
    
    load(fullfile(way,'estimate_spm12.mat'));
    stats_directory=fullfile(pairsampleTtest_directory, analysis_type{aa}, group(gg).name, pairs_name{cc});
    statsfile=fullfile(stats_directory,'SPM.mat');
    matlabbatch{1}.spm.stats.fmri_est.spmmat={statsfile};
    
    spm_jobman('run',matlabbatch);
    
    clear matlabbatch;
    
    load(fullfile(way,'contrast_manager_spm12.mat'));
    
    stats_directory=fullfile(pairsampleTtest_directory, analysis_type{aa}, group(gg).name, pairs_name{cc});
    statsfile=fullfile(stats_directory,'SPM.mat');
    matlabbatch{1}.spm.stats.con.spmmat={statsfile};
    
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Pre_gt_Post';         % t Contrast (f contrast = fcon)
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Post_gt_Pre';         % t Contrast (f contrast = fcon)
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
    
    spm_jobman('run',matlabbatch);
    
    end;    %cc (inputting each contrast (e.g. con_0001.img))
    
    end;    %gg (inputting each group type (e.g. ALL_11SJ))
    
    end;    %aa (inputting each analysis type (e.g. GLM))

end;    %end 'p' pair sample t-test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make model: 2nd level Factorial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'f')
    
    if ~exist(fullfile(withingroups_stats_directory,'Factorial'),'dir'); mkdir(fullfile(withingroups_stats_directory,'FlexibleFactorial')); end;
    factorial_directory = fullfile(withingroups_stats_directory,'Factorial');
    
    factor_1_name = 'Condition';
    factor_1_level = 4;
    factor_2_name = 'Session';
    factor_2_level = 3;
    contrast_num = factor_1_level*factor_2_level;
    level_array = [1 1 1 2 2 2 3 3 3 4 4 4; 1 2 3 1 2 3 1 2 3 1 2 3];
    
    for aa = 1%:length(analysis_type)            %change it if want just GLM, just Parametric, Just Bayesian etc
    
        if ~exist(fullfile(factorial_directory,analysis_type{aa}),'dir'); 
            mkdir(fullfile(factorial_directory,analysis_type{aa}));    %Make GLM/Parametric etc directory within onesamplettest
        end
            
    for gg=1:length(group)   
        
        if ~exist(fullfile(factorial_directory,analysis_type{aa}, group(gg).name),'dir');
        mkdir(fullfile(factorial_directory,analysis_type{aa}, group(gg).name));    %Make directory for each set of subjects (e.g. 'All_11SJ')
        end;
        

        load(fullfile(way,'second_level_factorial_analysis_spm12.mat'));
        
        output_directory = fullfile(factorial_directory,analysis_type{aa},group(gg).name);
        matlabbatch{1}.spm.stats.factorial_design.dir = {output_directory};
        
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = factor_1_name;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = factor_1_level;
        %matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;         %Independence
        %matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0;     %Variance Equal/Unequal
        %matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;        %Grand Mean Scaling (for PET/VBM)
        %matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;       %ANCOVA (for PET/VBM)
        
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = factor_2_name;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = factor_2_level;
        
        %matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;     %Generate contrasts: 1 = yes
    
    for cc = 1:contrast_num   
        
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cc).levels = level_array(:,cc);
        
    for ss=1:length(group(gg).subj)  
        
        contrast_end=['0' num2str(cc)]; if cc > 9; contrast_end=contrast_end(2:end); end;
        subj_end = ['0' num2str(group(gg).subj(ss))]; if length(subj_end) > 2; subj_end=subj_end(2:end); end;
        
        insert_contrast = [contrast_begin contrast_end '.nii'];
        insert_subj = [subj_begin subj_end];
        
        subject_directory = fullfile(D, insert_subj, 'Stats', analysis_type{aa});
        
        P=cellstr(spm_select('FPList', subject_directory, insert_contrast));
        
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell{ss} = P{1};
        
            clear P;

    end; %ss (inputting each subjects contrast)
        
    end; %cc (inputting particular contrast into cell)  
        
    %Other Possible Inputs
    %matlabbatch{1}.spm.stats.factorial_design.cov.c;                    %Covariate factor
    %matlabbatch{1}.spm.stats.factorial_design.cov.cname;                %Covariate name
    %matlabbatch{1}.spm.stats.factorial_design.cov.iCFi;                 %Covariate interactions (default = none)
    %matlabbatch{1}.spm.stats.factorial_design.cov.iCC;                  %Covariate centering (default = overall mean)
    %matlabbatch{1}.spm.stats.factorial_design.multi_cov.files;          %Multiple covariate files
    %matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCFI;           %Multiple covariate interactions
    %matlabbatch{1}.spm.stats.factorial_design.multi_cov.iCC;            %Multiple covariate centering

    %matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;   %Threshold masking (default = none)
    %matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;           %Implicit mask (default = yes [1])
    %matlabbatch{1}.spm.stats.factorial_design.masking.em = {};          %Explicit mask 

    %matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;           %Global calculation (default = omit; used for PET)
    %matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;   %Global normalisation (default = no grand mean scaling; used for PET)
    %matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;          %Global normalisation (default = none; used for PET)
    
 
    %keyboard                  %Used in bug-testing script   

    spm_jobman('run',matlabbatch);
    
    clear matlabbatch;
    
    load(fullfile(way,'estimate_spm12.mat'));
    stats_directory=fullfile(factorial_directory, analysis_type{aa}, group(gg).name);
    statsfile=fullfile(stats_directory,'SPM.mat');
    matlabbatch{1}.spm.stats.fmri_est.spmmat={statsfile};
    
    spm_jobman('run',matlabbatch);
    
    clear matlabbatch;
    

    end;    %gg (inputting each group type (e.g. ALL_11SJ))
    
    end;    %aa (inputting each analysis type (e.g. GLM))

end;    %end 'f' flexible factorial

return