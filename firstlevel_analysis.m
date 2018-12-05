function firstlevel_analysis(todo,D, sh, subj)
% Dec 2018 redone for spm12
% Based on script by Adam McNamara, edited by Philip Dean

% MAKE AND ESTIMATE MODEL
% Make model using experimental logfiles/timings
% Estimate model to create SPM file for first level analysis
% Create contrasts
% Results manager

% INPUT ARGUMENTS:

% 'todo'
% m = make model (Basic GLM)
% p = make model (Parametric)
% b = make model (Bayesian)
% e = estimate model    NB: Need to change this dependent on what model made (see above)
% c = contrasts manager NB: Need to change this dependent on what model made (see above)
% r = results report    NB: Need to change this dependent on what model made (see above)

% if todo not given defaults todo = 'me'

% 'D'
% This is the Directory, e.g.  'E:\MRI\BECi_Study\Data\Subject_01'

% So could call script as:
% firstlevel_analysis('me','E:\MRI\BECi_Study\Data\Subject_01')
% or, if just want to do model: 
% firstlevel_analysis('m','E:\MRI\BECi_Study\Data\Subject_01')

% 'sh', 'subj' 
% These are inputs from the "multisubject_analysis.m" script
% sh constains: sh.studypath; sh.imagepath; sh.behavpath
% subj contains: subj.path; subj.task; subj.log_sess1; subj.log_sess2; subj.log_sess3; subj.response_button

% Global Variables
spm('Defaults', 'FMRI');        % Reset SPM defaults for fMRI (not sure necessary - safety catch?)
global defaults;                % Reset Global defaults (not sure why needed?)

if ~exist('todo','var'); todo='me'; end;             % if nothing entered in "todo" bracket, then this is the default action

way='E:\MRI\BECi_Study\scripts\batch_files';            % Path to the "jobs"/batch files needed

TR = 3;                             % Bunched acqusition (2s acquire, 1s gap for EEG)             
nslices_fMRI = 25;                  % Number of slices
M = [0 0 0 0 0 0];                  % Movement Parameter in-fill for contrasts
numsess = 1;                        % Number of sessions

% Contrasts to be used. Design is:
% M is movement regressors, C is the constant for each session
% [Rest_sess1 0Back_sess1 2Back_sess1 4-Back_sess1 M Rest_sess2 etc then Rest_sess3 etc then C C C] 
% F Contrasts
% T Contrasts


tic                                 % start clock timing how long analysis takes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make model: Basic GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'m')

    % if folder for stats doesnt exist, create the folder
    if ~exist(fullfile(D,'Stats\GLM'),'dir'); cd(D); mkdir('Stats\GLM'); end;     
    stats_directory=fullfile(D,'Stats\GLM');
    
    load(fullfile(way,'first_level_spm12.mat'));
    
    %%%%%% Global setup for model
    matlabbatch{1}.spm.stats.fmri_spec.dir = {stats_directory}; %Output Directory
    %matlabbatch{1}.spm.stats.fmri_spec.timing.units = secs;    %TIMING secs/scans
    %matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 3;          %TR
    %matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;     %Microtime resolution
    %matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;     %Microtime Onset
    
    
    %%%%%% Filling in scans and conditions for each session
    for ss=1:numsess; %number of sessions
        directory=fullfile(D,['sess' num2str(ss)]);
        
        % Loading up pre-processed scans 
        P=cellstr(spm_select('FPList', directory,'^swuf.*\.nii$'));

            for ii=1:size(P, 1);
                matlabbatch{1}.spm.stats.fmri_spec.sess(ss).scans{ii} = P{ii};
            end
   
        clear P;
  
        
    %%%%%% Loading up onset/duration from behavioural logfile using script "behavioural_data_onset_duration" 
        global ev;
        ev = {};
      
        %if ss == 1;
            ev = behavioural_data_onset_duration((fullfile(sh.behavpath, subj.path, ['sess' num2str(ss)], subj.log_sess1)), subj.response_button);   %finding the behavioural file for participant
        %elseif ss == 2;
        %    ev = behavioural_data_onset_duration((fullfile(sh.behavpath, subj.path, ['sess' num2str(ss)], subj.log_sess2)), subj.response_button);   
        %else
        %    ev = behavioural_data_onset_duration((fullfile(sh.behavpath, subj.path, ['sess' num2str(ss)], subj.log_sess3)), subj.response_button);   
        %end
        
         % Onsets/Duration 
        n0=ev.blockstartduration_task_adjusted(find(ev.blockstartduration_task_adjusted(:,2) == 2),:);
        n2=ev.blockstartduration_task_adjusted(find(ev.blockstartduration_task_adjusted(:,2) == 3),:);
        n4=ev.blockstartduration_task_adjusted(find(ev.blockstartduration_task_adjusted(:,2) == 4),:);
    
   %%%%%% Filling in onsets and durations for each condition
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(1).name = 'Rest';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(1).onset =(ev.blockstartduration_rest(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(1).duration =(ev.blockstartduration_rest(:,end)./10000);
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).name = '0-Back';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).onset = (n0(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).duration = (n0(:,end)./10000);
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).name = '2-Back';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).onset = (n2(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).duration = (n2(:,end)./10000);
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).name = '4-Back';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).onset = (n4(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).duration = (n4(:,end)./10000);
    
     %%%%%% Loading in movement parameters as multiple regressor
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).multi_reg = cellstr(spm_select('FPList', directory,'^rp_f.*\.txt$'));
    
     %%%%%% Other possible inputs to Session:
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;       % Time Modulation
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod;           % Parametric Modulation
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;       % Orthogonalise Modulations (1=Yes)
    
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};                  % For Multiple Conditions (.mat file)             
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = [];                % Regressors to be regressed out of data             
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;              % High Pass filter (for filtering out scanner drift/noise)
    
    end;
    
     %%%%%% Other Possible inputs
    %matlabbatch{1}.spm.stats.fmri_spec.fact = [];                       % Factorial Design
    %matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];   % Basis Function & Derivatives
    %matlabbatch{1}.spm.stats.fmri_spec.volt = 1;                   % Model Interactions(Volterra)
    %matlabbatch{1}.spm.stats.fmri_spec.global = 'None';              % Global Normalisation
    %matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;              % Masking Threshold
    %matlabbatch{1}.spm.stats.fmri_spec.mask = {''};                       % Explicit Mask
    %matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';                % Serial Correlations

% keyboard                  %Used in bug-testing script   
  spm_jobman('run',matlabbatch);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make model: Parametric 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'p')

    % if folder for stats doesnt exist, create the folder
    if ~exist(fullfile(D,'Stats\Parametric'),'dir'); cd(D); mkdir('Stats\Parametric'); end;     
    stats_directory=fullfile(D,'Stats\Parametric');
    
    load(fullfile(way,'first_level_spm12.mat'));
    
    %%%%%% Global setup for model
    matlabbatch{1}.spm.stats.fmri_spec.dir = {stats_directory}; %Output Directory
    %matlabbatch{1}.spm.stats.fmri_spec.timing.units = secs;    %TIMING secs/scans
    %matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 3;          %TR
    %matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;     %Microtime resolution
    %matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;     %Microtime Onset
    
    
    %%%%%% Filling in scans and conditions for each session
    for ss=1:numsess; %number of sessions
        directory=fullfile(D,['sess' num2str(ss)]);
        
        % Loading up pre-processed scans 
        P=cellstr(spm_select('FPList', directory,'^swuf.*\.nii$'));

            for ii=1:size(P, 1);
                matlabbatch{1}.spm.stats.fmri_spec.sess(ss).scans{ii} = P{ii};
            end
   
        clear P;
  
        
     %%%%%% Loading up onset/duration and precent correct (PARAMETRIC) from behavioural logfile using script "behavioural_data_onset_duration"
        global ev;
        ev = {};
      
        %if ss == 1;
            ev = behavioural_data_onset_duration((fullfile(sh.behavpath, subj.path, ['sess' num2str(ss)], subj.log_sess1)), subj.response_button);   %finding the behavioural file for participant
        %elseif ss == 2;
        %    ev = behavioural_data_onset_duration((fullfile(sh.behavpath, subj.path, ['sess' num2str(ss)], subj.log_sess2)), subj.response_button);   
        %else
        %    ev = behavioural_data_onset_duration((fullfile(sh.behavpath, subj.path, ['sess' num2str(ss)], subj.log_sess3)), subj.response_button);   
        %end
        
        % Onsets/Duration 
        n0=ev.blockstartduration_task_adjusted(find(ev.blockstartduration_task_adjusted(:,2) == 2),:);
        n2=ev.blockstartduration_task_adjusted(find(ev.blockstartduration_task_adjusted(:,2) == 3),:);
        n4=ev.blockstartduration_task_adjusted(find(ev.blockstartduration_task_adjusted(:,2) == 4),:);
        % Percent Correct 
        p0=ev.percent_correct(find(ev.percent_correct(:,2) == 2),:);
        p2=ev.percent_correct(find(ev.percent_correct(:,2) == 3),:);
        p4=ev.percent_correct(find(ev.percent_correct(:,2) == 4),:);
    
     %%%%%% Filling in onsets and durations for each condition   
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(1).name = 'Rest';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(1).onset =(ev.blockstartduration_rest(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(1).duration =(ev.blockstartduration_rest(:,end)./10000);
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).name = '0-Back';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).onset = (n0(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).duration = (n0(:,end)./10000);
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).name = '2-Back';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).onset = (n2(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).duration = (n2(:,end)./10000);
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).name = '4-Back';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).onset = (n4(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).duration = (n4(:,end)./10000);
    
    %%%%%% Filling in percent correct for each task block (PARAMETRIC, Total % Correct)
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).pmod.name = '0-Back Correct';   % Parametric Modulation
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).pmod.poly = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).pmod.param = p0(:,end-2);           
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).pmod.name = '2-Back Correct';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).pmod.poly = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).pmod.param = p2(:,end-2);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).pmod.name = '4-Back Correct';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).pmod.poly = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).pmod.param = p4(:,end-2);
        
     %%%%%% Loading in movement parameters as multiple regressor
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).multi_reg = cellstr(spm_select('FPList', directory,'^rp_f.*\.txt$'));
    
     %%%%%% Other possible inputs to Session:
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;       % Time Modulation
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;       % Orthogonalise Modulations (1=Yes)
    
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};                  % For Multiple Conditions (.mat file)             
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = [];                % Regressors to be regressed out of data             
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;              % High Pass filter (for filtering out scanner drift/noise)
    
    end;
    
     %%%%%% Other Possible inputs
    %matlabbatch{1}.spm.stats.fmri_spec.fact = [];                       % Factorial Design
    %matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];   % Basis Function & Derivatives
    %matlabbatch{1}.spm.stats.fmri_spec.volt = 1;                   % Model Interactions(Volterra)
    %matlabbatch{1}.spm.stats.fmri_spec.global = 'None';              % Global Normalisation
    %matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;              % Masking Threshold
    %matlabbatch{1}.spm.stats.fmri_spec.mask = {''};                       % Explicit Mask
    %matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';                % Serial Correlations

% keyboard                  %Used in bug-testing script   
  spm_jobman('run',matlabbatch);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make model: Bayesian 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'b')
    
    % if folder for stats doesnt exist, create the folder
    if ~exist(fullfile(D,'Stats\Bayesian'),'dir'); cd(D); mkdir('Stats\Bayesian'); end;     
    stats_directory=fullfile(D,'Stats\Bayesian');
    
    load(fullfile(way,'first_level_spm12.mat'));

    %%%%%% Global setup for model
    matlabbatch{1}.spm.stats.fmri_spec.dir = {stats_directory}; %Output Directory
    %matlabbatch{1}.spm.stats.fmri_spec.timing.units = secs;    %TIMING secs/scans
    %matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 3;          %TR
    %matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;     %Microtime resolution
    %matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;     %Microtime Onset
    
    
    %%%%%% Filling in scans and conditions for each session
    for ss=1:numsess; %number of sessions
        directory=fullfile(D,['sess' num2str(ss)]);
        
        % Loading up pre-processed scans 
        P=cellstr(spm_select('FPList', directory,'^wuf.*\.nii$'));      % NON-SMOOTHED DATA USED

            for ii=1:size(P, 1);
                matlabbatch{1}.spm.stats.fmri_spec.sess(ss).scans{ii} = P{ii};
            end
   
        clear P;
  
        
      %%%%%% Loading up onset/duration and precent correct (PARAMETRIC) from behavioural logfile using script "behavioural_data_onset_duration"
        global ev;
        ev = {};
      
        %if ss == 1;
            ev = behavioural_data_onset_duration((fullfile(sh.behavpath, subj.path, ['sess' num2str(ss)], subj.log_sess1)), subj.response_button);   %finding the behavioural file for participant
        %elseif ss == 2;
        %    ev = behavioural_data_onset_duration((fullfile(sh.behavpath, subj.path, ['sess' num2str(ss)], subj.log_sess2)), subj.response_button);   
        %else
        %    ev = behavioural_data_onset_duration((fullfile(sh.behavpath, subj.path, ['sess' num2str(ss)], subj.log_sess3)), subj.response_button);   
        %end
        
        % Onsets/Duration 
        n0=ev.blockstartduration_task_adjusted(find(ev.blockstartduration_task_adjusted(:,2) == 2),:);
        n2=ev.blockstartduration_task_adjusted(find(ev.blockstartduration_task_adjusted(:,2) == 3),:);
        n4=ev.blockstartduration_task_adjusted(find(ev.blockstartduration_task_adjusted(:,2) == 4),:);
    
     %%%%%% Filling in onsets and durations for each condition     
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(1).name = 'Rest';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(1).onset =(ev.blockstartduration_rest(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(1).duration =(ev.blockstartduration_rest(:,end)./10000);
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).name = '0-Back';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).onset = (n0(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(2).duration = (n0(:,end)./10000);
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).name = '2-Back';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).onset = (n2(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(3).duration = (n2(:,end)./10000);
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).name = '4-Back';
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).onset = (n4(:,end-1)./10000);
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).cond(4).duration = (n4(:,end)./10000);
    
     %%%%%% Loading in movement parameters as multiple regressor
        matlabbatch{1}.spm.stats.fmri_spec.sess(ss).multi_reg = cellstr(spm_select('FPList', directory,'^rp_f.*\.txt$'));
    
     %%%%%% Other possible inputs to Session:
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;       % Time Modulation
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod;           % Parametric Modulation
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;       % Orthogonalise Modulations (1=Yes)
    
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};           % For Multiple Conditions (.mat file)             
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = [];           % Regressors to be regressed out of data             
        %matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;              % High Pass filter (for filtering out scanner drift/noise)
    
    end;
    
     %%%%%% Other Possible inputs
    %matlabbatch{1}.spm.stats.fmri_spec.fact = [];                  % Factorial Design
    %matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];   % Basis Function & Derivatives
    %matlabbatch{1}.spm.stats.fmri_spec.volt = 1;                   % Model Interactions(Volterra)
    %matlabbatch{1}.spm.stats.fmri_spec.global = 'None';            % Global Normalisation
    %matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;              % Masking Threshold
    %matlabbatch{1}.spm.stats.fmri_spec.mask = {''};                % Explicit Mask
    %matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';              % Serial Correlations
    
    
% keyboard                  %Used in bug-testing script   
  spm_jobman('run',matlabbatch);

end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'e')
    
    stats_directory=fullfile(D,'Stats\GLM');
    %stats_directory=fullfile(D,'Stats\Parametric');            % change these dependent on what analysis you have done
    %stats_directory=fullfile(D,'Stats\Bayesian');              % NEEDS BAYESIAN ESTIMATION
    
    statsfile=fullfile(stats_directory,'SPM.mat');
    
    load(fullfile(way,'estimate_spm12.mat'));
    %load(fullfile(way,'estimate_bayesian_spm12.mat'));         % BAYESIAN ESTIMATION
        
    %What SPM.mat file to open and estimate
    matlabbatch{1}.spm.stats.fmri_est.spmmat={statsfile};

    %%%%% Other Options
    %matlabbatch{1}.spm.stats.fmri_est.write_residuals=0;       % Write Residuals: 0=No
    %matlabbatch{1}.spm.stats.fmri_est.method.Classical=1;      % Method = classical (opposed to Bayesian)
    
    %%%%% Bayesian Options
    %matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.space.volume.block_type = 'Slices';
    %matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.signal = 'UGL';
    %matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.ARP = 3;
    %matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.noise.UGL = 1;
    %matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.LogEv = 'No';
    %matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.anova.first = 'No';
    %matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.anova.second = 'Yes';
    %matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.gcon = ; Simple contrasts. name & convec
    
% keyboard                  Used in bug-testing script   
  spm_jobman('run',matlabbatch);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contrast Manager
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'c')
    
    stats_directory=fullfile(D,'Stats\GLM');
    statsfile=fullfile(stats_directory,'SPM.mat');
    
    load(fullfile(way,'contrast_manager_spm12.mat'));
        
    %What SPM.mat file to open
    matlabbatch{1}.spm.stats.con.spmmat={statsfile};
    
    % What contrasts to put in
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'ME_Rest';         % t Contrast (f contrast = fcon)
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'sess';       
    %Reproduce across sessions: 'none' dont replicate; 'sess' create per session; 'repl' replicate; 'both' replicate and create
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'ME_0Back';         % t Contrast
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'sess';     
    
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'ME_2Back';         % t Contrast
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'sess';
    
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'ME_4Back';         % t Contrast
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1];
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'sess'; 
    
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = '0Back_gt_R';       % t Contrast
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [-1 1 0 0];
    matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'sess'; 
    
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = '2Back_gt_R';       % t Contrast
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [-1 0 1 0];
    matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'sess';
    
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = '4Back_gt_R';       % t Contrast
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [-1 0 0 1];
    matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'sess';
    
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = '2Back_gt_0Back';   % t Contrast
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 -1 1 0];
    matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'sess'; 
    
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = '4Back_gt_0Back';   % t Contrast
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [0 -1 0 1];
    matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'sess'; 
    
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = '0Back_gt_2Back';  % t Contrast
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [0 1 -1 0];
    matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'sess'; 
    
    matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = '0Back_gt_4Back';  % t Contrast
    matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [0 1 0 -1];
    matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'sess';
    
    matlabbatch{1}.spm.stats.con.delete = 1;        % delete existing contrasts (0=No, 1=Yes)
    
    % Contrast in terms of condition or regressor instead of columns in
    % design matrix. Allows creation of contrast automatically even if some
    % columns not always present (e.g. parametric modulations)
    %matlabbatch{1}.spm.stats.con.consess{3}.tconsess.name - '';     % T Contrast (cond/sess based)
    %matlabbatch{1}.spm.stats.con.consess{3}.tconsess.sessions = [1 2 3];    % Which sessions this contrast should be over
    %matlabbatch{1}.spm.stats.con.consess{3}.tconsess.coltype.colconds.conweight = [];
    %matlabbatch{1}.spm.stats.con.consess{3}.tconsess.coltype.colconds.colcond = [];
    %matlabbatch{1}.spm.stats.con.consess{3}.tconsess.coltype.colconds.colbf = [];
    %matlabbatch{1}.spm.stats.con.consess{3}.tconsess.coltype.colconds.colmod = [];
    %matlabbatch{1}.spm.stats.con.consess{3}.tconsess.coltype.colconds.colmodord = [];
    
% keyboard                  Used in bug-testing script   
  spm_jobman('run',matlabbatch);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results Report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'r')
    
    stats_directory=fullfile(D,'Stats\GLM');
    statsfile=fullfile(stats_directory,'SPM.mat');
    
    threshtype = 'FWE';     % Test type
    p_value = 0.0500;       % P Value
    vox_extent = 0;         % Number of voxels per cluster
    
    load(fullfile(way,'results_report_spm12.mat'));
    
    %What SPM.mat file to open
    matlabbatch{1}.spm.stats.results.spmmat={statsfile};

    % Which contrast files to report
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 1 0Back gt Rest';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 13;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype;    % Test type
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;           % P Value
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;        % Number of voxels per cluster
    %matlabbatch{1}.spm.stats.results.conspec(1).conjunction = 1;   
    %matlabbatch{1}.spm.stats.results.conspec(1).mask.none = 1;              % No masking
    
    %matlabbatch{1}.spm.stats.results.units = 1;            % 1=Volumetric 2D/3D; 
    %matlabbatch{1}.spm.stats.results.export{1}.ps = 1;     % Export results as postscript
    %Can also be exported as: Thresholded SPM, All Clusters (binary/n-ary),
    %eps, pdf, jpeg, png, tiff, Matlab figure, CSV file, Excel spreadsheet,
    %NIDM
    
 %keyboard                 % Used in bug-testing script   
  spm_jobman('run',matlabbatch);
  
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 2 0Back gt Rest';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 14;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;

  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 3 0Back gt Rest';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 15;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;
    
  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 1 2Back gt Rest';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 16;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent; 
    
  spm_jobman('run',matlabbatch);
        
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 2 2Back gt Rest';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 17;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;
    
  spm_jobman('run',matlabbatch);
        
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 3 2Back gt Rest';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 18;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;
    
  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 1 4Back gt Rest';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 19;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;
    
  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 2 4Back gt Rest';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 20;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;
    
  spm_jobman('run',matlabbatch);
        
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 3 4Back gt Rest';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 21;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;
    
  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 1 2Back gt 0Back';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 22;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent; 
    
  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 2 2Back gt 0Back';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 23;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;
    
  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 3 2Back gt 0Back';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 24;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent; 
    
  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 1 4Back gt 0Back';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 25;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent; 
    
  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 2 4Back gt 0Back';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 26;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;
    
  spm_jobman('run',matlabbatch);
    
    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = 'Session 3 4Back gt 0Back';
    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = 27;
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = threshtype; 
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = p_value;    
    matlabbatch{1}.spm.stats.results.conspec(1).extent = vox_extent;

  spm_jobman('run',matlabbatch);
    
  
end 

   toc                                          % stop clock timing how long analysis takes
return  

