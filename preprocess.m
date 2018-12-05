function preprocess(todo,D)
% Dec 2018 Redone for spm12
% Based on script by Adam McNamara, edited by Philip Dean

% INITIAL PREPROCESSING:
% Dicom transfer to get structural image (d: f*.nii & s*.nii)
% [Slice timing if necessary (a: af*.nii)]
% Realign and unwarp functional data (r: uf*.nii)
% Segment structural data (s: c1*.nii; c2*.nii; c3*.nii; ms*.nii; y_s*.nii)
% [Optional Skull-strip bias-corrected structural (ms*.nii) using imcalc]
% Coregister invididual structural data to mean functional and apply to all functional data (c: headers changed) 
% Normalise functional data using deformation field (y_s*nii) (n: wuf*.nii)
% Smooth normalised functional data (g: swuf*.nii)

% INPUT ARGUMENTS:

%%%%% 'todo'
% if todo not given, then it defaults to 'drsbcnog'
%           d = dicom transfer
%           a = slice timing (af*.nii)
%           r = realignment
%           s = segment
%           b = Skull strip bias-corrected Brain using Imcalc
%           c = coregistration to structural
%           n = Normalise (functional)
%           o = Normalise Structural
%           g = Smoothing with Gaussian kernal

% OTHER POSSIBLY USEFUL COMMANDS:
%           x = Delete uf*.img and wuf*.img to save disk space (only do if performed CheckReg to see if preprocessing OK) 
%           p = print movement parameters to pdf 
%           (NB these are already saved as .ps postscript files in format e.g. spm_2017Jan25.ps along with other preprocessing figures) 

%%%%% 'D'
% This is the Directory, e.g.  'E:\MRI\BECi_Study\Data\Subject_01'

% So could call script as:
%           preprocess('drsbcnog','E:\MRI\BECi_Study\Data\Subject_01')
% or, if just want to do dicom transfer: 
%           preprocess('d','E:\MRI\BECi_Study\Data\Subject_01')

% Global Variables
spm('Defaults', 'FMRI');        % Reset SPM defaults for fMRI (not sure necessary - safety catch?)
global defaults;                % Reset Global defaults (not sure why needed?)

if ~exist('todo','var'); todo='drsbcnog'; end;      % if nothing entered in "todo" bracket, then this is the default action

way='E:\MRI\BECi_Study\scripts\batch_files';    % Path to the "jobs"/batch files needed

TR = 3;                             % Bunched acquisition (2s acquire, 1s gap for EEG)             
nslices_fMRI = 25;                  % Number of slices
sliceorder = [];                    % Left blank here, but can be used to specify slice order for slice time correction
vxl_fmri = [3 3 3];                 % fMRI resolution 3x3x3 with 1mm gap (3x3x4)
vxl_str = [1 1 1];                  % Structural resolution

tic                                 % start clock timing how long analysis takes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dicom transfer     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'d')
    % if folders for each import dont exist, create these folders
    if ~exist(fullfile(D,'sess1'),'dir'); cd(D); mkdir('sess1'); end;     
%    if ~exist(fullfile(D,'sess2'),'dir'); cd(D); mkdir('sess2'); end;      % add these if more than one session
%    if ~exist(fullfile(D,'sess3'),'dir'); cd(D); mkdir('sess3'); end; 
    if ~exist(fullfile(D,'structural'),'dir'); cd(D); mkdir('structural'); end;

% use num_scans to get data to import (see nums_scans function below)
[T]=num_scans(D);

% import fMRI and structural
for tt=1:length(T)
  load(fullfile(way,'dicom_spm12.mat'));
    if length(T(tt).files) > 123;
      [cr,ap]=fileparts(D);  
        for ii=1:size(T(tt).files,2);
           matlabbatch{1}.spm.util.import.dicom.data{ii}=fullfile(fileparts(fileparts(D)),'Raw_Data',ap, T(tt).files{ii});%put the scans in
        end;
      fprintf('\n Number of Volumes = %d\n',size(matlabbatch{1}.spm.util.import.dicom.data,1));
      matlabbatch{1}.spm.util.import.dicom.outdir{1}=fullfile(D,T(tt).scantype);
% keyboard                  Used in bug-testing script
      spm_jobman('run',matlabbatch);
    end;
end;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slice-timing correction             NOT USED IN THIS ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'a')

load(fullfile(way,'slice_timing_spm12.mat'));

matlabbatch{1}.spm.temporal.st.nslices = nslices_fMRI;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TR - (TR/nslices_fMRI);
matlabbatch{1}.spm.temporal.st.so = [1:2:n_slices_fMRI 2:2:nslices_fMRI]; % interleaved bottom up
matlabbatch{1}.spm.temporal.st.refslice = 1; %Reference slice is first slice

P=cellstr(spm_select('FPList', directory,'^f.*\.nii$'));

for ii=1:size(P, 1);
      matlabbatch{1}.spm.temporal.st.scans{ii}=P{ii};%put the scans in
end;

fprintf('\nTR = %d',matlabbatch{1}.spm.temporal.st.tr);
fprintf('\nnslices = %d',matlabbatch{1}.spm.temporal.st.nslices);

% keyboard                  %Used in bug-testing script
spm_jobman('run',matlabbatch);

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Realignment: Realign & Unwarp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'r')
    
  load(fullfile(way,'realign_unwarp_spm12.mat'));           % change this to e.g. "realign_unwarp_2sess_spm12.mat" if more than one session
    
  for ss=1:3; %number of sessions
    directory=fullfile(D,['sess' num2str(ss)]);  

    P=cellstr(spm_select('FPList', directory,'^f.*\.nii$'));

    for ii=1:size(P, 1);
     matlabbatch{1}.spm.spatial.realignunwarp.data(ss).scans{ii} = P{ii};
    end
   
  clear P;
  end;
% keyboard                  %Used in bug-testing script
  spm_jobman('run',matlabbatch);
  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'s') 

  load(fullfile(way,'segment_biascorrected_spm12.mat'));
  struct_directory=fullfile(D,'structural');
  
  matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(spm_select('FPList', struct_directory,'^s.*\.nii$'));
  
% matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];   Save Bias Corrected Image
% matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];      Forward Deformation
% keyboard                  %Used in bug-testing script
  spm_jobman('run',matlabbatch);

end; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Skull Strip ImCalc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'b') 

  load(fullfile(way,'imcalc_brainstrip_spm12.mat'));
  struct_directory=fullfile(D,'structural');
  
  matlabbatch{1}.spm.util.imcalc.input(1) = cellstr(spm_select('FPList', struct_directory,'^c1s.*\.nii$'));  % i1: GM Segment
  matlabbatch{1}.spm.util.imcalc.input(2) = cellstr(spm_select('FPList', struct_directory,'^c2s.*\.nii$'));  % i2: WM Segment
  matlabbatch{1}.spm.util.imcalc.input(3) = cellstr(spm_select('FPList', struct_directory,'^c3s.*\.nii$'));  % i3: CSF Segment
  matlabbatch{1}.spm.util.imcalc.input(4) = cellstr(spm_select('FPList', struct_directory,'^ms.*\.nii$'));   % i4: Bias Corrected Image
  
  matlabbatch{1}.spm.util.imcalc.outdir = {struct_directory};     % Output directory
% keyboard                  %Used in bug-testing script
  spm_jobman('run',matlabbatch);
  
end;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coregistration: Estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'c') 

  load(fullfile(way,'coregister_est_spm12.mat'));

  struct_directory=fullfile(D,'structural');
  matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(spm_select('FPList', struct_directory,'^Brain.*\.nii$'));
  
  directory=fullfile(D,'sess1');
  matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(spm_select('FPList', directory,'^meanuf.*\.nii$'));
  
  P=cellstr(spm_select('FPList', directory,'^uf.*\.nii$'));        % input session 1 scans
 
%  directory=fullfile(D,'sess2');                                   % add these if more than one session
%  P=[P;[cellstr(spm_select('FPList', directory,'^uf.*\.nii$'))]];  % add session 2 scans
%  directory=fullfile(D,'sess3');
%  P=[P;[cellstr(spm_select('FPList', directory,'^uf.*\.nii$'))]];  % add session 3 scans
  
  for ii=1:size(P, 1);
     matlabbatch{1}.spm.spatial.coreg.estimate.other{ii} = P{ii};
  end

  clear P;
% keyboard                  %Used in bug-testing script
  spm_jobman('run',matlabbatch);
  
end;  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalization: Write
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'n')   
    
  load(fullfile(way,'normalise_write_333_spm12.mat'));                  % Change batch file if resolution not 3x3x3 or change in script below

% matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3]     Change to [1 1 1] for structural
  
  struct_directory=fullfile(D,'structural');
  matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList', struct_directory,'^y_s.*\.nii$'));
  
  directory=fullfile(D,'sess1');
  P=cellstr(spm_select('FPList', directory,'^uf.*\.nii$'));        % input session 1 scans
  
%  directory=fullfile(D,'sess2');                                   % add these if more than one session
%  P=[P;[cellstr(spm_select('FPList', directory,'^uf.*\.nii$'))]];  % add session 2 scans
%  directory=fullfile(D,'sess3');
%  P=[P;[cellstr(spm_select('FPList', directory,'^uf.*\.nii$'))]];  % add session 3 scans

  
  for ii=1:size(P, 1);
     matlabbatch{1}.spm.spatial.normalise.write.subj.resample{ii} = P{ii};
  end

  clear P;
% keyboard                  %Used in bug-testing script
  spm_jobman('run',matlabbatch);
  
end;

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalization (structural): Write
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'o')   
    
  load(fullfile(way,'normalise_write_111_spm12.mat'));                  % Change batch file if resolution not 3x3x3 or change in script below

% matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1]     Change to [3 3 3] for functional

  struct_directory=fullfile(D,'structural');
  matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList', struct_directory,'^y_s.*\.nii$'));
  matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(spm_select('FPList', struct_directory,'^Brain.*\.nii$'));
% keyboard                  Used in bug-testing script 
  spm_jobman('run',matlabbatch);
  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'g')
   
    load(fullfile(way,'smooth_spm12.mat'));

% Other possible things to change:
% matlabbatch{1}.spm.spatial.smooth.fhwm = [8 8 8]      % Smoothing level
% matlabbatch{1}.spm.spatial.smooth.im = 0              % Implicit mask (0 = no, 1 = yes)

    
  directory=fullfile(D,'sess1');                                    % add these if more than one session
  P=cellstr(spm_select('FPList', directory,'^wuf.*\.nii$'));        % input session 1 scans
  
%  directory=fullfile(D,'sess2');
%  P=[P;[cellstr(spm_select('FPList', directory,'^wuf.*\.nii$'))]];  % add session 2 scans
%  directory=fullfile(D,'sess3');
%  P=[P;[cellstr(spm_select('FPList', directory,'^wuf.*\.nii$'))]];  % add session 3 scans

  
  for ii=1:size(P, 1);
     matlabbatch{1}.spm.spatial.smooth.data{ii} = P{ii};
  end

  clear P;
% keyboard                  %Used in bug-testing script   
  spm_jobman('run',matlabbatch);
  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Delete uf*.im and wuf*.img files to save space
% WARNING: DELETES FILES PERMANENTLY - TAKE CARE!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'x')   
 
 % Preprocessing files:
 % f*.nii             KEEP
 % uf*.nii
 %      rp_f*.txt     KEEP
 %      meanuf*.nii   KEEP
 % wuf*.nii
 % swuf*.nii          KEEP

 here = pwd;            %remember where I am before I start to work
 
 directory_sess1=fullfile(D,'sess1');
 %directory_sess2=fullfile(D,'sess2');
 %directory_sess3=fullfile(D,'sess3');

% keyboard                  %Used in bug-testing script  
 
 cd(directory_sess1);
 delete ('uf*.*')       % to save diskspace (take care!)
 delete ('wuf*.*')      % to save diskspace (take care!)

% cd(directory_sess2);
% delete ('uf*.*')       % to save diskspace (take care!)
% delete ('wuf*.*')      % to save diskspace (take care!)
% cd(directory_sess3);
% delete ('uf*.*')       % to save diskspace (take care!)
% delete ('wuf*.*')      % to save diskspace (take care!)
 
 cd (here);             % go back to where I was
 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Graph and Print Movement Regressors to PDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strfind(todo,'p')   
 
 here = pwd;            %remember where I am before I start to work
 
 directory_sess1=fullfile(D,'sess1');
 %directory_sess2=fullfile(D,'sess2');
 %directory_sess3=fullfile(D,'sess3');
 
 [path,sj] = fileparts(D);

 movesess1=dlmread(spm_select('FPList', directory_sess1, '^rp_f.*\.txt$'));
% movesess2=dlmread(spm_select('FPList', directory_sess2, '^rp_f.*\.txt$'));
% movesess3=dlmread(spm_select('FPList', directory_sess3, '^rp_f.*\.txt$'));
 
 scrsz = get(groot,'ScreenSize');
 figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)]);
 
 trans1 = subplot(2,1,1);
    plot(movesess1(:,1), 'b');
    hold on; 
    plot(movesess1(:,2), 'g');
    plot(movesess1(:,3), 'r');
    title([sj ' Translation']);
    xlabel ('image sess1');
    ylabel ('mm');
    legend ('x', 'y', 'z', 'Location', 'northeastoutside');
 
% trans2 = subplot(3,2,3);
%    plot(movesess2(:,1), 'b');
%    hold on; 
%    plot(movesess2(:,2), 'g');
%    plot(movesess2(:,3), 'r');
%    title([sj ' Translation']);
%    xlabel ('image sess2');
%    ylabel ('mm');
%    legend ('x', 'y', 'z', 'Location', 'northeastoutside');
%  trans3 = subplot(3,2,5);
%    plot(movesess3(:,1), 'b');
%    hold on; 
%    plot(movesess3(:,2), 'g');
%    plot(movesess3(:,3), 'r');
%    title([sj ' Translation']);
%    xlabel ('image sess3');
%    ylabel ('mm');
%    legend ('x', 'y', 'z', 'Location', 'northeastoutside');
    
  rot1 = subplot(2,1,2);
    plot(movesess1(:,4), 'b');
    hold on; 
    plot(movesess1(:,5), 'g');
    plot(movesess1(:,6), 'r');
    title([sj ' Rotation']);
    xlabel ('image sess1');
    ylabel ('degrees');
    legend ('pitch', 'roll', 'yaw', 'Location', 'northeastoutside');  
  
%  rot2 = subplot(3,2,4);
%    plot(movesess2(:,4), 'b');
%    hold on; 
%    plot(movesess2(:,5), 'g');
%    plot(movesess2(:,6), 'r');
%    title([sj ' Rotation']);
%    xlabel ('image sess2');
%    ylabel ('degrees');
%    legend ('pitch', 'roll', 'yaw', 'Location', 'northeastoutside');
%  rot3 = subplot(3,2,6);
%    plot(movesess3(:,4), 'b');
%    hold on; 
%    plot(movesess3(:,5), 'g');
%    plot(movesess3(:,6), 'r');
%    title([sj ' Rotation']);
%    xlabel ('image sess3');
%    ylabel ('degrees');
%    legend ('pitch', 'roll', 'yaw', 'Location', 'northeastoutside');

%    maxmin_mm_values = [max(movesess1(:,1)), min(movesess1(:,1)), max(movesess1(:,2)), min(movesess1(:,2)), max(movesess1(:,3)), min(movesess1(:,3)); 
%    max(movesess2(:,1)), min(movesess2(:,1)), max(movesess2(:,2)), min(movesess2(:,2)), max(movesess2(:,3)), min(movesess2(:,3));
%    max(movesess3(:,1)), min(movesess3(:,1)), max(movesess3(:,2)), min(movesess3(:,2)), max(movesess3(:,3)), min(movesess3(:,3))];
 
    maxmin_mm_values = [max(movesess1(:,1)), min(movesess1(:,1)), max(movesess1(:,2)), min(movesess1(:,2)), max(movesess1(:,3)), min(movesess1(:,3))];
 

  cd (D);
  
  print([sj '_Head_Movement'],'-dpdf','-fillpage');
  csvwrite([sj '_Head_Movement_MaxMin.csv'], maxmin_mm_values);
  
 cd (here);             % go back to where I was
 
end;

   toc                                          % stop clock timing how long analysis takes
return  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OTHER FUNCTIONS USED BY SCRIPT ABOVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [T]=num_scans(D)
% For DICOM Import
% Returns most likely sorting procedure for your sessions
% Looks into Raw Data, and classifies what data you have dependent on the number of scans in it
% Looks in more than one folder if necessary (multi-session, 1 folder only protocol below)

% In this case:
% Scan Number: 3 = Localiser
% Scan Number: 176 = structural
% Scan Number: >176 = fMRI Data (sess1)

% These labels also used to create folder in Dicom Import as needed
% Change dependent on your setup

cd(D);
[cr,ap]=fileparts(pwd);

%%%%%% Look for Raw Data
d=dir(fullfile('..\..\Raw_Data',ap));       %This assumes you have a setup as described in the "multisubject_analysis" file

scan_mem_d=0;
c_d=0;

for jj=1:20;T(jj).files={};end;

for ss = 3:length(d);
  if strcmp(d(ss).name(end-3:end),'.IMA');
      f=find(double(d(ss).name) == 46);
      scan_d=str2num(d(ss).name(f(3)+1:f(4)-1));
  if scan_mem_d ~= scan_d; c_d=1; scan_mem_d=scan_d; else; c_d= c_d+1; end;
      T(scan_d).files{c_d}=d(ss).name;
  end    
end;
      
c_d=1; 
for jj=1:length(T); 
  if ~isempty(T(jj).files); 
      t(c_d)=T(jj); c_d=c_d+1; 
  end; 
end; 

T=t;
c_d=1;
for jj=1:length(T); 
    if length(T(jj).files) == 3; T(jj).scantype='localizer';  end;
    if length(T(jj).files) == 176; T(jj).scantype='structural';  end;
    if length(T(jj).files) > 176; T(jj).scantype=['sess' num2str(c_d)]; c_d=c_d+1;  end;
end