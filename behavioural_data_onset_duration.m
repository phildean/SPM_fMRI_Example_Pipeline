function [ev] = behavioural_data_onset_duration(filename,response_button)
% Updated Dec 2018
% Based on script by Adam McNamara, edited by Philip Dean

% READS BEHAVIOURAL FILE TO GET TIMING OF BLOCKS
% (NB: This is for block design fMRI study. Event related needs this for
% each individual stimuli)


% First level analysis needs block start and block duration (or event start/duration) relative to start of fMRI scan.

% Also find timing of other factors, especially rest screen (and intro
% screen, feedback, stimuli and response)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define variables to be used in script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ev = {};                            % This is the ev variable that is reported out by the script

%%%% These subvariable arrays are made when going through the logfile below
ev.blockstartend = [];              % This is a sub-variable in ev which notes start and end of block: [Block, Block_Type, StartEnd, Time]
ev.stimulus = [];                   % Stimuli in Block: [Block, Block_Type, Actual_Number, Is_Target, Time]
ev.stimulus_other = [];             % Other Stimuli (outside block): [Block, Block_Type, Stimuli_Type, Time]
ev.response = [];                   % Response in Block: [Block, Exp_Block, Block_Type, Response_Code, Is_Target, Time]
ev.response_other = [];             % Response outside block: [Block, Block_Type, Response_Code, Time]
ev.pulse_timing = [];               % Timing of MRI Pulse: [Block, Block_Type, Time]

%%%% These are variables used in making the arrays above
Block = 0;                      % Counts up number of blocks
Exp_Block = 0;                  % Counts up number of task (non-rest) blocks. 
BlockType = 0;                  % 1=Rest; 2=0-Back; 3=2-Back; 4=4-Back;  0=Beginning of experiment
startend = 2;                   % 0=In Block; 1=End of Block & Between Block; 2=Beginning of Experiment
first_time_only = 0;            % Make sure that only counts first MRI pulse as start of session
stim_other_type = 0;            % 1=nBack Intro Screen, 2=Rest Intro Screen, 3=Intro Rest, 4=Intro 0-Back, 5=Intro 2-Back, 6=Intro 4-Back

response_on = 0;                % Determines whether to record response and prevents a double response to one stimuli 0=record response; 1=don't record response (response outside Block, and prevents a double response to one stimuli) 
is_target = 0;                  % Whether stimuli is "target" or not: 0=non-target; 1=target

%%%% These subvariable arrays created from above arrays after going through logfile
ev.blockstartduration = [];         % This is in the format onset/duration as needed by 1st level analysis [Block, Block_Type, StartEnd, Time, Duration]
ev.blockstartduration_rest = [];    % [Block, Block_Type, StartEnd, Time, Duration]
ev.blockstartduration_task = [];    % [Block, Block_Type, StartEnd, Time, Duration]
ev.blockstartduration_task_adjusted = [];% Adjusted for start of block (time_adjust_2 below) 

ev.percent_correct = [];            % USED IN FIRST LEVEL PARAMETRIC ANALYSIS [Block, Block_Type, Total % Correct, Target % Correct, Non-Target % Correct]
singleblock = [];               % used in percent correct calculation

%%%% These are variables used in making the arrays above
num_correct_T = 0;              % Number Correct: Target
num_correct_N = 0;              % Number Correct: Non-Target
num_correct = 0;                % Total Number Correct

%%%% Constants for your experiment type. used in percent correct calculation
trial_per_block = 15;
target_stimuli_block = 5;
nontarget_stimuli_block = 10;

time_adjust_1 = 60000;              % seconds between intro screen and start of block in first experiment
time_adjust_2 = 21500;           % seconds between Pulse MR screen and start of block in main experiment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename);          % opens the logfile

begin = 0;                      %identifies first trial and switches to collecting data

%%%%%%%%%%%% Code works out if at end of logfile by finding 20 consecutive "gaps"
%%%%%%%%%%%% in logfile
keeplooking = 0;                % ignores the occasional gap in the file but notices the end of the file
while keeplooking < 20;         % i.e more than 20 empty returns signals end of file
   a = fscanf(fid,'%s',1);      %this just looks for the next bit of info and reads as string.
   if strcmp(a,''); keeplooking = keeplooking+1; else keeplooking = 0; end; %counts number of times there is nothing returned
%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for responses and stimuli in logfile 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pulses (MRI Pulse, for EEG-MRI only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(a,'Pulse')
    Code = fscanf(fid,'%s',1);          % Next column with entry is the "code" column (100 for Pulse)
    Time = fscanf(fid,'%s',1);          % Next column with entry is "Time" column
    CodePulse = str2num(Code);
    if first_time_only == 0;
        if CodePulse == 100; 
            sess_start=str2num(Time);       % Works out when session started from first MRI Pulse
            first_time_only = 1;            % Makes sure only 
        end
    else
        if CodePulse == 100;    
            ev.pulse_timing = [ev.pulse_timing;[Block BlockType [str2num(Time)-sess_start]]];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(a,'Response') 
     Code = fscanf(fid,'%s',1);
     Time = fscanf(fid,'%s',1);
     CodeNum = str2num(Code);
%       if  CodeNum == 5;               % This is the mouse button press that usually starts the experiment
%       sess_start=str2num(Time);       % Works out when session started from first press of mouse button (CodeNum 5)
%       sess= sess + 1;
%       end  
     if exist('sess_start');
         if startend == 0;                  % Responses during blocks
             if response_on == 0;           % Response only recorded if response_on = 0 (prevents recording double hit of button)
               ev.response = [ev.response;[Block Exp_Block BlockType CodeNum is_target [str2num(Time)-sess_start]]];
               response_on = 1;             % means next response will not be counted unless response_on reset to 0 (prevents recording double hit of button)
             elseif response_on == 1;
             end  
         else                               % startend=2 (beginning of block); startend=1 (responses between blocks)    
             if response_on == 0;         
               ev.response_other = [ev.response_other;[Block BlockType CodeNum [str2num(Time)-sess_start]]];
               response_on = 1;             
             elseif response_on == 1;
             end
         end    
     end
 end;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
 if strcmp(a,'Picture')
    Code = fscanf(fid,'%s',1);
    Code2 = fscanf(fid,'%s',1);
    
    switch Code2                         % used to look at each instance of what is in this position in the logfile in turn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%           Start Block Screens          %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
     case 'Block'                           %%%% Rest Block
        TrialCode = fscanf(fid,'%s',1);
        TrialCode2 = fscanf(fid,'%s',1);
       
        Time = fscanf(fid,'%s',1);
       
        Block = Block + 1;
        BlockType = 1;
        startend = 0;                       % Start of Block = 0
        response_on = 1;                    % to prevent first button press being registered in % correct
         
        ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];
         
        
     case 'MR'                              %%%% Post MR Pulse 0-Back; Post MR Pulse 2-Back; Post MR Pulse 4-Back;
        Code3 = fscanf(fid,'%s',1);         
        Code4 = fscanf(fid,'%s',1);
        TrialCode = fscanf(fid,'%s',1);
        TrialCode2 = fscanf(fid,'%s',1);
        TrialCode3 = fscanf(fid,'%s',1);
        TrialCode4 = fscanf(fid,'%s',1);
        
        Time = fscanf(fid,'%s',1);
        %Time_Adjust = Time + 21660;    % Adjust as the info screen lasts 2.15s (2.166 in logfile)
        
        if strcmp(Code4,'0-Back') ==  1     %%%% Post MR Pulse 0-Back
            Block = Block + 1;
            Exp_Block = Exp_Block + 1;
            BlockType = 2;
            startend = 0;                       % Start of Block = 0
            response_on = 1;                    % to prevent first button press being registered in % correct
  
            ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];
       
        elseif strcmp(Code4,'2-Back') ==  1     %%%% Post MR Pulse 2-Back
            Block = Block + 1;
            Exp_Block = Exp_Block + 1;
            BlockType = 3;
            startend = 0;                       % Start of Block = 0
            response_on = 1;                    % to prevent first button press being registered in % correct
  
            ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];    
        
        elseif strcmp(Code4,'4-Back') ==  1     %%%% Post MR Pulse 0-Back THIS NEVER OCCURS IN LOGFILE!
            Block = Block + 1;
            Exp_Block = Exp_Block + 1;
            BlockType = 4;
            startend = 0;                       % Start of Block = 0
            response_on = 1;                    % to prevent first button press being registered in % correct
  
            ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];   
            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%           End Block Screens           %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
     case 'Rest'                             %%%% Intro Rest; EndBl Rest
        TrialCode = fscanf(fid,'%s',1);
        TrialCode2 = fscanf(fid,'%s',1);
            
        Time = fscanf(fid,'%s',1);
            
        if strcmp(TrialCode,'EndBl') ==  1  %%%% EndBl Rest
            startend = 1;                       % End of Block = 1
            response_on = 0;
        
            ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];
        
        elseif strcmp(TrialCode,'Intro') ==  1
            stim_other_type = 3;
            response_on = 1;                % means next response is not recorded
            ev.stimulus_other = [ev.stimulus_other;[Block BlockType stim_other_type [str2num(Time)-sess_start]]];
        
        end   
        
     case '0-Back'                                %%%% EndBl 0-Back
        TrialCode = fscanf(fid,'%s',1);
        TrialCode2 = fscanf(fid,'%s',1);
        
        Time = fscanf(fid,'%s',1);
        
        startend = 1;                       % End of Block = 1
        response_on = 0;
        
        ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];
       
        
     case '2-Back'                                %%%% EndBl 2-Back
        TrialCode = fscanf(fid,'%s',1);
        TrialCode2 = fscanf(fid,'%s',1);

        Time = fscanf(fid,'%s',1);
        
        startend = 1;                       % End of Block = 1
        response_on = 0;
        
        ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];
       
        
     case '4-Back'                                %%%% EndBl 4-Back 
        TrialCode = fscanf(fid,'%s',1);
        TrialCode2 = fscanf(fid,'%s',1);

        Time = fscanf(fid,'%s',1);
        
        startend = 1;                       % End of Block = 1
        response_on = 0;
        
        ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%   Before Experiment Stimuli (and before block for n-Back)   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     case 'Intro'                           %%%% nBack Intro Screen; Rest Intro Screen
        Code3 = fscanf(fid,'%s',1);
        TrialCode = fscanf(fid,'%s',1);     %%%% nBack or Rest
        TrialCode2 = fscanf(fid,'%s',1);    %%%% Intro
        TrialCode3 = fscanf(fid,'%s',1);    %%%% Screen

        if strcmp(TrialCode,'nBack') ==  1
            Time = fscanf(fid,'%s',1);     
            stim_other_type = 1;
            response_on = 1;                % means next response is not recorded
            ev.stimulus_other = [ev.stimulus_other;[Block BlockType stim_other_type [str2num(Time)-sess_start]]];
        elseif strcmp(TrialCode,'Rest') ==  1
            Time = fscanf(fid,'%s',1);     
            stim_other_type = 2;
            response_on = 1;                % means next response is not recorded
            ev.stimulus_other = [ev.stimulus_other;[Block BlockType stim_other_type [str2num(Time)-sess_start]]];
        end
        
     case '0'                               %%%% Intro 0 n-Back
        Code3 = fscanf(fid,'%s',1);
        TrialCode = fscanf(fid,'%s',1);     
        TrialCode2 = fscanf(fid,'%s',1);    
        TrialCode3 = fscanf(fid,'%s',1);  
        
        Time = fscanf(fid,'%s',1); 
        
        stim_other_type = 4;
        response_on = 1;                % means next response is not recorded
        ev.stimulus_other = [ev.stimulus_other;[Block BlockType stim_other_type [str2num(Time)-sess_start]]];
        
        %%%%%%% Code used for those logfiles without "Post MR Pulse n-Back"
        %Block = Block + 1;
        %BlockType = 2;
        %startend = 0;                       % Start of Block = 0
        %ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];
        
     case '2'                               %%%% Intro 2 n-Back
        Code3 = fscanf(fid,'%s',1);
        TrialCode = fscanf(fid,'%s',1);     
        TrialCode2 = fscanf(fid,'%s',1);    
        TrialCode3 = fscanf(fid,'%s',1);
        
        Time = fscanf(fid,'%s',1); 
        
        stim_other_type = 5;
        response_on = 1;                % means next response is not recorded
        ev.stimulus_other = [ev.stimulus_other;[Block BlockType stim_other_type [str2num(Time)-sess_start]]];
        
        %%%%%%% Code used for those logfiles without "Post MR Pulse n-Back"
        %Block = Block + 1;
        %BlockType = 3;
        %startend = 0;                       % Start of Block = 0
        %ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];
        
     case '4'                               %%%% Intro 4 n-Back
        Code3 = fscanf(fid,'%s',1);
        TrialCode = fscanf(fid,'%s',1);     
        TrialCode2 = fscanf(fid,'%s',1);    
        TrialCode3 = fscanf(fid,'%s',1);  
        
        Time = fscanf(fid,'%s',1); 
        
        stim_other_type = 6;
        response_on = 1;                % means next response is not recorded
        ev.stimulus_other = [ev.stimulus_other;[Block BlockType stim_other_type [str2num(Time)-sess_start]]];
        
        %%%%%%% Code used for those logfiles without "Post MR Pulse n-Back"
        %Block = Block + 1;
        %BlockType = 4;
        %startend = 0;                       % Start of Block = 0
        %ev.blockstartend = [ev.blockstartend;[Block BlockType startend [str2num(Time)-sess_start]]];
        
%%%%%%%%%% All Stimuli in the block (all the numbers 1-9)       
        
     otherwise                              %%%% All Stimuli in the Block       
        
        TrialCode = fscanf(fid,'%s',1);
        TrialCode2 = fscanf(fid,'%s',1);
        Number = fscanf(fid,'%s',1);
        Target = fscanf(fid,'%s',1);
        Time = fscanf(fid,'%s',1);
        NumberNum = str2num(Number);
        
        switch Target
       
        case 'Yes'                          %%%% Number is a Target
           is_target = 1;
           response_on = 0;         % Next response counted
           
           ev.stimulus = [ev.stimulus;[Block BlockType NumberNum is_target [str2num(Time)-sess_start]]];
          
     
        case 'No'                           %%%% Number is Non-Target
           is_target = 0;
           response_on = 0;         % Next response counted
           
           ev.stimulus = [ev.stimulus;[Block BlockType NumberNum is_target [str2num(Time)-sess_start]]];
       
        otherwise                           %%%% Is Not a Yes (Target) or No (Non-Target) Stimuli - Ignores it. 
       
        end;
	 end;
	 
 end;  % END Looking for Event Type "Picture"
     
 end;  % END "while keep looking<20", i.e. end looking at logfile                                    
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Looking through Logfile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%% AND NOW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modification of existing arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Create Block Start Duration array %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

     ev.blockstartduration=ev.blockstartend(find(ev.blockstartend(:,3) == 0),:); % find start of each block
     ev.blockstartduration=[ev.blockstartduration ev.blockstartend(find(ev.blockstartend(:,3) == 1),end)-ev.blockstartduration(:,end)]; % add duration at end of array
     %%%% [Block, Block_Type, StartEnd, Time, Duration]
     
     ev.blockstartduration(4,2) = 4;        % Change "block_type" to 4 (from 3) for 4-Back
     ev.blockstartduration(6,2) = 4;
     ev.blockstartduration(12,2) = 4; 
     ev.blockstartduration(13,2) = 4; 
     ev.blockstartduration(19,2) = 4; 
     ev.blockstartduration(22,2) = 4; 
     ev.blockstartduration(29,2) = 4; 
     ev.blockstartduration(33,2) = 4; 
     ev.blockstartduration(34,2) = 4; 
     ev.blockstartduration(38,2) = 4; 
         
     ev.blockstartduration_rest = ev.blockstartduration(find(ev.blockstartduration(:,2) == 1),:);
     ev.blockstartduration_task = ev.blockstartduration(find(ev.blockstartduration(:,2) == 2 | ev.blockstartduration(:,2) == 3 | ev.blockstartduration(:,2) == 4),:);
     
     ev.blockstartduration_task_adjusted = ev.blockstartduration_task;
     ev.blockstartduration_task_adjusted(:,4) = (ev.blockstartduration_task_adjusted(:,4)+ time_adjust_2);
     ev.blockstartduration_task_adjusted(:,5) = (ev.blockstartduration_task_adjusted(:,5)- time_adjust_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Create Percent Correct Array %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     
     if response_button == 1;
         if isempty(ev.response) == 1;                      %if no responses at all
         
         else 
            for ii = 1:size(ev.blockstartduration_task_adjusted, 1) 
                ev.singleblock = ev.response((find(ev.response(:,2) == ii)), :);    % find only block ii (block 1-Block)
                num_correct_T = size(find(ev.singleblock(:,4) == 1 & ev.singleblock(:,5) == 1)); % Num correct target: response = 2 and is_target = 1
                num_correct_N = size(find(ev.singleblock(:,4) == 2 & ev.singleblock(:,5) == 0)); % Num correct non-target: response = 1 and is_target = 0
                num_correct = num_correct_T(1,1) + num_correct_N(1,1);      % add these two together

                ev.percent_correct = [ev.percent_correct;[ii ev.singleblock(1,3) ((num_correct/trial_per_block)*100) ((num_correct_T(1,1)/target_stimuli_block)*100) ((num_correct_N(1,1)/nontarget_stimuli_block)*100)]];
                %%%% [ExpBlock, Block_Type, Total % Correct, Target % Correct, Non-Target % Correct]
            end;
         end
            
     elseif response_button == 2;  
         if isempty(ev.response) == 1;                      %if no responses at all
         
         else    
            for ii = 1:size(ev.blockstartduration_task_adjusted, 1)                               % For each of the Blocks         
                ev.singleblock = ev.response((find(ev.response(:,2) == ii)), :);    % find only block ii (block 1-Block)
                num_correct_T = size(find(ev.singleblock(:,4) == 2 & ev.singleblock(:,5) == 1)); % Num correct target: response = 2 and is_target = 1
                num_correct_N = size(find(ev.singleblock(:,4) == 1 & ev.singleblock(:,5) == 0)); % Num correct non-target: response = 1 and is_target = 0
                num_correct = num_correct_T(1,1) + num_correct_N(1,1);      % add these two together
            
                ev.percent_correct = [ev.percent_correct;[ii ev.singleblock(1,3) ((num_correct/trial_per_block)*100) ((num_correct_T(1,1)/target_stimuli_block)*100) ((num_correct_N(1,1)/nontarget_stimuli_block)*100)]];
                %%%% [ExpBlock, Block_Type, Total % Correct, Target % Correct, Non-Target % Correct]
            end;   
         end      
     end
     
     ev.percent_correct(3,2) = 4;       % Change "block_type" to 4 (from 3) for 4-Back
     ev.percent_correct(5,2) = 4;
     ev.percent_correct(10,2) = 4;
     ev.percent_correct(11,2) = 4;
     ev.percent_correct(15,2) = 4;
     ev.percent_correct(18,2) = 4;
     ev.percent_correct(22,2) = 4;
     ev.percent_correct(26,2) = 4;
     ev.percent_correct(27,2) = 4;
     ev.percent_correct(29,2) = 4; 