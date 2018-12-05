# SPM fMRI Example Pipeline
Example pipeline for preprocessing, reading in behavioural files, first level and second level analysis

This folder contains:

	Powerpoint summary of pipeline used (Example_SPM_Pipeline)
	Matlab Scripts
	SPM Batch Files

The scripts are organised such that everything runs from the "_multisubject_analysis_" script

This script calls all the other scripts, running through subjects and folders to do this. 

In order:

	First do the preprocessing ("_preprocess_" script) 
	After this is completed, you can run first level analysis ("_firstlevel_analysis_")
	This uses the script "_behavioural_data_onset_duration_" to load in onsets and durations of stimuli 
	After this second level analysis ("_secondlevel_analysis_") can be done.

The scripts will not run "off the shelf" perfectly, and will need some modification for your data pathways and design (e.g. change file paths). This will definitely be the case for the script to read behavioural data ("_behavioural_data_onset_duration_"), but will also be likely for the first and second level analyses. 

It is also only an example pipeline (although based on a simplification of one used in published analyses). Therefore, there may be a better design or paradigm that fits your specific data. However, it is hopefully a good place to start when learning SPM analysis, to run preliminary analysis on data, and to adapt to fit your data as you become comfortable with it. 

There are comments within the scripts to help understand the script and what it is doing, and where to modify it. And where possible, the modifiable variables are usually at the top of the script. 

You may also notice that the script basically opens up an SPM batch file (kept in the "_batch_files_" folder) and fills it in with the data and variables for the study, then runs it. This can be hard coded in SPM (ie no need for batch files, just directly run the analysis), but filling in batch files from script gives a balance between coding and direct interface input that might be more intuitive when learning SPM analysis. 

You can load the batch files and see what they look like in SPM interface, then load them in the script and see what variables are contained within to have an idea of what the program is doing. 
 
	
