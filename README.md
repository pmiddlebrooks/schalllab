# schalllab
all schall lab related matlab analysis source code (replaces matlab/ folder)

## Folders:
### behavior
General functions related to behavioral performance. Specific task functions call these functions.

### ccm (choice countermanding)
Behavioral and neural analyses functions to analyze performance and brain activity during choice countermanding task.

### cmd (countermanding)
Behavioral and neural analyses functions to analyze performance and brain activity during simple (one target per trial) countermanding task.

### cmtb
Cognitive modeling toolbox. Just one file to calculate distribution function.

### create_datafile_pgm
Files written to translate data from Plexon to matlab readable format. 

### dPCA
A demixed principle component analysis code:
Wieland Brendel & Christian Machens, published at NIPS 2011 "Demixed
Principal Component Analysis", code@http://sourceforge.net/projects/dpca/

### figure_setup
Functions related to creating figures

### gng
Behavioral and neural analysis functions to analyze performance and brain activity during go/no-go task.

### local_data
Local repository for data storage. Gets called by analyses functions.

### maskbet
Metacognition task functions

### matlab_cod_bbz
Code created by Bram Zandbelt for specialized purposes.

### matlab_file_exchange_tools
Code downloaded from matlab file exchange.

### mem
Behavioral and neural analyses functions to analyze performance and brain activity during memory-guided saccade task.

### metacognition
Code to run a metacognitive task on the Eyelink human eye-tracking set-up

### modeling
Functions specific to cognitive modeling

### neural
General functions related to neuronal analyses. Specific task functions call these functions.

### sam
Stochastic accumulator model toolbox. Created by Bram Zandbelt. Models stochastic accumulators for countermanding tasks.

### toolboxes
Various functions and/or toolboxes downloaded from external sources, called by other functions.

### trial_visualization
Function(s) to visualize aspects of individual trials



## Functions/scripts in root:
#### batch.m
Scripts (junk) to run batches of other functions

#### cell_to_mat.m
Converts translated files' cell arrays to doubles if needed. Gets called by load_data.m

#### data_file_path.m
Determines where date files are

#### eeg_electrode_map.m
Maps monkey eeg electrode placements a la human coordinates

#### find_eeg_sessions.m
Finds data files with eeg data

#### find_sessions_data.m
Finds sessions with a particular data type (lfp, spike, eeg, etc)

#### find_task_sessions.m
Finds sessions associate with a particular task

#### get_environment.m
Determines whether you're working on your local computer or on ACCRE

#### load_data.m
Loads a single session dataset

#### local_data_path.m
Defines root folder to look for or save files on local computer

#### local_figure_path.m
Defines root folder to save figures on local computer

#### psychometric_function.m
Unknown. Variant of Weibull

#### scratch.m
Junk scripts to run sets of functions

#### subject_data_path.m
Defines where to look for data files locally or on Teba

#### task_session_array
Returns a set of sessions depending on what set you ask for.


