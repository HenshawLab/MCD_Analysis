close all
clear all
clc

matlabrc;
addpath('peripherals');

%% Instructions:
% Keep in the subfolder "peripherals", as the subcodes are kept there
% Does not need to be with the main data directory, I have just included an
% example of it

% Folder Layout:
% Main Directory
    % Bio Replicates 
        % Technical Replicates

% Experiment layout
% Top of image - ASW
% Middle of image - Cells
% Bottom of image - Stimulant

% Set the variable "chem_location" to top or bottom at the start to account
% for this

%% Setup - change these values

% Inputs
MainDir = 'D:/WakeTest/'; % Important to end with the "/"
ExpName = 'ExperimentName'; % Name you want on the plot titles/folder name
WorkingDir = pwd; % Working Directory, where we are currently
imgextension = '*.tif'; % Image format (.tiff, .tif, .jpg etc.)

BioReps = [1]; % Bio replicates. For example, if you did 4 replicates but only want the first, second, and fourth, set this to be [1,2,4]
Reps = [1,2]; % Technical replicates, same style as BioReps
NBio = length(BioReps); NRep = length(Reps);

NLoop = 75; % Number of cycles each channel is imaged for
NPos = 6; % Number of channels (positions) used
concentrations = [200,49,3.6,0.12,0.0012,0];

imaging_period = 10; % Rough timestep (seconds) between measurements. Based off my previous guesses
chem_location = 'bottom'; % Choices: 'bottom' , 'top'. Location of chemotstimulant relative to the channel

% Choices for preprocessing/analysis - so you don't have to run the whole thing for a
% specific section. Default: true
RUN_BACKGROUND = true;
RUN_CROPPING = false;
RUN_PRETRACKPARAMETERS = false;

RUN_PARTICLELOCATION = true;
RUN_ANALYSIS = true;

% Set the particle to be 'bright' or 'dark' dependent on your imaging setup
% Default: dark
particle_type = 'dark';

% Image naming parameters
filenaming = 'custom'; % 'custom' for images as provided, 'default' if enough leading zeros to put files in correct order

% Width of accumulation region (from boundaries, in microns)
accum_width = 200;

BinW = 25; % Bin width (microns) for heatmap
Exclusion = 75; % Width (+/-) that is excluded from centre of heatmaps for PLOTTING PURPOSES ONLY

Mag = 6; % TOTAL Magnification
PixSize = 3.45; % Pixel size of your camera CHECK
PixToMum = PixSize/Mag; % Conversion of pixels to microns

%% Outputs - don't change these

OutputMain = 'MCD_Analysis/'; mkdir(OutputMain);
Output_Background = [OutputMain 'BackgroundImages/']; 
Output_Cropping = [OutputMain 'CroppingLimits/'];
Output_Pretrack = [OutputMain 'PretrackingParameters/'];
mkdir(Output_Background); mkdir(Output_Cropping); mkdir(Output_Pretrack);
OutputDir = [OutputMain ExpName '/']; mkdir(OutputDir);
PreProcDir = [OutputDir 'PreProcessing/']; % Place where preproc will be saved
mkdir(PreProcDir);
ExperimentOutDir = [OutputMain '/' ExpName '/'];
mkdir(ExperimentOutDir);

% Output directories
mkdir([ExperimentOutDir 'Bacteria_Positions/']);
mkdir([ExperimentOutDir 'AccumulationCurves/']);

FigDir = [ExperimentOutDir 'Figures/'];
PNGDir = [FigDir 'PNGS/']; 
mkdir(FigDir); mkdir(PNGDir);
mkdir([FigDir 'Heatmaps/']); mkdir([PNGDir 'Heatmaps']);
BetaFigDir = [OutputMain 'Beta/Figs/'];
BetaDir = [OutputMain 'Beta/'];
BetaPNGDir = [OutputMain 'Beta/Figs/PNGS/'];

mkdir(BetaDir); mkdir(BetaFigDir); mkdir(BetaPNGDir);

%% Run codes

% Wake_MCD_PreAnalysis;
% Wake_MCD_Analysis;
% Wake_MCD_BetaPlotting;
Wake_MCD_MakeSpreadsheets;
