%% MAIN MISSION SCRIPT
%  Author: Group 2542
%  Date:   1/7/2026
%
%  DESCRIPTION:
%  This script initializes the environment, adds necessary folders to the 
%  path, and executes the complete simulation pipeline sequentially:
%  - Constants & Data Loading
%  - Phase Analysis & Grid Refinement
%  - 3D Grid Search
%  - Gravity Assist Optimization
%  - Trajectory Reconstruction & Visualization

%% 0. INITIALIZATION & PATH MANAGEMENT
clear; clc; close all;

fprintf('--------------------------------------------------\n');
fprintf('  INITIALIZING MISSION SIMULATION\n');
fprintf('--------------------------------------------------\n');

restoredefaultpath; 
addpath(genpath(pwd)); 
fprintf('> Path configured correctly.\n');

%% 1. CONFIG PARAMETERS
fprintf('> Loading all the parameters...\n');
[constants,data,constraints,Leg1,Leg2] = config();

%% 2a. EXECUTE PHASE 1a: GRID REFINEMENT
fprintf('> Phase Analysis to refine the grid\n');
ass1_resonanceTimeAnalysis

%% 3. EXECUTE PHASE 2: GRID SEARCH (LEG 1)
fprintf('> Running Grid Search...\n');
ass1_3DgridSearch 
% !!! to change time window for 3D grid search: 
%           ---> open "ass1_3DgridSearch.m" and uncomment the required window

%% 4. EXECUTE PHASE 3: OPTIMIZATION
fprintf('> Optimizing the solution...\n');
ass1_refinev2

%% 5. VISUALIZATION & OUTPUT
fprintf('> Generating Plots and Reports...\n');
