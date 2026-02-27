# Interplanetary Mission Analysis Project
**Jupiter to Asteroid via Mars Gravity Assist**

## Project Overview
This project performs a comprehensive trajectory analysis for an 
interplanetary mission starting from Jupiter, performing a Gravity Assist 
(Flyby) at Mars, and arriving at a target Asteroid. 
The simulation includes synodic analysis, Lambert arc solutions, 
porkchop plotting, and trajectory optimization.

## How to Run
The project is designed to be executed from a single entry point that 
manages the environment and path configurations automatically.

1. Open MATLAB.
2. Open the file **MAIN.m**.
3. Click the **Run** button.

The script will automatically add all necessary subfolders to the 
MATLAB path, load the constants, and execute the analysis phases 
sequentially.

## Configuration & Advanced Options
### Changing the 3D Grid Search Window
By default, the 3D optimization phase scans a specific time window. 
If you wish to explore different departure/arrival windows:

1. Open the script **ass1_3DGridSearch.m**.
2. Locate the "Time Window Settings" section (at the top of the file).
3. **Comment/Uncomment** the lines corresponding to the desired date range.
4. Save the file and re-run MAIN.m (or run ass1_3DGridSearch.m individually if checking that specific phase).

## Requirements
* MATLAB R2024b or later.
* Optimization Toolbox (recommended for fmincon/fzero).

---
*Author: Caccia Riccardo*