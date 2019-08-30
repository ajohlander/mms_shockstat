% SHOCKSTATMMS master file to repo
%
%   ----------------- OPERATIONS -----------------
%       
%   ----- Operations on SITL txt lists ------
%       - Calculate shock parameters: select a txt list and calculate various
%       shock parameters for the events. Relies on reading a large amount
%       of MMS data and is therefore slow. The function relies on a SITL
%       text list and a time interval list in csv format.The function 
%       typically writes a mat file containing saved data.
%       - Set time intervals for up- and downstream: Visually set up- and
%       downstream intervals for every event. Writes a csv file with time
%       intervals.
%       - Get line numbers of events with stable OMNI data: Automatic
%       selection of events with stable solar wind conditions. Writes new
%       SITL txt list. 
%       
%   ----- Operations on MAT lists ------
%       -Load shock parameters: loads mat file with certain constraints.
%       -Plot various shock parameters: plots a large amount of figures
%       doing various statistics.
%       -Construct bow shock model: Fits a conic section to all spacecraft
%       positions.
%
%   ----- Plots from SITL txt lists ------
%       All functions save figures as pngs
%       -MMS vs OMNI plots: comparison between MMS and OMNI data.
%       -Ion burst overview plots
%       -Electron burst overview plots
%       -Ion composition overview plots: uses HPCA data.
% 
%  ----- Other ------
%       - Write shock parameters to CSV file: takes a mat file and writes
%       certain parameters to a csv file
%
%
%   ----------------- FILES -----------------
%   The programs mainly deals with three different types of data files:
%
%       - SITL files in txt format: Contains lines copied from SITL
%       reports, preferrably where segments are merged into one line. Each
%       line constitutes one shock crossing.
%
%       - Time interval files in cvs format: Contains time intervals
%       specifying one upstrean and one downstream time interval for each
%       shock crossing. This file can only be constructed by visually
%       selecting up- and downstream (option in SHOCKSTATMMS).
%
%       - Files containing calculated shock parameters in mat format:
%       Calculated using the other two types of files. Data analysis is
%       best done on these files.




%% always do this
u = irf_units;

%% Select which program to call
optStr = sprintf(['\nPlease select an option: \n',...
    '1: Help \n',...
    '2: Calculate shock parameters \n',...
    '3: Set time intervals for up- and downstream \n',...
    '4: Get line numbers of events with stable OMNI data \n',...
    '5: Load shock parameters \n',...
    '6: Plot various shock parameters \n',...
    '7: Construct bow shock model \n',...
    '8: MMS vs OMNI plots \n',...
    '9: Ion burst overview plots \n',...
    '10: Electron burst overview plots\n',...
    '11: Ion composition burst overview plots (experimental) \n',...
    '12: Write shock parameters to CSV file']);

fprintf([optStr,'\n'])

scriptFileNameArr = {'',...
    'getShockParameters',...
    'setUpDownTimes',...
    'getStableEvents',...
    'loadShockParameters',...
    'plotShockParameters',...
    'constructBSmodel',...
    'mmsVomni',...
    'plotShockOverview',...
    'plotShockOverview',...
    'plotShockOverview',...
    'writeShockParametersToFile'};


optNum = irf_ask('Select option: [%]>','optNum',1);
fprintf('\n')

% display help if number 1 was selected
if optNum == 1
    help ShockStatMMS
    return
end

%% Ask if data should be loaded
% load data for plotShockParameters, loadShockParameters, constructBS,...
if optNum == 5 || optNum == 6 || optNum == 7 || optNum == 12 % must be loaded
    doLoadData = 1;
else
    doLoadData = 0;
end


%% If not writeshockparameters, select list file, sc number, start/stop line

if ~doLoadData && optNum ~= 12
    
    listFileArr = dir('*.txt');
    preSelectedFile = 1;
    for ii = 1:length(listFileArr)
        fprintf([num2str(ii),':     ',listFileArr(ii).name,'\n'])
        % make shock_list default file
        if strcmp(listFileArr(ii).name,'shock_list.txt')
            preSelectedFile = ii;
        end
    end
    
    listFileInd = irf_ask('Select file: [%]>','listFileInd',preSelectedFile);
    listFileName = listFileArr(listFileInd).name;
    
    if optNum ~= 6 && optNum ~= 7 && optNum ~= 9
        % sc number
        fprintf('\n')
        ic = irf_ask('Select spacecraft: [%]>','ic',2);
    end
    
    % which line in list to start at
    fprintf('\n')
    startLine = irf_ask('Start at line: [%]>','startLine',1);
    stopLine = irf_ask('Stop at line: [%]>','stopLine',20000);
    
end



%% Select mat file to load if requested

if doLoadData && optNum ~= 5
    if optNum == 6
        colorMode = irf_ask('Select color scheme for plots(1:Dark, 2:Light): [%]>','colorMode',2);
    end
    fprintf('\n')
    loadShockParameters
end


%% choose type of overview plot if requested

if optNum == 9
    plotType = 'ion';
elseif optNum == 10
    plotType = 'electron';
elseif optNum == 11
    plotType = 'ioncomp';
end


%% Call the requested script
fprintf('\n')
eval(scriptFileNameArr{optNum})
    
    
