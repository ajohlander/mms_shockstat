% Master file for repo, only call this script
% Asks user what to do
%

% always do this
u = irf_units;

%% Select which program to call
optStr = sprintf(['\nPlease select an option: \n',...
    '1: TBD \n',...
    '2: Ion burst overview plots \n',...
    '3: Electron burst overview plots \n',...
    '4: Ion composition burst overview plots \n',...
    '5: Calculate shock parameters \n',...
    '6: Plot various shock parameters \n',...
    '7: Construct bow shock model \n',...
    '8: MMS vs OMNI plots \n',...
    '9: Get line numbers of events with stable OMNI data \n',...
    '10: Set time intervals for up- and downstream']);

fprintf([optStr,'\n'])

scriptFileNameArr = {'',...
    'plotShockOverview',...
    'plotShockOverview',...
    'plotShockOverview',...
    'getShockParameters',...
    'plotShockParameters',...
    'constructBSmodel',...
    'mmsVomni',...
    'getStableEvents',...
    'setUpDownTimes'};


optNum = irf_ask('Select option: [%]>','optNum',1);
fprintf('\n')

%% Ask if data should be loaded
% only valid for shockParams and BSmodel
if optNum == 6 || optNum == 7 % must be loaded
    doLoadData = 1;
else
    doLoadData = 0;
end


%% Select list file, sc number, start/stop line

if ~doLoadData
    
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

if doLoadData
    loadShockParameters
end


%% choose type of overview plot if requested

if optNum == 2
    plotType = 'ion';
elseif optNum == 3
    plotType = 'electron';
elseif optNum == 4
    plotType = 'ioncomp';
end


%% Call the requested script
fprintf('\n')
eval(scriptFileNameArr{optNum})
    
    
