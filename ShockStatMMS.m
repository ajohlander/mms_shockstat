% Master file for repo, only call this script
% Asks user what to do
%

%% Select which program to call
optStr = sprintf(['\nPlease select an option: \n',...
    '1: TBD \n',...
    '2: Ion burst overview plots \n',...
    '3: Electron burst overview plots \n',...
    '4: Ion composition burst overview plots \n',...
    '5: Plot/calculate shock parameters \n',...
    '6: Construct bow shock model \n',...
    '7: MMS vs OMNI plots \n',...
    '8: Get line numbers of events with stable OMNI data']);

fprintf([optStr,'\n'])

scriptFileNameArr = {'',...
    'plotShockOverview',...
    'plotShockOverview',...
    'plotShockOverview',...
    'shockParameters',...
    'constructBSmodel',...
    'mmsVomni',...
    'getStableEvents'};


optNum = irf_ask('Select option: [%]>','optNum',1);
fprintf('\n')

%% Ask if data should be loaded
% only valid for shockParams and BSmodel
if optNum == 5
    doLoadData = irf_ask('Load data? (0:no, 1:yes) [%]>','doLoadData',0);
elseif optNum == 6 % must be loaded
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
    
    if optNum ~= 8
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
    
    dataFileArr = dir('*.mat');
    
    if isempty(dataFileArr)
        irf.log('c','No data files')
        
    else
        preSelectedFile = 1;
        for ii = 1:length(dataFileArr)
            fprintf([num2str(ii),':     ',dataFileArr(ii).name,'\n'])
            % make shockParams default file
            if strcmp(dataFileArr(ii).name,'shockParams.mat')
                preSelectedFile = ii;
            end
        end
        fprintf('\n')
        dataFileInd = irf_ask('Select file: [%]>','listFileName',preSelectedFile);
        
        load(dataFileArr(dataFileInd).name)
    end
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
    
    
