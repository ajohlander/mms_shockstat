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
    '9: Get line numbers of events with stable OMNI data']);

fprintf([optStr,'\n'])

scriptFileNameArr = {'',...
    'plotShockOverview',...
    'plotShockOverview',...
    'plotShockOverview',...
    'getShockParameters',...
    'plotShockParameters',...
    'constructBSmodel',...
    'mmsVomni',...
    'getStableEvents'};


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
    
    % Select data file
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
    
    %% colors for plots
    colorMode = irf_ask('Color mode (1:Fancy, 2:Boring): [%]>','colorMode',2);
    
    switch colorMode
        case 1 % darker, cooler
            cmap = 'strangeways';
            axcol = [1,1,1]*.4;
            figcol = [1,1,1]*.2;
            textcol = [1,1,1]*.95;
            daycol = [1,1,1]*.95;
            col1 = [253,232,159]/255;
            col2 = [211,64,82]/255;
            col3 = [53,151,103]/255;
            bincol = [157,214,166]/255;
            daysidecol = textcol;
        case 2 % lighter, boring
            cmap = 'strangeways';
            axcol = [1,1,1];
            figcol = [1,1,1];
            textcol = [1,1,1]*.05;
            daycol = [1,1,1];
            col1 = [1,1,1]*.5;
            col2 = [211,64,82]/255;
            col3 = [53,151,103]/255;
            % bincol = [157,214,166]/255;
            bincol = [100,100,100]/255;
            daysidecol = [1,1,1];
    end
    
    colmat = [col1;col2;col3];
    
    
    
    
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
    
    
