
% The boolean arrays are 1 for accepted criteria, -1 for rejected criteria,
% and 0 for not reached


% number of events (should be more than actual events)
N = 1000;

% ask wether to write new text file
writeToTextFile = irf_ask('Write new text file with accepted lines= (0:no, 1:yes) [%]>','writeToTextFile',0);

if writeToTextFile
    fileName = irf_ask('File name: [%]>','fileName','shock_list_accepted');
end

% Ask for parameters
intMinutes = irf_ask('How many minutes of interval (+-)? [%]>','intMinutes',5);

% |B|
magAbsStdLim = irf_ask('Standard deviation limit for |B| (in nT) [%]>','magAbsStdLim',20);
magAbsStdBool = zeros(1,N);
magAbsLim = irf_ask('Maximum deviation limit for |B| (in nT) [%]>','magAbsLim',20);
magAbsBool = zeros(1,N);

% B direction
magDegStdLim = irf_ask('Standard deviation limit for B direction (in deg) [%]>','magDegStdLim',10);
magDegStdBool = zeros(1,N);
magDegLim = irf_ask('Maximum deviation limit for B direction (in deg) [%]>','magDegLim',30);
magDegBool = zeros(1,N);

% |V|
velAbsStdLim = irf_ask('Standard deviation limit for |V| (in km/s) [%]>','velAbsStdLim',200);
velAbsStdBool = zeros(1,N);
velAbsLim = irf_ask('Maximum deviation limit for |V| (in km/s) [%]>','velAbsLim',50);
velAbsBool = zeros(1,N);

% V direction
velDegStdLim = irf_ask('Standard deviation limit for V direction (in deg) [%]>','velDegStdLim',50);
velDegStdBool = zeros(1,N);
velDegLim = irf_ask('Maximum deviation limit for V direction (in deg) [%]>','velDegLim',20);
velDegBool = zeros(1,N);

% n
nStdLim = irf_ask('Standard deviation limit for n (in /cc) [%]>','nStdLim',20);
nStdBool = zeros(1,N);
nLim = irf_ask('Maximum deviation limit for n (in /cc) [%]>','nLim',50);
nBool = zeros(1,N);

% Ma
MaStdLim = irf_ask('Standard deviation limit for SW Ma [%]>','MaStdLim',3);
MaStdBool = zeros(1,N);
MaLim = irf_ask('Maximum deviation limit for SW Ma [%]>','MaLim',6);
MaBool = zeros(1,N);

% overall
lineBool = zeros(1,N);


%% read event list
tline = 1;
lineNum = 0;

fid = fopen(listFileName);

% loop through skipped lines
for ii = 1:startLine-1
    lineNum = lineNum+1;
    tline = fgets(fid);
end



%% create and open new text file if requested

if writeToTextFile
    % create text file
    irf.log('w','Creating new text file. Make sure relevant data is not overwritten.')
    system(['touch ',fileName,'.txt'])
    
    % open text file
    fidWrite = fopen([fileName,'.txt'],'w');
end

%% super trooper looper
accLineNums = zeros(1,1e3); % a "large" array with line numbers
nLine = 1; % count accepted line numbers
while tline ~= -1
    lineNum = lineNum+1;
    disp(['Current line number: ',num2str(lineNum)])
    
    %% end loop if requested
    if lineNum >= stopLine; disp('Reached stop line, exiting...'); break; end
    
    %% read line from file
    tline = fgets(fid);
    
    if tline(1) == -1 || ~strcmp(tline(1),'2')
        disp('done!')
        break;
    end
    tintStr = [tline(1:10),'T',tline(12:19),'/',tline(23:32),'T',tline(34:41)];
    
    % all other strings if needed
    restStr = tline(44:end);
    
    % find commas
    idcomma = strfind(restStr,',');
    % find start paranthesis
    idpar = strfind(restStr,'(');
    
    fomStr = restStr(1:idcomma(1)-1);
    sitlStr = restStr(idcomma(1)+2:idpar(1)-1);
    
    %% get time interval
    tint = irf.tint(tintStr);
    
    % add some time
    tint2 = tint(1)+intMinutes*60*[-1,1]+diff(tint.epochUnix)/2;
    
    %% read omni data
    
    ff= irf_get_data(tint2,'bx,by,bz,vx,vy,vz,n,Ma','omni_min');
    
    if isempty(ff)
        irf.log('c','failed to read omni data, skipping...')
        continue;
    else
        omniTime = irf_time(ff(:,1),'epoch>epochTT');
        Bomni = irf.ts_vec_xyz(omniTime,ff(:,2:4));
        % correct velocity for abberation
        Vomni = irf.ts_vec_xyz(omniTime,ff(:,5:7)+[0,29.8,0]);
        Nomni = irf.ts_scalar(omniTime,ff(:,8));
        Maomni = irf.ts_scalar(omniTime,ff(:,9));
    end
    
    %for good measures, remove old ff
    clear ff
    
    
    %% get average upstream omni values
    
    Bu = nanmean(Bomni.tlim(tint).data,1);
    if isnan(Bu(1))
        disp('No magnetic field data, skipping,...')
        continue;
    end
    Vu = nanmean(Vomni.tlim(tint).data,1);
    nu = nanmean(Nomni.tlim(tint).data,1);
    Mau = nanmean(Maomni.tlim(tint).data,1);
    if isnan(Vu(1)) || isnan(nu)
        disp('No plasma data, skipping,...')
        continue;
    end
    
    
    
    
    %% Get differences and deviations
    % calculate angles from "upstream" values
    
    thB = Bomni.abs;
    thB.data = acosd((Bomni.data*Bu')./(Bomni.abs.data*norm(Bu)));
    
    thV = Vomni.abs;
    thV.data = acosd((Vomni.data*Vu')./(Vomni.abs.data*norm(Vu)));
    
    % angle values deviation
    thBstd = nanstd(thB.data);
    thVstd = nanstd(thV.data);
    thBdiff = max(thB.data);
    thVdiff = max(thV.data);
    
    % scalar values deviation
    absBstd = nanstd(abs(Bomni.abs.data-norm(Bu)));
    absVstd = nanstd(abs(Vomni.abs.data-norm(Vu)));
    Nstd = nanstd(abs(Nomni.abs.data-norm(nu)));
    Mastd = nanstd(abs(Maomni.abs.data-norm(Mau)));
    
    absBdiff = max(abs(Bomni.abs.data-norm(Bu)));
    absVdiff = max(abs(Vomni.abs.data-norm(Vu)));
    Ndiff = max(abs(Nomni.abs.data-norm(nu)));
    Madiff = max(abs(Maomni.abs.data-norm(Mau)));
    
    
    %% Check if criterias are fulfilled
    
    isOK = 1;
    
    fprintf('\n')
    disp('Checking OMNI data to limits specified...')
    
    % B
    fprintf(sprintf('B magnitude standard deviation: '))
    if absBstd<magAbsStdLim
        fprintf(sprintf('OK \n'))
        magAbsStdBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        magAbsStdBool(lineNum) = -1;
    end
    
    fprintf(sprintf('B magnitude maximum deviation: '))
    if absBdiff<magAbsLim
        fprintf(sprintf('OK \n'))
        magAbsBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        magAbsBool(lineNum) = -1;
    end
    
    
    fprintf(sprintf('B angle standard deviation: '))
    if thBstd<magDegStdLim
        fprintf(sprintf('OK \n'))
        magDegStdBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        magDegStdBool(lineNum) = -1;
    end
    
    fprintf(sprintf('B angle maximum deviation: '))
    if thBdiff<magDegLim
        fprintf(sprintf('OK \n'))
        magDegBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        magDegBool(lineNum) = -1;
    end
    
    % V
    fprintf(sprintf('V magnitude standard deviation: '))
    if absVstd<velAbsStdLim
        fprintf(sprintf('OK \n'))
        velAbsStdBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        velAbsStdBool(lineNum) = -1;
    end
    
    fprintf(sprintf('V magnitude maximum deviation: '))
    if absVdiff<velAbsLim
        fprintf(sprintf('OK \n'))
        velAbsBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        velAbsBool(lineNum) = -1;
    end
    
    
    fprintf(sprintf('V angle standard deviation: '))
    if thVstd<velDegStdLim
        fprintf(sprintf('OK \n'))
        velDegStdBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        velDegStdBool(lineNum) = -1;
    end
    
    fprintf(sprintf('V angle maximum deviation: '))
    if thVdiff<velDegLim
        fprintf(sprintf('OK \n'))
        velDegBool(lineNum) = 1;
    else
        isOK = 0;
        velDegBool(lineNum) = -1;
    end
    
    
    % N
    fprintf(sprintf('N standard deviation: '))
    if Nstd<nStdLim
        fprintf(sprintf('OK \n'))
        nStdBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        nStdBool(lineNum) = -1;
    end
    
    fprintf(sprintf('N maximum deviation: '))
    if Ndiff<nLim
        fprintf(sprintf('OK \n'))
        nBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        nBool(lineNum) = -1;
    end
    
    
    % Ma
    fprintf(sprintf('Ma standard deviation: '))
    if Mastd<MaStdLim
        fprintf(sprintf('OK \n'))
        MaStdBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        MaStdBool(lineNum) = -1;
    end
    
    fprintf(sprintf('Ma maximum deviation: '))
    if Madiff<MaLim
        fprintf(sprintf('OK \n'))
        MaBool(lineNum) = 1;
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
        MaBool(lineNum) = -1;
    end
    
    fprintf('\n')
    
    if isOK
        accLineNums(nLine) = lineNum;
        nLine = nLine+1;
    end
    
    lineBool(lineNum) = isOK*2-1;
    
    
    %% Write to text file if requested
    
    if writeToTextFile && isOK
        fprintf(fidWrite,'%50s',tline);
    end
    
end

if writeToTextFile
    fclose(fidWrite);
    disp('File written!')
end


%% Print results (TODO)
fprintf('\n')

disp('Criteria')
% should also print all selection criterias
fprintf(['|B| diff std limit: \t\t ',num2str(magAbsStdLim),' nT \n'])
fprintf(['|B| diff max limit: \t\t ',num2str(magAbsLim),' nT \n\n'])

fprintf(['B angle diff std limit: \t ',num2str(magDegStdLim),' deg \n'])
fprintf(['B angle diff max limit: \t ',num2str(magAbsLim),' deg \n\n'])

fprintf(['|V| diff std limit: \t\t ',num2str(velAbsStdLim),' km/s \n'])
fprintf(['|V| diff max limit: \t\t ',num2str(velAbsLim),' km/s \n\n'])

fprintf(['V angle diff std limit: \t ',num2str(velDegStdLim),' deg \n'])
fprintf(['V angle diff max limit: \t ',num2str(velAbsLim),' deg \n\n'])

fprintf(['N diff std limit: \t\t ',num2str(nStdLim),' cm^-3 \n'])
fprintf(['N diff max limit: \t\t ',num2str(nLim),' cm^-3 \n\n'])

fprintf(['Ma diff std limit: \t\t ',num2str(MaStdLim),' \n'])
fprintf(['Ma diff max limit: \t\t ',num2str(MaLim),' \n\n'])

fprintf('\n')
disp('Resulting line numbers:')
disp(num2str(accLineNums(accLineNums~=0)'))


%% Make bar graph of rejected/accepted lines

rejectFraction = [numel(find(magAbsStdBool==-1)),numel(find(magAbsBool==-1)),...
    numel(find(magDegStdBool==-1)),numel(find(magDegBool==-1)),...
    numel(find(magAbsStdBool==-1 | magAbsBool==-1 | magDegStdBool==-1 | magDegBool==-1)),... % combined magnetic field
    numel(find(velAbsStdBool==-1)),numel(find(velAbsBool==-1)),...
    numel(find(velDegStdBool==-1)),numel(find(velDegBool==-1)),...
    numel(find(velAbsStdBool==-1 | velAbsBool==-1 | velDegStdBool==-1 | velDegBool==-1)),... % combined magnetic field
    numel(find(nStdBool==-1)),numel(find(nBool==-1)),...
    numel(find(nStdBool==-1 | nBool==-1)),... % combined density
    numel(find(MaStdBool==-1)),numel(find(MaBool==-1)),...
    numel(find(MaStdBool==-1 | MaBool==-1)),... % combined Ma
    0,numel(find(lineBool==-1))]/numel(find(magAbsStdBool));

%criteriaNames = {'B std','B','B deg std','B deg','B tot','V std','V','V deg std','V deg','V tot','n std','n'};
criteriaNames = {'$\delta B$','$\Delta B$',...
    '$\delta\theta_B$','$\Delta\theta_B$',...
    '$B$ Total',...
    '$\delta V$','$\Delta V$',...
    '$\delta\theta_V$','$\Delta\theta_V$',...
    '$V$ Total',...
    '$\delta n$','$\Delta n$',...
    '$n$ Total',...
    '$\delta M_A$','$\Delta M_A$',...
    '$M_A$ Total',...
    '','Total'};

fig = figure;
hca = axes(fig);

hbar = bar(hca,rejectFraction*100);

hca.XTick = 1:length(rejectFraction);
hca.XTickLabel = criteriaNames;
hca.XTickLabelRotation = 45;

hbar.FaceColor = [211,64,82]/255;
hbar.LineWidth = 1.3;
hca.XTickLabelRotation = 45;

hca.TickLabelInterpreter = 'latex';


hca.Position(2) = 0.2;
hca.Position(4) = 0.75;
hca.FontSize = 14;
hca.LineWidth = 1.3;

hca.YLim = [0,100];


%% smaller version of the bar plot just for Bangle and Ma


rejectFraction2 = [numel(find(magDegStdBool==-1)),numel(find(magDegBool==-1)),...
    numel(find(MaStdBool==-1)),numel(find(MaBool==-1)),...
    numel(find(lineBool==-1))]/numel(find(magAbsStdBool));

criteriaNames2 = {'$\delta\theta_B$','$\Delta\theta_B$',...
    '$\delta M_A$','$\Delta M_A$',...
    'Total'};

fig = figure;
hca = axes(fig);

hbar = bar(hca,rejectFraction2*100);

hca.XTickLabel = criteriaNames2;
hca.XTickLabelRotation = 45;

hbar.FaceColor = [211,64,82]/255;
hbar.LineWidth = 1.3;
hca.XTickLabelRotation = 45;

hca.TickLabelInterpreter = 'latex';


hca.Position(2) = 0.2;
hca.Position(4) = 0.75;
hca.FontSize = 14;
hca.LineWidth = 1.3;

hca.YLim = [0,100];

