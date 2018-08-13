

% ask wether to write new text file
writeToTextFile = irf_ask('Write new text file with accepted lines= (0:no, 1:yes) [%]>','writeToTextFile',0);

if writeToTextFile
    fileName = irf_ask('File name: [%]>','fileName','shock_list_accepted');
end

% Ask for parameters
intMinutes = irf_ask('How many minutes of interval? [%]>','intMinutes',15);

magAbsStdLim = irf_ask('Standard deviation limit for |B| (in nT) [%]>','magAbsStdLim',2);
magDegStdLim = irf_ask('Standard deviation limit for B direction (in deg) [%]>','magDegStdLim',10);
magAbsLim = irf_ask('Maximum deviation limit for |B| (in nT) [%]>','magAbsLim',2);
magDegLim = irf_ask('Maximum deviation limit for B direction (in deg) [%]>','magDegLim',45);

velAbsStdLim = irf_ask('Standard deviation limit for |V| (in km/s) [%]>','velAbsStdLim',2);
velDegStdLim = irf_ask('Standard deviation limit for V direction (in deg) [%]>','velDegStdLim',5);
velAbsLim = irf_ask('Maximum deviation limit for |V| (in km/s) [%]>','magAbsLim',50);
velDegLim = irf_ask('Maximum deviation limit for V direction (in deg) [%]>','magDegLim',20);

nStdLim = irf_ask('Standard deviation limit for n (in /cc) [%]>','nStdLim',2);
nLim = irf_ask('Maximum deviation limit for n (in /cc) [%]>','nLim',5);


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
    
    ff= irf_get_data(tint2,'bx,by,bz,vx,vy,vz,n','omni_min');
    
    if isempty(ff)
        irf.log('c','failed to read omni data, skipping...')
        continue;
    else
        omniTime = irf_time(ff(:,1),'epoch>epochTT');
        Bomni = irf.ts_vec_xyz(omniTime,ff(:,2:4));
        % correct velocity for abberation
        Vomni = irf.ts_vec_xyz(omniTime,ff(:,5:7)+[0,29.8,0]);
        Nomni = irf.ts_scalar(omniTime,ff(:,8));
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
    absVstd = nanstd(abs(Bomni.abs.data-norm(Bu)));
    Nstd = nanstd(abs(Nomni.abs.data-norm(nu)));
    
    absBdiff = max(abs(Bomni.abs.data-norm(Bu)));
    absVdiff = max(abs(Bomni.abs.data-norm(Bu)));
    Ndiff = max(abs(Nomni.abs.data-norm(nu)));
    
    
    %% Check if criterias are fulfilled
    
    isOK = 1;
    
    fprintf('\n')
    disp('Checking OMNI data to limits specified...')
    
    % B
    fprintf(sprintf('B magnitude standard deviation: '))
    if absBstd<magAbsStdLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    fprintf(sprintf('B magnitude maximum deviation: '))
    if absBdiff<magAbsLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    
    fprintf(sprintf('B angle standard deviation: '))
    if thBstd<magDegStdLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    fprintf(sprintf('B angle maximum deviation: '))
    if thBdiff<magDegLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    % V
    fprintf(sprintf('V magnitude standard deviation: '))
    if absVstd<velAbsStdLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    fprintf(sprintf('V magnitude maximum deviation: '))
    if absVdiff<velAbsLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    
    fprintf(sprintf('V angle standard deviation: '))
    if thVstd<velDegStdLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    fprintf(sprintf('V angle maximum deviation: '))
    if thVdiff<velDegLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    
    % N
    fprintf(sprintf('N standard deviation: '))
    if Nstd<nStdLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    fprintf(sprintf('N maximum deviation: '))
    if absVdiff<nLim
        fprintf(sprintf('OK \n'))
    else
        isOK = 0;
        fprintf(sprintf('NOT OK \n'))
    end
    
    fprintf('\n')
    
    if isOK
        accLineNums(nLine) = lineNum;
        nLine = nLine+1;
    end
    
    
    %% Write to text file if requested
    
    if writeToTextFile && isOK
        fprintf(fidWrite,'%50s',tline);
    end
    
end

if writeToTextFile
    fclose(fidWrite);
end


%% Print results (TODO)

% should also print all selection criterias

fprintf('\n')
disp(num2str(accLineNums(accLineNums~=0)'))



