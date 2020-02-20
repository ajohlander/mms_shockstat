%% Set filename
fprintf('Writing new cvs file with shock times and some parameters \n')
csvFileName = irf_ask('Name of new csv file: [%]','csvFileName','mms_shock_burst');

if ~strcmp(csvFileName(end-3:end),'.csv')
    csvFileName = [csvFileName,'.csv'];
end

if isfile(csvFileName)
    overwriteCSV = irf_ask('File already exists, overwrite? (0:no, 1:yes) [%]:','overwriteCSV',0);
    if ~overwriteCSV; disp('exiting...'); return; end
end


%% Prepare parameters

N = length(TV);

yearV = zeros(N,1);
monthV = zeros(N,1);
dayV = zeros(N,1);

hourV = zeros(N,1);
minV = zeros(N,1);
secV = zeros(N,1);


for ii = 1:N
    tt = irf_time(TV(ii),'epoch>vector');
    
    yearV(ii) = tt(1);
    monthV(ii) = tt(2);
    dayV(ii) = tt(3);
    hourV(ii) = tt(4);
    minV(ii) = tt(5);
    secV(ii) = tt(6);
end


%%

writetable(table(...
    yearV,monthV,dayV,hourV,minV,secV,dTV,round(RV(:,1),1),...
    round(RV(:,2),1),round(RV(:,3),1),...
    round(thBnV1,0),round(MaV1,0),...
    'VariableNames',{'year' 'month' 'day','hour','minute','second','duration',...
    'Rx','Ry','Rz','theta_Bn','AlfvenMach'}),csvFileName)



