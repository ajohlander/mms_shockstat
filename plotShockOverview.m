% Master file for all overview plots

%% set stuff depending on plotType

guessFilePath = ['../mms_shock_plots/',plotType,'/'];

guessFilePrefix = ['sh_',plotType,'_'];

%% read time intervals

rePlotAll = irf_ask('Redo all plots? (0:no, 1:yes) [%]>','rePlotAll',0);

filePath = irf_ask('Choose path for plots [%]>','filePath',guessFilePath);

filePrefix = irf_ask('Prefix of file name [%]>','filePrefix',guessFilePrefix);

if strcmp(plotType,'ioncomp')
    hpcaMode = irf_ask('HPCA mode (auto/srvy/brst) [%]>','hpcaMode','auto');
end

u = irf_units;

% normal bs model to use (cannot be slho)
shModel = 'farris';

% for reduced ion distributions
nMC = 2e2;
vlim = 800;
vg = linspace(-vlim,vlim,100);

% plot position parameters
axl = 0.125; % axis left
axu = 0.06;
axw = .71;
axh = .87;
cbl = axl+axw+.01;
cbw = .03;


%% Make plot
tline = 1;
lineNum = 0;

fid = fopen(listFileName);

% loop through skipped lines
for ii = 1:startLine-1
    lineNum = lineNum+1;
    tline = fgets(fid);
end

% super trooper looper
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
    
    %%
    descStr = restStr(idcomma(2)+2:end);
    if length(descStr)>50
        descStr = descStr(1:50);
    end
    
    %% get time interval and check if already done
    tint = irf.tint(tintStr);
    
    tstr = tint(1).toUtc;
    tstr([5,8,14,17]) = '';
    tstr = tstr(1:15);
    
    listing = dir(filePath);
    fnCell = {listing.name};
    fileName = [filePrefix,tstr,'.png'];
    
    if ismember(fileName,fnCell) && ~rePlotAll
        disp(['Skipping ',tstr])
        continue;
    end
    
    
    
    %% get data
    
    % ----- common data is B, ni, ne, R -----
    c_eval('B? = mms.get_data(''B_gse_fgm_brst_l2'',tint,?);')
    c_eval('B = B?;',ic)
    
    if isempty(B)
        disp(['no FGM data for ',tstr,', skipping...'])
        continue;
    end
    ni = mms.get_data('Ni_fpi_brst_l2',tint,ic);
    ne = mms.get_data('Ne_fpi_brst_l2',tint,ic);
    
    try
        R = mms.get_data('R_gse',tint);
    catch
        disp(['error reading sc position for ',tstr,', skipping...'])
        continue;
    end
    

    switch plotType
        case 'ion'
            % ----- ion specific data -----
            Vi = mms.get_data('Vi_gse_fpi_brst_l2',tint,ic);
            
            
            try
                iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
                iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);
                % ignore flux where count is 1 (also makes function faster)
                iPDist.data(iPDist.data<1.1*iPDistErr.data) = 0;
            catch
                disp(['error reading PDist for ',tstr,', skipping...'])
                continue;
            end
            
        case 'electron'
            % ----- electron specific data -----
            Ve = mms.get_data('Ve_gse_fpi_brst_l2',tint,ic);
            Te = mms.get_data('Te_gse_fpi_brst_l2',tint,ic);
            
            if isempty(Te)
                disp(['no electron moments for ',tstr,', skipping...'])
                continue;
            end
            
            E = mms.get_data('E_gse_edp_fast_l2',tint,ic);
            Eb = mms.get_data('E_gse_edp_brst_l2',tint,ic);
            TeFac = mms.rotate_tensor(Te,'fac',B);
            TePerp = (TeFac.yy+TeFac.zz)/2;
            TePar = TeFac.xx;
            
            try
                ePDist = mms.get_data('PDe_fpi_brst_l2',tint,ic);
            catch
                disp(['error reading PDist for ',tstr,', skipping...'])
                continue;
            end
            
        case 'ioncomp'
            % ----- ion composition specific data -----
            % HPCA distributions PSD
            hpcaModeTemp = hpcaMode;
            dataRead = 0;
            if strcmp(hpcaMode,'auto')
                hpcaMode = 'brst';
            end
            
            
            while ~dataRead
                c_eval(['hpDist = mms.db_get_ts(''mms?_hpca_',hpcaMode,'_l2_ion'',''mms?_hpca_hplus_phase_space_density'',tint);'],ic)
                c_eval(['he2pDist = mms.db_get_ts(''mms?_hpca_',hpcaMode,'_l2_ion'',''mms?_hpca_heplusplus_phase_space_density'',tint);'],ic)
                c_eval(['opDist = mms.db_get_ts(''mms?_hpca_',hpcaMode,'_l2_ion'',''mms?_hpca_oplus_phase_space_density'',tint);'],ic)
                
                c_eval(['hpca_azimuth = mms.db_get_ts(''mms?_hpca_',hpcaMode,'_l2_ion'',''mms?_hpca_azimuth_angles_degrees'',tint);'],ic)
                
                if isempty(hpDist) && strcmp(hpcaMode,'brst') % avoid infinite loop
                    hpcaMode = 'srvy';
                else
                    dataRead = 1;
                end    
            end
            
            if isempty(hpDist)
                disp('No HPCA data available, skipping...')
            end
            
            % placeholder
            hpca_elevation = transpose([123.750000000000;101.250000000000;78.7500000000000;56.2500000000000;33.7500000000000;11.2500000000000;11.2500000000000;33.7500000000000;56.2500000000000;78.7500000000000;101.250000000000;123.750000000000;146.250000000000;168.750000000000;168.750000000000;146.250000000000]);
            hpca_energy = [1.3550000;1.5718000;1.8428000;2.2221999;2.6015999;3.0894001;3.6314001;4.2818003;5.0406003;5.9620004;6.9917998;8.2384005;9.7559996;11.490399;13.550000;15.989000;18.861601;22.276201;26.232801;30.948200;36.530800;43.089001;50.785400;59.945202;70.676804;83.413803;98.373001;116.04220;136.85500;161.46181;190.45880;224.65901;264.98380;312.57138;368.72260;434.95502;513.05725;605.19720;713.86823;842.05121;993.32343;1171.6956;1382.1000;1630.2819;1923.0702;2268.4326;2675.7998;3156.2830;3723.1064;4391.7178;5180.4360;6110.7246;7208.1123;8502.5713;10029.493;11830.613;13955.199;16461.354;19417.529;22904.596;27017.889;31869.818;37593.121];     
            
            try
                iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
                iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);
                % ignore flux where count is 1 (also makes function faster)
                iPDist.data(iPDist.data<1.1*iPDistErr.data) = 0;
            catch
                disp(['error reading PDist for ',tstr,', skipping...'])
                continue;
            end
            
    end
    
    % check if all data is there
    if isempty(B1) || isempty(B2) || isempty(B3) || isempty(B4) || isempty(B) || ...
            isempty(ni) || isempty(ne) 
        disp(['no data for ',tstr,', skipping...'])
        continue;
    end
    
    
    %% get model bs normal
    scd = [];
    scd.Bu = zeros(1,3); scd.Vu = zeros(1,3); scd.nu = 0;
    scd.Bd = zeros(1,3); scd.Vd = zeros(1,3); scd.nd = 0;
    scd.R = R;
    nst = irf_shock_normal(scd);
    nvec = nst.n.(shModel);
    
    
    %% try to read omni data
    
    ff = irf_get_data(tint,'bx,by,bz,Ma,v,n','omni_min');
    
    if ~isempty(ff)
        Bu = nanmean(ff(:,2:4),1);
        if isnan(Bu); Bu = nan(1,3); end
        Ma = nanmean(ff(:,5));
        Vu = nanmean(ff(:,6));
        Nu = nanmean(ff(:,7));
        thBn = acosd(dot(Bu,nvec)/(norm(Bu)));
        if thBn>90 
            thBn = 180-thBn;
        end
    else
        disp('failed to read omni data, moving on with plot...')
        thBn = nan;
        Ma = nan;
        Bu = nan(1,3);
        Nu = nan;
        Vu = nan;
    end
    
    %for good measures, remove old ff
    clear ff
    
    
    %% set title strings
    meanRx = mean(R.gseR1(:,1)+R.gseR2(:,1)+R.gseR3(:,1)+R.gseR4(:,1))/4/u.RE*1e3;
    meanRy = mean(R.gseR1(:,2)+R.gseR2(:,2)+R.gseR3(:,2)+R.gseR4(:,2))/4/u.RE*1e3;
    meanRz = mean(R.gseR1(:,3)+R.gseR2(:,3)+R.gseR3(:,3)+R.gseR4(:,3))/4/u.RE*1e3;
    Rstr = [num2str(round(meanRx,1)),',',num2str(round(meanRy,1)),',',num2str(round(meanRz,1))];
    nstr = [num2str(round(nvec(1),2)),',',num2str(round(nvec(2),2)),',',num2str(round(nvec(3),2))];
    
    Bustr = [num2str(round(Bu(1),2)),',',num2str(round(Bu(2),2)),',',num2str(round(Bu(3),2))];
    Vustr = num2str(round(Vu));
    Nustr = num2str(round(Nu));
    thstr = num2str(round(thBn));
    Mastr = num2str(round(Ma));
    
    if isnan(Bu(1)); Bustr = '-'; end
    if isnan(Vu); Vustr = '-'; end
    if isnan(Nu); Nustr = '-'; end
    if isnan(thBn); thstr = '-'; end
    if isnan(Ma); Mastr = '-'; end
    
    %% reduce fi
    % fpi
    if strcmp(plotType,'ion') || strcmp(plotType,'ioncomp')
        f1Dn = iPDist.reduce('1D',nvec,'vg',vg,'nMC',nMC);
    end
    
    % hpca
    if strcmp(plotType,'ioncomp')
        % made up coordinates
        t1vec = cross([0,0,1],nvec)/norm(cross([0,0,1],nvec));
        t2vec = cross(nvec,t1vec)/norm(cross(nvec,t1vec));
        
        [hpred,~] = hpca_reduce_dist(hpDist,hpca_energy,hpca_azimuth,hpca_elevation,'nMC',nMC,'vg',vg,'xyz',[nvec;t1vec;t2vec],'m',1);
        [he2pred,~] = hpca_reduce_dist(he2pDist,hpca_energy,hpca_azimuth,hpca_elevation,'nMC',nMC,'vg',vg,'xyz',[nvec;t1vec;t2vec],'m',2);
    end
    
    
    %% plot

    switch plotType
        case 'ion'
            ionShockOverview
        case 'electron'
            eleShockOverview
        case 'ioncomp'
            ioncompShockOverview
    end
    
    %% add title
    
    ll = -0.14; % legend left
    lr = 1.14; % right
    
    irf_legend(h(1),descStr(1:end-1),[ll,1.35],'Fontsize',15,'interpreter','latex','horizontalalignment','left')
    irf_legend(h(1),['FOM: ',fomStr],[lr,1.35],'Fontsize',15,'interpreter','latex','horizontalalignment','right')
    
    irf_legend(h(1),['$\mathbf{R} = [',Rstr,']\,R_E$'],[ll,1.2],'Fontsize',15,'interpreter','latex','horizontalalignment','left')
    irf_legend(h(1),['$\mathbf{\hat{n}} = [',nstr,']$'],[lr,1.2],'Fontsize',15,'interpreter','latex','horizontalalignment','right')
   
    irf_legend(h(1),['$\mathbf{B}_u = [',Bustr,']$\,nT,'],[ll,1.045],'Fontsize',14,'interpreter','latex','horizontalalignment','left')
    irf_legend(h(1),['$N_u = ',Nustr,'$\,cm$^{-3}$, $V_u = ',Vustr,'$\,km/s, $M_A = ',Mastr,'$, $\theta_{Bn} = ',thstr,'^{\circ}$'],[lr,1.045],'Fontsize',14,'interpreter','latex','horizontalalignment','right')
    %% save figure
    % give it the time stamp of the start time
    irf_print_fig(h(1),[filePath,fileName(1:end-4)],'png')
    
    %% close figure
    close(h(1).Parent)
    
end
