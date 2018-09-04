% get and plot shock parameters for all events in list

%% Set parameters

% If data is read, ask to save it
if ~doLoadData
    saveParameters = irf_ask('Save parameters? (0=no, 1=yes) [%]>','saveParameters',0);
    
    dataNameGuess = [listFileName(1:find(ismember(listFileName,'_'))-1),'Params'];
    if saveParameters
        fileName = irf_ask('File name: [%]>','fileName',dataNameGuess);
    end
else
    saveParameters = 0;
end

eisMode = irf_ask('EIS mode (auto/srvy/brst) [%]>','eisMode','srvy');

% normal bs model to use (cannot be slho)
shModel = 'farris';

% number of events (should be more than actual events)
N = 1000;

%% initiate arrays
% only read data if data is not loaded
if ~doLoadData
    disp('Getting shock parameters')
    
    %TV = irf.time_array('2000-01-01T00:00:00Z',zeros(1,N));
    TV = zeros(N,1);
    dTV = zeros(N,1);
    MaV = zeros(N,1);
    VuV = zeros(N,1);
    thBnV = zeros(N,1);
    thBrV = zeros(N,1);
    thVnV = zeros(N,1);
    accEffV = zeros(N,1);
    EmaxV = zeros(N,1);
    hasEISV = zeros(N,1);
    RV = zeros(N,3);
    sigV = zeros(N,1);
    dstV = zeros(N,1);
    kpV = zeros(N,1);
    ssnV = zeros(N,1);
    s107V = zeros(N,1);
    lineNumV = zeros(N,1);
    
    
    %% Read line  
    
    tline = 1;
    
    fid = fopen(listFileName);
    
    lineNum = 0;
    % loop through skipped lines
    for ii = 1:startLine-1
        lineNum = lineNum+1;
        tline = fgets(fid);
    end
    
    count = 1;
    % super trooper looper
    while tline ~= -1
        
        %% read line from file
        tline = fgets(fid);
        lineNum = lineNum+1;
        
        if tline(1) == -1 || ~strcmp(tline(1),'2')
            disp('done!')
            break;
        end
        tintStr = [tline(1:10),'T',tline(12:19),'/',tline(23:32),'T',tline(34:41)];
        
        %% get time interval
        tint = irf.tint(tintStr);
        
        tstr = tint(1).toUtc;
        tstr([5,8,14,17]) = '';
        tstr = tstr(1:15);
        
        
        %% read data
        
        % ---------- FPI-DIS ----------
        % read iPDist for MMSic
        try
            iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
            iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);
            % ignore flux where count is 1
            iPDist.data(iPDist.data<1.1*iPDistErr.data) = 0;
        catch
            disp(['error reading PDist for ',tstr,', skipping...'])
            continue;
        end
        % must have FPI data
        if isempty(iPDist)
            disp(['no FPI data for ',tstr,', skipping...'])
            continue;
        end
        
        % ---------- SC position ----------
        % try to read sc position data
        try
            R = mms.get_data('R_gse',tint);
        catch
            disp(['error reading sc position for ',tstr,', skipping...'])
            continue;
        end
        
        
        % ---------- EPD-EIS ----------
        % read EIS data, keep going even if there is no data
        switch eisMode
            case 'srvy'
                c_eval('EISdpf! = mms.db_get_ts(''mms?_epd-eis_srvy_l2_phxtof'',''mms?_epd_eis_phxtof_proton_P4_flux_t!'',tint);',1:4,0:5)
                % assume all energy tables are the same
                [filepath,filename] = mms.get_filepath('mms2_epd-eis_srvy_l2_phxtof',tint(1));
                do = dataobj([filepath,filename]);
                c_eval('Eeis = 1e3*do.data.mms?_epd_eis_phxtof_proton_t0_energy.data;',ic) % in eV
                c_eval('dEeisMinus = 1e3*do.data.mms?_epd_eis_phxtof_proton_t0_energy_dminus.data;',ic) % in eV
                c_eval('dEeisPlus = 1e3*do.data.mms?_epd_eis_phxtof_proton_t0_energy_dplus.data;',ic) % in eV
                
            case 'brst' % burst is broken?
                error('EIS burst not finished, energy delta plus/minus missing')
                c_eval('EISdpf! = mms.db_get_ts(''mms?_epd-eis_brst_l2_phxtof'',''mms?_epd_eis_brst_phxtof_proton_P4_flux_t!'',tint);',1:4,0:5)
                % assume all energy tables are the same
                [filepath,filename] = mms.get_filepath('mms2_epd-eis_srvy_l2_phxtof',tint(1));
                do = dataobj([filepath,filename]);
                Eeis = 1e3*do.data.mms2_epd_eis_brst_phxtof_proton_t0_energy.data;
        end
        % hack to fix bug if tint goes over two days
        t1str = tint(1).toUtc; t2str = tint(2).toUtc; 
        
        % check if EIS data was read, continue either way
        if isempty(EISdpf0) % skip
            hasEIS = 0;% do not include EIS in acceff
        elseif strcmp(t1str(1:10),t2str(1:10)) % day is not the same, fix should be better
            hasEIS = 0;
        else % convert EIS data to PSD
            hasEIS = 1;% do include EIS in acceff
            mm = 1; % ion mass
            % energy table for each time step and energy index
            ETabEis = repmat(Eeis,EISdpf0.length,1);
            % copy objects
            c_eval('EISpsd? = EISdpf?;',0:5)
            % sort of magic conversion :(
            c_eval('EISpsd?.data = EISdpf?.data/1e12*mm^2*0.53707./ETabEis;',0:5)
            % mysterious correction factor :(
            c_eval('EISpsd?.data = EISpsd?.data*1e-3;',0:5)
            % average EIS psd from all detectors (same time stamps)
            EISpsd = (EISpsd0+EISpsd1+EISpsd2+EISpsd3+EISpsd4+EISpsd5)/6;
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
        ff2 = irf_get_data(tint+[-1,1]*3.6e3,'dst,kp,ssn,f10.7','omni2');
        
        if ~isempty(ff)
            Bu = nanmean(ff(:,2:4),1);
            if isnan(Bu); Bu = nan(1,3); end
            Ma = nanmean(ff(:,5));
            Vu = nanmean(ff(:,6));
            Nu = nanmean(ff(:,7));
            thBn = acosd(dot(Bu,nvec)/(norm(Bu)));
            
            % other info
            dst = nanmean(ff2(:,2));
            kp = nanmean(ff2(:,3));
            ssn = nanmean(ff2(:,4));
            s107 = nanmean(ff2(:,5));
            
            
            if thBn>90
                thBn = 180-thBn;
            end
            % guess the solar wind is in x direction (~4 deg wrong)
            thVn = acosd(nvec(1));
            if thVn>90
                thVn = 180-thVn;
            end
            
            thBr = acosd(Bu(1)/norm(Bu));
            if thBr>90
                thBr = 180-thBr;
            end
            
        else
            disp('failed to read omni data, moving on with plot...')
            thBn = nan;
            Ma = nan;
            Bu = nan(1,3);
            Nu = nan;
            Vu = nan;
                        
            % other info
            dst = nan;
            kp = nan;
            ssn = nan;
            s107 = nan;
        end
        
        %for good measures, remove old ff
        clear ff
        
        
        %% Merge FPI and EIS into one data product 
        % This should really not be done for alternating energy tables
        
        % FPI energy matrix
        emat = double(iPDist.energy); % in eV
        
        if hasEIS
            % initiate matrices (unknown size so use cells)
            Fcomb = cell(1,EISpsd.length);
            Ecomb = cell(1,EISpsd.length);
            dEcombMinus = cell(1,EISpsd.length);
            dEcombPlus = cell(1,EISpsd.length);
            
            
            % time difference between EIS measurements, assume constant
            dt = median(diff(EISpsd.time.epochUnix));
            for it = 1:EISpsd.length
                disp([num2str(it),'/',num2str(EISpsd.length)])
                
                % time of this time step
                t1 = EISpsd.time(it).epochUnix;
                % eis psd array for one time step
                Feis = double(EISpsd.data(it,:));
                % time indicies of FPI that fall within EIS time
                idFpi = find(iPDist.time.epochUnix>=t1 & iPDist.time.epochUnix<t1+dt);
                % average over fpi times
                Ffpi = nanmean(double(iPDist.convertto('s^3/m^6').omni.data(idFpi,:)),1);
                % FPI energy in [eV] (not good if esteptable is used)
                Efpi = mean(emat(idFpi,:),1);
                % delta energy of FPI [eV] ()
                dEfpi = double(iPDist.ancillary.delta_energy_plus(1,:))*2; 
                
                
                % ---- combine data ----
                % now, use entire energy range of FPI, could be chaged
                % later
                % find index of first EIS data point (last lower edge under FPI max E)
                idEisFirst = find(Eeis-dEeisMinus<Efpi(end)+dEfpi(end)/2,1,'last');
                % get "first" dE minus of EIS part (can be negative)
                dEeisMinus_first = Eeis(idEisFirst)-(Efpi(end)+dEfpi(end)/2);
                % combined dE minus
                dEcombMinus{it} = [dEfpi/2,dEeisMinus_first,dEeisMinus(idEisFirst+1:end)];
                % combined dE plus
                dEcombPlus{it} = [dEfpi/2,dEeisPlus(idEisFirst:end)];
                % combined E
                Ecomb{it} = [Efpi,Eeis(idEisFirst:end)];                
                % combined psd
                Fcomb{it} = [Ffpi,Feis(idEisFirst:end)];
                
            end
        end
        
        %% get energy flux of ions with E>10Esw
        % particle mass
        M = u.mp;
        
        % solar wind energy in J
        Esw = .5*M*(Vu*1e3)^2;
        
        % number of energy bins
        if hasEIS; nE = length(Fcomb{1}); else; nE = length(iPDist.depend{1}); end
        
        % number of time steps
        if hasEIS; nT = length(Fcomb); else; nT = iPDist.length; end
        
        % total energy flux
        EF = zeros(iPDist.length,1);
        % energetic energy flux
        EFen = zeros(iPDist.length,1);
        
        % ------- time loop :D --------
        for it = 1:nT
            disp([num2str(it),'/',num2str(nT)])
            
            % 1d data matrix of PSD for time index it [s^3/m^6]
            if hasEIS
                Fpsd = Fcomb{it}; 
            else
                Fpsd = double(iPDist.convertto('s^3/m^6').omni.data(it,:)); 
            end
            
            % energy in [J]
            if hasEIS; E = Ecomb{it}*u.e; else; E = emat(it,:)*u.e; end
            
            % eelta energy [J]
            if hasEIS
                dE = (dEcombMinus{it}+dEcombPlus{it})*u.e;
            else
                % assumes energy table is the same for all time steps
                % also assumes energy_delta_plus/minus are the same
                dE = double(iPDist.ancillary.delta_energy_plus(1,:))*2*u.e;
            end
            
            % average energy flux per energy level [J/m^3]
            dEF = E.*Fpsd.*dE;
            
            % total energy flux
            EF(it) = 4*pi*sqrt(2/M^3)*sum(dE.*sqrt(E.^3).*Fpsd);
            
            % --- energetic energy flux ---
            % first find last lower bin edge which is less than 10Esw
            idll = find(E-dE/2<10*Esw,1,'last');
            % if idll is empty, then probably 10Esw is greater than
            % instrument limit, maybe deal with?
            % then get "first" dE of energetic part
            dEen_first = E(idll)+dE(idll)/2-10*Esw;
            % array of dEs of energetic part
            dEen = [dEen_first,dE(idll+1:end)];
            % array of Es of energetic part
            Een = E(idll:end);
            % array of psd of energetic part
            Fpsden = Fpsd(idll:end);
            
            % energetic energy flux
            EFen(it) = 4*pi*sqrt(2/M^3)*sum(dEen.*sqrt(Een.^3).*Fpsden);
            
            
        end
        % 2 alternatives, not sure which is best
        % accEff = mean(EFen)/mean(EF);
        accEff = mean(EFen)/(Nu*Esw*1e6);
        
        
        %% set values
        % fill arrays
        TV(count) = tint(1).epochUnix;
        dTV(count) = diff(tint.epochUnix);
        MaV(count) = Ma;
        VuV(count) = Vu;
        thBnV(count) = thBn;
        thVnV(count) = thVn;
        thBrV(count) = thBr;
        accEffV(count) = accEff;
        % Emax is not perfect for alternating energy table but who cares?
        EmaxV(count) = max(E);
        hasEISV(count) = hasEIS;
        % mean of all four
        RV(count,:) =mean((R.gseR1(:,1:3)+R.gseR2(:,1:3)+R.gseR3(:,1:3)+R.gseR4(:,1:3))/4)/u.RE*1e3;
        % compression factor of models
        sigV(count) = nst.info.sig.(shModel);
        % dst index
        sigV(count) = nst.info.sig.(shModel);
        
        dstV(count) = dst;
        kpV(count) = kp;
        ssnV(count) = ssn;
        s107V(count) = s107;
        
        lineNumV = lineNum;
        
        disp(['Actually completed one, count = ',num2str(count),' lineNumber = ',num2str(lineNum)])
        count = count+1;
        
        %% end loop if requested
        if lineNum >= stopLine; disp('Reached stop line, exiting...'); break; end
        
        
    end
    
    disp('Done calculating parameters!')
    
end

%% Clean arrays
dTV = dTV(TV~=0);
MaV = MaV(TV~=0);
VuV = VuV(TV~=0);
thBnV = thBnV(TV~=0);
thBrV = thBrV(TV~=0);
thVnV = thVnV(TV~=0);
accEffV = accEffV(TV~=0,:);
EmaxV = EmaxV(TV~=0,:);
hasEISV = hasEISV(TV~=0,:);
RV = RV(TV~=0,:);
sigV = sigV(TV~=0);
dstV = dstV(TV~=0);
kpV = kpV(TV~=0);
ssnV = ssnV(TV~=0);
s107V = s107V(TV~=0);
lineNumV = lineNumV(TV~=0);

N = numel(TV(TV~=0));

TV = TV(TV~=0);


%% Save parameters if requested

if saveParameters
    disp('Saving parameters...')
    save(fileName,'dTV','MaV','VuV','thBnV','thVnV','accEffV','EmaxV','hasEISV','RV','sigV','TV','N','dstV','kpV','ssnV','s107V','thBrV','lineNumV')
    disp('saved!')
end
