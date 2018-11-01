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

% normal bs models to use (see irf_shock_normal for explanation)
shModel = {'farris','slho','per','fa4o','fan4o','foun'};

% number of events (should be more than actual events)
N = 1000;

%% initiate arrays
% only read data if data is not loaded
if ~doLoadData
    disp('Getting shock parameters')
    
    % single values
    TV = zeros(N,1); % start time
    dTV = zeros(N,1); % interval duration
    VuV = zeros(N,1); % upstream speed
    thBrV = zeros(N,1); % magnetic field/radial angle
    betaiV = zeros(N,1); % ion plasma beta
    accEffV = zeros(N,1); % acceleraton efficiency
    accEffFpiV = zeros(N,1); % acceleraton efficiency (only FPI)
    accEffAltV = zeros(N,1); % acceleraton efficiency (alternative)
    accEffAltFpiV = zeros(N,1); % acceleraton efficiency (alternative, only FPI)
    EmaxV = zeros(N,1); % maximum energy used [eV]
    EfpiMaxV = zeros(N,1); % maximum energy of FPI-DIS [eV]
    hasEISV = zeros(N,1); % boolean if EIS exists or not
    RV = zeros(N,3); % sc position
    dstV = zeros(N,1); % DST index
    kpV = zeros(N,1); % Kp index
    ssnV = zeros(N,1); % sunspot number
    s107V = zeros(N,1); % F10.7 index
    aeV = zeros(N,1); % AE index
    lineNumV = zeros(N,1); % line number in text file
    
    % structures (depend on shock model)
    nvecV = [];
    MaV = []; % Alfven Mach number
    MfV = []; % magnetosonic/fast Mach number
    thBnV = []; % shock angle
    thVnV = []; % flow incidence angle
    % loop through shModel array
    for jj = 1:length(shModel)
        nvecV.(shModel{jj}) = zeros(N,3); % shock normal vector
        MaV.(shModel{jj}) = zeros(N,1); % Alfven Mach number
        MfV.(shModel{jj}) = zeros(N,1); % magnetosonic/fast Mach number
        thBnV.(shModel{jj}) = zeros(N,1); % shock angle
        thVnV.(shModel{jj}) = zeros(N,1); % flow incidence angle
        sigV.(shModel{jj}) = zeros(N,1); % bow shock compression factor
    end
    
    
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
                try
                    [filepath,filename] = mms.get_filepath('mms2_epd-eis_srvy_l2_phxtof',tint(1));
                    do = dataobj([filepath,filename]);
                    c_eval('Eeis = 1e3*do.data.mms?_epd_eis_phxtof_proton_t0_energy.data;',ic) % in eV
                    c_eval('dEeisMinus = 1e3*do.data.mms?_epd_eis_phxtof_proton_t0_energy_dminus.data;',ic) % in eV
                    c_eval('dEeisPlus = 1e3*do.data.mms?_epd_eis_phxtof_proton_t0_energy_dplus.data;',ic) % in eV
                catch
                    disp('error reading EIS data, continuing... ')
                end
                
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
        elseif ~strcmp(t1str(1:10),t2str(1:10)) % day is not the same, fix should be better
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
         
        
        %% try to read omni data
                
        % initialize empty data structure for irf_shock_normal
        scd = [];
        scd.Bu = zeros(1,3); scd.Vu = zeros(1,3); scd.nu = 0;
        scd.Bd = zeros(1,3); scd.Vd = zeros(1,3); scd.nd = 0;
        scd.R = R;
        
        ff = irf_get_data(tint,'bx,by,bz,Ma,vx,vy,vz,n,T','omni_min');
        ff2 = irf_get_data(tint+[-1,1]*3.6e3,'dst,kp,ssn,f10.7,ae','omni2');
        
        if ~isempty(ff)
            Bu = nanmean(ff(:,2:4),1);
            if isnan(Bu); Bu = nan(1,3); end
            Ma = nanmean(ff(:,5)); % Alfven Mach number 
            Vu = nanmean(ff(:,6:8),1); % velocity [km/s]
            if isnan(Vu); Vu = nan(1,3); end
            Nu = nanmean(ff(:,9)); % density [/cc]
            Tu = nanmean(ff(:,10))/u.e*u.kB; % temperature [eV]
            
            % put real values in structure
            scd.Bu = Bu; scd.Vu = Vu; scd.nu = Nu;
            % get model bs normals
            nst = irf_shock_normal(scd);
            % nvec is a structure with fields defined by shModel
            nvec = [];
            for jj = 1:length(shModel)
                nvec.(shModel{jj}) = nst.n.(shModel{jj});
            end
            
            % calculate shock angle
            % thBn is a structure with fields defined by shModel
            thBn = [];
            for jj = 1:length(shModel)
                thBn.(shModel{jj}) = acosd(dot(Bu,nvec.(shModel{jj}))/(norm(Bu)));
                if thBn.(shModel{jj})>90
                    thBn.(shModel{jj}) = 180-thBn.(shModel{jj});
                end
            end
            
            % calculate sw incidence angle
            % thVn is a structure with fields defined by shModel
            thVn = [];
            for jj = 1:length(shModel)
                thVn.(shModel{jj}) = acosd(dot(Vu,nvec.(shModel{jj}))/(norm(Vu)));
                if thVn.(shModel{jj})>90
                    thVn.(shModel{jj}) = 180-thVn.(shModel{jj});
                end
            end
            
            % calculate B to radial angle
            thBr = acosd(Bu(1)/norm(Bu));
            if thBr>90
                thBr = 180-thBr;
            end
            
            % ion beta
            betai = Nu*1e6*Tu*u.e/(norm(Bu*1e-9)^2/(2*u.mu0));
            
            % calculate fast Mach number
            vA = norm(Vu)/Ma*1e3; % Alfven speed [m/s]
            % assume ion and electron temperatures are equal
            cs = sqrt(4*Tu*u.e/u.mp); % sound speed [m/s]
            % magnetosonic speed perp to B
            cms = sqrt(cs^2+vA^2);
            
            % loop through shModel array
            vms = [];
            Mf = [];
            for jj = 1:length(shModel)
                % magnetosonic group speed along normal
                vms.(shModel{jj}) = sqrt(cms^2/2+sqrt(cms^4/4-vA^2*cs^2*cosd(thBn.(shModel{jj}))^2));
                % fast Mach number of stationary shock (not solar wind)
                Mf.(shModel{jj}) = dot(-Vu,nvec.(shModel{jj}))/vms.(shModel{jj})*1e3;
            end
            
            % other info
            dst = nanmean(ff2(:,2));
            kp = nanmean(ff2(:,3));
            ssn = nanmean(ff2(:,4));
            s107 = nanmean(ff2(:,5));
            ae = nanmean(ff2(:,6));
            
        else
            disp('failed to read omni data, moving on with plot...')
            thBn = nan;
            Ma = nan;
            Bu = nan(1,3);
            Nu = nan;
            Vu = nan(1,3);
            Tu = nan;
                        
            % other info
            dst = nan;
            kp = nan;
            ssn = nan;
            s107 = nan;
            ae = nan;
        end
        
        %for good measures, remove old ff
        %clear ff
        
        
        %% Merge FPI and EIS into one data product 
        % This should really not be done for alternating energy tables
        
        % FPI energy matrix
        emat = double(iPDist.energy); % in eV
        
        FfpiMat = double(iPDist.convertto('s^3/m^6').omni.data);
        
        if hasEIS
            % initiate matrices (unknown size so use cells)
            Fcomb = cell(1,EISpsd.length);
            Ecomb = cell(1,EISpsd.length);
            dEcombMinus = cell(1,EISpsd.length);
            dEcombPlus = cell(1,EISpsd.length);
            
            % time difference between EIS measurements, assume constant
            dt = median(diff(EISpsd.time.epochUnix));
            for it = 1:EISpsd.length
                
                % time of this time step
                t1 = EISpsd.time(it).epochUnix;
                % eis psd array for one time step
                Feis = double(EISpsd.data(it,:));
                % time indicies of FPI that fall within EIS time
                idFpi = find(iPDist.time.epochUnix>=t1 & iPDist.time.epochUnix<t1+dt);
                % average over fpi times
                Ffpi = nanmean(FfpiMat(idFpi,:),1);
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
        
        % record FPI max energy [eV]
        EfpiMaxV(count) = max(Efpi);
        
        %% get energy flux of ions with E>10Esw
        % particle mass
        M = u.mp;
        
        % solar wind energy in [J]
        Esw = .5*M*(norm(Vu)*1e3)^2;
        
        % number of energy bins
        if hasEIS; nEEis = length(Fcomb{1}); else; nEEis = 0; end
        nEFpi = length(iPDist.depend{1});
        
        % number of time steps
        if hasEIS; nTEis = length(Fcomb); else; nTEis = 0; end
        nTFpi = iPDist.length;
        
        % total energy flux
        EF = zeros(nTEis,1);
        % energetic energy density
        EFen = zeros(nTEis,1);
        % energetic energy density using only FPI
        EFenFpi = zeros(nTFpi,1);
        % total ion energy density using only FPI
        EFFpi = zeros(nTFpi,1);
        
        % ------- time loop :D --------
        % For FPI+EIS
        if hasEIS
            for it = 1:nTEis
              
                % 1d data matrix of PSD for time index it [s^3/m^6]
                Fpsd = Fcomb{it};
                
                % energy in [J]
                E = Ecomb{it}*u.e;
                
                % delta energy [J]
                dE = (dEcombMinus{it}+dEcombPlus{it})*u.e;
                
                % average energy flux per energy level [J/m^3]
                %dEF = E.*Fpsd.*dE;
                
                % total energy density [J/m^3]
                EF(it) = 4*pi*sqrt(2/M^3)*sum(dE.*sqrt(E.^3).*Fpsd);
                
                % --- energetic energy density ---
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
                % array of Es of energetic part using only FPI (assume nE=32)
                EenFpi = E(idll:32);
                % array of Es using only FPI (assume nE=32)
                EFpi = E(1:32);
                % array of psd of energetic part
                Fpsden = Fpsd(idll:end);
                
                % energetic energy density
                EFen(it) = 4*pi*sqrt(2/M^3)*sum(dEen.*sqrt(Een.^3).*Fpsden);

            end
        end
        
        % For FPI only
        for it = 1:nTFpi
            
            % 1d data matrix of PSD for time index it [s^3/m^6]
            FpsdFpi = FfpiMat(it,:);
            
            % energy in [J]
            EFpi = emat(it,:)*u.e;
            
            % delta energy [J]
            % assumes energy table is the same for all time steps
            % also assumes energy_delta_plus/minus are the same
            dEFpi = double(iPDist.ancillary.delta_energy_plus(1,:))*2*u.e;
            
            % average energy density per energy level [J/m^3]
            %dEF = E.*Fpsd.*dE;
            
            % total energy density [J/m^3]
            EF(it) = 4*pi*sqrt(2/M^3)*sum(dEFpi.*sqrt(EFpi.^3).*FpsdFpi);
            

            
            % --- energetic energy density ---
            % first find last lower bin edge which is less than 10Esw
            idll = find(EFpi-dEFpi/2<10*Esw,1,'last');
            % if idll is empty, then probably 10Esw is greater than
            % instrument limit, maybe deal with?
            % then get "first" dE of energetic part
            dEen_first = EFpi(idll)+dE(idll)/2-10*Esw;
            
            
            % array of dEs of energetic part using only FPI (assume nE=32)
            dEenFpi = [dEen_first,dEFpi(idll+1:end)];
            
            % array of psd of energetic part using only FPI (assume nE=32)
            FpsdenFpi = FpsdFpi(idll:end);
            
            % energetic energy density using only FPI
            EFenFpi(it) = 4*pi*sqrt(2/M^3)*sum(dEenFpi.*sqrt(EenFpi.^3).*FpsdenFpi);
            % total ion energy density using only FPI
            EFFpi(it) = 4*pi*sqrt(2/M^3)*sum(dEFpi.*sqrt(EFpi.^3).*FpsdFpi);
            
            
        end
        
        % 2 ways, the alternative way is not really correct
        accEff = mean(EFen)/mean(EF);
        accEffAlt = mean(EFen)/(Nu*Esw*1e6);
        accEffFpi = mean(EFenFpi)/mean(EFFpi);
        accEffAltFpi = mean(EFenFpi)/(Nu*Esw*1e6);
        
        %% set values (fill arrays)
        % single values
        TV(count) = tint(1).epochUnix;
        dTV(count) = diff(tint.epochUnix);
        VuV(count) = norm(Vu);
        betaiV(count) = betai;
        thBrV(count) = thBr;
        accEffV(count) = accEff;
        accEffAltV(count) = accEffAlt;
        accEffFpiV(count) = accEffFpi;
        accEffAltFpiV(count) = accEffAltFpi;
        % Emax is not perfect for alternating energy table but who cares?
        % E max including EIS [eV]
        EmaxV(count) = max(E)/u.e;
        hasEISV(count) = hasEIS;
        % mean of all four
        RV(count,:) =mean((R.gseR1(:,1:3)+R.gseR2(:,1:3)+R.gseR3(:,1:3)+R.gseR4(:,1:3))/4)/u.RE*1e3;
        % dst index
        dstV(count) = dst;
        kpV(count) = kp;
        ssnV(count) = ssn;
        s107V(count) = s107;
        aeV(count) = ae;
        
        % structures
        for jj = 1:length(shModel)
            nvecV.(shModel{jj})(count,:) = nvec.(shModel{jj});
            MaV.(shModel{jj})(count) = Ma*cosd(thVn.(shModel{jj})); % of shock (not solar wind)
            MfV.(shModel{jj})(count) = Mf.(shModel{jj}); % of shock (not solar wind)
            thBnV.(shModel{jj})(count) = thBn.(shModel{jj});
            thVnV.(shModel{jj})(count) = thVn.(shModel{jj});
            % compression factor of models
            sigV.(shModel{jj})(count) = nst.info.sig.(shModel{jj});
        end
        
        lineNumV(count) = lineNum;
        
        disp(['Actually completed one, count = ',num2str(count),' lineNumber = ',num2str(lineNum)])
        count = count+1;
        
        %% end loop if requested
        if lineNum >= stopLine; disp('Reached stop line, exiting...'); break; end
        
        
    end
    
    disp('Done calculating parameters!')
    
end

%% Clean arrays
% single values
dTV = dTV(TV~=0);

VuV = VuV(TV~=0);

thBrV = thBrV(TV~=0);
betaiV = betaiV(TV~=0);
accEffV = accEffV(TV~=0,:);
accEffAltV = accEffAltV(TV~=0,:);
accEffFpiV = accEffFpiV(TV~=0,:);
accEffAltFpiV = accEffAltFpiV(TV~=0,:);
EmaxV = EmaxV(TV~=0,:);
EfpiMaxV = EfpiMaxV(TV~=0,:);
hasEISV = hasEISV(TV~=0,:);
RV = RV(TV~=0,:);
dstV = dstV(TV~=0);
kpV = kpV(TV~=0);
ssnV = ssnV(TV~=0);
s107V = s107V(TV~=0);
aeV = aeV(TV~=0);
lineNumV = lineNumV(TV~=0);

% structures
for jj = 1:length(shModel)
    nvecV.(shModel{jj}) = nvecV.(shModel{jj})(TV~=0,:);
    MaV.(shModel{jj}) = MaV.(shModel{jj})(TV~=0);
    MfV.(shModel{jj}) = MfV.(shModel{jj})(TV~=0);
    thVnV.(shModel{jj}) = thVnV.(shModel{jj})(TV~=0);
    thBnV.(shModel{jj}) = thBnV.(shModel{jj})(TV~=0);
    sigV.(shModel{jj}) = sigV.(shModel{jj})(TV~=0);
end

% now N is number of events
N = numel(TV(TV~=0));

% finally clean time array
TV = TV(TV~=0);


%% Save parameters if requested

if saveParameters
    disp('Saving parameters...')
    save(fileName,'dTV','MaV','MfV','VuV','thBnV','thBrV','thVnV','betaiV','accEffV','accEffAltV','accEffFpiV','accEffAltFpiV','EmaxV','EfpiMaxV','hasEISV','RV','sigV','TV','N','dstV','kpV','ssnV','s107V','aeV','lineNumV','nvecV','shModel')
    disp('saved!')
end
