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

eisMode = irf_ask('EIS mode (auto/srvy/brst/none) [%]>','eisMode','none');

% normal bs models to use (see irf_shock_normal for explanation)
shModel = {'farris','slho','per','fa4o','fan4o','foun'};

% number of events (should be more than actual events)
N = 1000;


%% read and get up and downstream times table

% Select data file to read
tintFileArr = dir('*.csv');

if isempty(tintFileArr)
    irf.log('c','No data files')
    
else
    preSelectedFile = 1;
    for ii = 1:length(tintFileArr)
        fprintf([num2str(ii),':     ',tintFileArr(ii).name,'\n'])
    end
    fprintf('\n')
    tintFileInd = irf_ask('Select tint file: [%]>','dataFileInd',preSelectedFile);
    tintFileName = tintFileArr(tintFileInd).name;
end

tintTab = readtable([pwd,'/',tintFileName]);

%% initiate arrays
% only read data if data is not loaded
if ~doLoadData
    disp('Getting shock parameters')
    
    % single values
    TV = zeros(N,1); % start time
    dTV = zeros(N,1); % interval duration
    TuV = zeros(N,1); % start time of upstream interval
    dTuV = zeros(N,1); % upstream interval duration
    TdV = zeros(N,1); % start time of downstream interval
    dTdV = zeros(N,1); % downstream interval duration
    VuV = zeros(N,3); % upstream velocity
    VuLV = zeros(N,3); % local upstream velocity
    VdV = zeros(N,3); % downstream velocity
    BuV = zeros(N,3); % upstream magnetic field vector
    BuLV = zeros(N,3); % local upstream magnetic field vector
    BdV = zeros(N,3); % downstream magnetic field vector
    NuV = zeros(N,1); % upstream number density
    NuLV = zeros(N,1); % local upstream number density
    NdV = zeros(N,1); % downstream number density
    thBrV = zeros(N,1); % magnetic field/radial angle
    betaiV = zeros(N,1); % ion plasma beta
    TiV = zeros(N,1); % upstream ion temperature
    fdV = zeros(N,32); % spherical mean of ion psd (only FPI)
    fuV = zeros(N,32); % spherical mean of ion psd (only FPI)
    EV = zeros(N,32);
    dEV = zeros(N,32);
    %fOmniComb = cell(N,1); % spherical mean of ion sph (combined FPI & EIS)
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
    imfIdV = zeros(N,1); % spacecraft index for B data omni
    swIdV = zeros(N,1); % spacecraft index for ion data omni
    
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
        
        %% get time intervals
        tint = irf.tint(tintStr);
        
        tstr = tint(1).toUtc;
        tstr([5,8,14,17]) = '';
        tstr = tstr(1:15);
        
        tintu = irf_time(tintTab.Var1(lineNum)+[tintTab.Var2(lineNum);tintTab.Var3(lineNum)],'epoch>epochTT');
        tintd = irf_time(tintTab.Var1(lineNum)+[tintTab.Var4(lineNum);tintTab.Var5(lineNum)],'epoch>epochTT');
        
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
            
            case 'none'
                % do not read any data
                EISdpf0 = [];
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
        
        % other data for local parameters (brst for good measure)
        Vi = mms.get_data('Vi_gse_fpi_brst_l2',tint,ic);
        Ni = mms.get_data('Ni_fpi_brst_l2',tint,ic);
        B = mms.get_data('B_gse_fgm_brst_l2',tint,1);
         
        
        %% try to read omni data
                
        % initialize empty data structure for irf_shock_normal
        scd = [];
        scd.Bu = zeros(1,3); scd.Vu = zeros(1,3); scd.nu = 0;
        scd.Bd = zeros(1,3); scd.Vd = zeros(1,3); scd.nd = 0;
        scd.R = R;
        
        ff = irf_get_data(tint,'bx,by,bz,Ma,vx,vy,vz,n,T,imfid,swid','omni_min');
        ff2 = irf_get_data(tint+[-1,1]*3.6e3,'dst,kp,ssn,f10.7,ae','omni2');
        
        if ~isempty(ff)
            Bu = nanmean(ff(:,2:4),1);
            if isnan(Bu); Bu = nan(1,3); end
            Ma = nanmean(ff(:,5)); % Alfven Mach number 
            Vu = nanmean(ff(:,6:8),1); % velocity [km/s]
            if isnan(Vu); Vu = nan(1,3); end
            Nu = nanmean(ff(:,9)); % density [/cc]
            Tu = nanmean(ff(:,10))/u.e*u.kB; % temperature [eV]
            imfId = nanmean(ff(:,11)); % not integers mean that sc switches
            swId = nanmean(ff(:,12)); 
            
            
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
            imfId = nan;
            swId = nan;
                        
            % other info
            dst = nan;
            kp = nan;
            ssn = nan;
            s107 = nan;
            ae = nan;
        end
        
        
        %% Get local up- and downstream parameters
        
        BuL = nanmean(B.tlim(tintu).data);
        BdL = nanmean(B.tlim(tintd).data);
        
        NuL = nanmean(Ni.tlim(tintu).data);
        NdL = nanmean(Ni.tlim(tintd).data);
        
        VuL = nanmean(Vi.tlim(tintu).data);
        VdL = nanmean(Vi.tlim(tintd).data);
        
        
        %% get psd of up- and downstream ions
        
        % downstream distribution
        iPDistDown = iPDist.tlim(tintd);
        iPDistUp = iPDist.tlim(tintu);
        
        % FPI energy matrix
        emat = double(iPDistDown.energy); % in eV
        % energy in [eV] (assumes all energy tables are the same as the first one)
        EFpi = emat(1,:);
        % delta energy of FPI [eV]
        dEFpi = double(iPDist.ancillary.delta_energy_plus(1,:))*2;
        
        % record FPI max energy [eV]
        EfpiMaxV(count) = EFpi(end)+dEFpi(end)/2;
        
        % get the psd from fpi in SI units
        FfpiMatDown = double(iPDistDown.convertto('s^3/m^6').omni.data);
        FfpiMatUp = double(iPDistUp.convertto('s^3/m^6').omni.data);
        % just save the average up- and downstream distruibution
        FdPsdFpi = mean(FfpiMatDown);
        FuPsdFpi = mean(FfpiMatUp);

        %% set values (fill arrays)
        % single values
        TV(count) = tint(1).epochUnix;
        dTV(count) = diff(tint.epochUnix);
        TuV(count) = tintu(1).epochUnix;
        dTuV(count) = diff(tintu.epochUnix);
        TdV(count) = tintd(1).epochUnix;
        dTdV(count) = diff(tintd.epochUnix);
        VuV(count,:) = Vu;
        VuLV(count,:) = VuL;
        VdV(count,:) = VdL;
        BuV(count,:) = Bu;
        BuLV(count,:) = BuL;
        BdV(count,:) = BdL;
        NuV(count) = Nu;
        NuLV(count) = NuL;
        NdV(count) = NdL;
        betaiV(count) = betai;
        thBrV(count) = thBr;
        fdV(count,:) = FdPsdFpi;
        fuV(count,:) = FuPsdFpi;
        EV(count,:) = EFpi;
        dEV(count,:) = dEFpi;
        % hasEISV(count) = hasEIS;
        % mean of all four
        RV(count,:) = mean((R.gseR1(:,1:3)+R.gseR2(:,1:3)+R.gseR3(:,1:3)+R.gseR4(:,1:3))/4)/u.RE*1e3;
        % dst index
        dstV(count) = dst;
        kpV(count) = kp;
        ssnV(count) = ssn;
        s107V(count) = s107;
        aeV(count) = ae;
        imfIdV(count) = imfId;
        swIdV(count) = swId;
        
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
TuV = TuV(TV~=0);
dTuV = dTuV(TV~=0);
TdV = TdV(TV~=0);
dTdV = dTdV(TV~=0);

VuV = VuV(TV~=0,:);
VdV = VdV(TV~=0,:);
VuLV = VuLV(TV~=0,:);
BuV = BuV(TV~=0,:);
BuLV = BuLV(TV~=0,:);
BdV = BdV(TV~=0,:);
NuV = NuV(TV~=0,:);
NuLV = NuLV(TV~=0,:);
NdV = NdV(TV~=0,:);
thBrV = thBrV(TV~=0);
betaiV = betaiV(TV~=0);
fdV = fdV(TV~=0,:);
fuV = fuV(TV~=0,:);
EV = EV(TV~=0,:);
dEV = dEV(TV~=0,:);
EfpiMaxV = EfpiMaxV(TV~=0,:);
hasEISV = hasEISV(TV~=0,:);
RV = RV(TV~=0,:);
dstV = dstV(TV~=0);
kpV = kpV(TV~=0);
ssnV = ssnV(TV~=0);
s107V = s107V(TV~=0);
aeV = aeV(TV~=0);
lineNumV = lineNumV(TV~=0);
swIdV = swIdV(TV~=0);
imfIdV = imfIdV(TV~=0);

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
    save(fileName,'dTV','TuV','dTuV','TdV','dTdV','MaV','MfV','VuV','VuLV','VdV','BuV','BuLV','BdV','NuV','NuLV','NdV','thBnV','thBrV','thVnV','betaiV','EV','dEV','fdV','fuV','EfpiMaxV','hasEISV','RV','sigV','TV','N','dstV','kpV','ssnV','s107V','aeV','lineNumV','nvecV','shModel','imfIdV','swIdV')
    disp('saved!')
end
