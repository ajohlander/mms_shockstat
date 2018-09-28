% get and plot shock parameters for all events in list

%% Set parameters

u = irf_units;

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
    RV = zeros(N,3);
    sigV = zeros(N,1);
    dstV = zeros(N,1);
    kpV = zeros(N,1);
    ssnV = zeros(N,1);
    s107V = zeros(N,1);
    lineNumberV = zeros(N,1);
    
    
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
        
        
        %% read some data
        
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
        
        if isempty(iPDist)
            disp(['no FPI data for ',tstr,', skipping...'])
            continue;
        end
        
        try
            R = mms.get_data('R_gse',tint);
        catch
            disp(['error reading sc position for ',tstr,', skipping...'])
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
        
        
        %% get energy flux of ions with E>10Esw
        
        nV = 32;
        nAz = 32;
        nEle = 16;
        
        M = u.mp;
        
        % solar wind energy in J
        Esw = .5*M*(Vu*1e3)^2;
        
        emat = double(iPDist.energy); % in eV
        
        dV = zeros(1,nV);
        % total energy flux
        EF = zeros(iPDist.length,1);
        % energetic energy flux
        EFen = zeros(iPDist.length,1);
        for it = 1:iPDist.length
            
            % 3d data matrix for time index it
            F3d = double(squeeze(double(iPDist.data(it,:,:,:)))*1e12); % s^3/m^6
            
            energy = emat(it,:);
            v = sqrt(2*energy*u.e/M); % m/s
            
            % azimuthal angle
            phi = double(iPDist.depend{2}(it,:)); % in degrees
            phi = phi*pi/180; % in radians
            
            % elevation angle
            th = double(iPDist.depend{3}); % polar angle in degrees
            th = th-90; % elevation angle in degrees
            th = th*pi/180; % in radians
            
            % diffs
            dV(2:end) = diff(v); dV(1) = dV(2); % quick and dirty
            dPhi = abs(median(diff(phi))); % constant
            dTh = abs(median(diff(th))); % constant
            
            % 3D matrices for instrumental bin centers
            TH = repmat(th,nV,1,nAz);       % [phi,th,v]
            TH = permute(TH,[1,3,2]);       % [v,phi,th]
            VEL = repmat(v,nAz,1,nEle);     % [phi,v,th]
            VEL = permute(VEL,[2,1,3]);     % [v,phi,th]
            DV = repmat(dV,nAz,1,nEle);     % [phi,v,th]
            DV = permute(DV,[2,1,3]);       % [v,phi,th]
            E = .5*M*VEL.^2; % energy [J]
            
            dE = E.*F3d.*VEL.^2.*cos(TH)*dPhi*dTh.*DV; % [J/m^3]
            
            % total energy flux
            EF(it) = sum(sum(sum(dE)));
            % energetic energy flux
            EFen(it) = sum(sum(sum(dE(E>10*Esw))));
            
            
            
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
        EmaxV(count) = max(energy);
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
RV = RV(TV~=0,:);
sigV = sigV(TV~=0);

dstV = dstV(TV~=0);
kpV = kpV(TV~=0);
ssnV = ssnV(TV~=0);
s107V = s107V(TV~=0);



TV = TV(TV~=0);


%% Save parameters if requested

if saveParameters
    disp('Saving parameters...')
    save(fileName,'dTV','MaV','VuV','thBnV','thVnV','accEffV','EmaxV','RV','sigV','TV','N','dstV','kpV','ssnV','s107V','thBrV')
    disp('saved!')
end

%% colors for plots
cmap = 'strangeways';
axcol = [1,1,1]*.4;
figcol = [1,1,1]*.2;
textcol = [1,1,1]*.95;
col1 = [253,232,159]/255;
col2 = [211,64,82]/255;


%% some bin edges
dthBin = 10;
thBinEdges = 0:dthBin:90;

%% approved data indices
vlim = 500;
dind = (~isnan(MaV) & VuV<vlim);

dTV = dTV(dind);
MaV = MaV(dind);
VuV = VuV(dind);
thBnV = thBnV(dind);
thBrV = thBrV(dind);
thVnV = thVnV(dind);
accEffV = accEffV(dind);
EmaxV = EmaxV(dind);
RV = RV(dind,:);
sigV = sigV(dind);

dstV = dstV(dind);
kpV = kpV(dind);
ssnV = ssnV(dind);
s107V = s107V(dind);

TV = TV(dind);

Nevents = numel(dind);

% angle between earth-sun line and sc position in xy plane
[alphaV,~] = cart2pol(RV(:,1),RV(:,2),RV(:,3));
alphaV = alphaV*180/pi; % degrees

% angle between earth-sun line and sc position
[phiV,~] = cart2pol(RV(:,1),sqrt(RV(:,2).^2+RV(:,3).^2));
phiV = phiV*180/pi; % degrees


%% Plot simple position
plotShockPos

%% plot parameter space
fig = figure;
hca = axes(fig);

scatter(hca,thBnV,MaV.*cosd(thVnV),400,'.')

hca.XLim = [0,90];
hca.YLim(1) = 0;

hca.Box = 'on';

ylabel(hca,'$M_A$','fontsize',15,'interpreter','latex')
xlabel(hca,'$\theta_{Bn}$ [$^{\circ}$]','fontsize',15,'interpreter','latex')

hca.LineWidth = 1.2;
hca.FontSize = 14;


%% Plot acceleration efficiency as a function of shock angle
fig = figure;
hca = axes(fig);

% plot events with Ma as color
scatter(hca,thBnV,accEffV*100,400,MaV.*cosd(thVnV),'.')
%scatter(hca,thBnV,accEffV*100,200,[1,1,1]*.2,'.')
%scatter(hca,thBnV(VuV<vlim),accEffV(VuV<vlim)*100,400,VuV(VuV<vlim),'.')
hold(hca,'on')
hca.XLim = [0,90];
hca.YLim(1) = 0;
plot(hca,45*[1,1],hca.YLim,'--','color',textcol,'linewidth',1.2)

sh_cmap(hca,cmap)
hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
%grid(hca,'on')

hcb = colorbar(hca);
hcb.Color = textcol;
hca.CLim = [3,20];

hca.Box = 'on';

ylabel(hca,'Acceleration efficiency [$\%$]','Fontsize',15,'interpreter','latex')
xlabel(hca,'$\theta_{Bn}$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')
ylabel(hcb,'$M_A$','Fontsize',15,'interpreter','latex')

title(hca,'Energy flux of ions with $E>10E_{sw}$ measured by MMS-FPI','Fontsize',15,'interpreter','latex','color',textcol)
hleg = irf_legend(hca,['$N = ',num2str(Nevents),'$'],[0.98,0.98],'Fontsize',15,'interpreter','latex','color',textcol);
hleg.BackgroundColor = hca.Color;

hca.LineWidth = 1.2;
hca.FontSize = 14;
hcb.LineWidth = 1.2;


%% same as above but with avg with proper Baysian limits
fig = figure;
hca = axes(fig);
hold(hca,'on')
hsc = scatter(hca,thBnV,accEffV*100,200,col2,'.');

% set significance
beta = .95;

idTh = discretize(thBnV,thBinEdges);
accEffAvg = zeros(1,length(thBinEdges)-1);
accEffStd = zeros(1,length(thBinEdges)-1);

accEffErrUp = zeros(1,length(thBinEdges)-1);
accEffErrDown = zeros(1,length(thBinEdges)-1);

for ii = 1:length(thBinEdges)-1
    % mean and std
    mu =  nanmean(accEffV(idTh==ii));
    sig = nanstd(accEffV(idTh==ii),1);
    
    % after a lot(!) of math
    accEffErrUp(ii) = norminv(beta*normcdf(mu/sig)-normcdf(mu/sig)+1)*sig;
    accEffErrDown(ii) = norminv(1-beta*normcdf(mu/sig))*sig;
    
    accEffAvg(ii) = mu;
    accEffStd(ii) = sig;
    
    % clean arrays
    accEffErrUp(isnan(accEffErrUp)) = 0;
    accEffErrDown(isnan(accEffErrDown)) = 0;
    
end

errorbar(hca,thBinEdges(1:end-1)+dthBin/2,accEffAvg*100,accEffErrDown*100,accEffErrUp*100,'-o','color',col1,'linewidth',2.2)
xErr = thBinEdges(1:end-1)+dthBin/2;

[smErrX,smUpY] = anjo.smooth_line(xErr,accEffErrUp);
[~,smDownY] = anjo.smooth_line(xErr,accEffErrDown);
[~,smAvgY] = anjo.smooth_line(xErr,accEffAvg);

%hfill = fill(hca,[xErr,fliplr(xErr)],([accEffAvg,fliplr(accEffAvg)]+[accEffErrUp,fliplr(accEffErrDown)])*100,col1);
hfill = fill(hca,[smErrX,fliplr(smErrX)],([smAvgY,fliplr(smAvgY)]+[smUpY,fliplr(smDownY)])*100,col1);
hfill.EdgeColor = 'none';
hfill.FaceAlpha = .5;
%uistack(hfill,'bottom');

hca.XLim = [0,90];
hca.YLim = [0,15];

uistack(hsc,'top');

plot(hca,45*[1,1],hca.YLim,'--','color',textcol,'linewidth',1.2)

hca.Box = 'on';

ylabel(hca,'Acceleration efficiency [$\%$]','Fontsize',15,'interpreter','latex')
xlabel(hca,'$\theta_{Bn}$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

title(hca,'Energy flux of ions with $E>10E_{sw}$ measured by MMS-FPI','Fontsize',15,'interpreter','latex','color',textcol)
hleg = irf_legend(hca,['$N = ',num2str(Nevents),'$'],[0.98,0.98],'Fontsize',15,'interpreter','latex','color',textcol);
irf_legend(hca,[num2str(beta*100),'\% significance'],[0.98,0.92],'Fontsize',15,'interpreter','latex','color',textcol);
hleg.BackgroundColor = hca.Color;

hca.LineWidth = 1.2;
hca.FontSize = 14;


%% Plot acceleration efficiency as a function of angle to sun-earth line
fig = figure;
hca = axes(fig);

% plot events with Ma as color
scatter(hca,phiV,accEffV*100,400,MaV.*cosd(thVnV),'.')
hold(hca,'on')
hca.XLim = [0,90];
hca.YLim(1) = 0;

sh_cmap(hca,cmap)
hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
%grid(hca,'on')

hcb = colorbar(hca);
hcb.Color = textcol;
hca.CLim = [0,20];

hca.Box = 'on';

ylabel(hca,'Acceleration efficiency [$\%$]','Fontsize',15,'interpreter','latex')
xlabel(hca,'$\alpha$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')
ylabel(hcb,'$M_A$','Fontsize',15,'interpreter','latex')

title(hca,'Energy flux of ions with $E>10E_{sw}$ measured by MMS-FPI','Fontsize',15,'interpreter','latex','color',textcol)
hleg = irf_legend(hca,['$N = ',num2str(Nevents),'$'],[0.98,0.98],'Fontsize',15,'interpreter','latex','color',textcol);
hleg.BackgroundColor = hca.Color;

hca.LineWidth = 1.2;
hca.FontSize = 14;
hcb.LineWidth = 1.2;



%% plot time

fig = figure;
hca = axes(fig);

irf_plot([TV,alphaV],'.','markersize',20)
hold(hca,'on')
plot(hca.XLim,[0,0],'k--','linewidth',1.8)
plot(hca.XLim,45*[1,1],'--','color',[1,1,1]*.4,'linewidth',1.2)
plot(hca.XLim,-45*[1,1],'--','color',[1,1,1]*.4,'linewidth',1.2)
plot(hca.XLim,90*[1,1],'--','color',[1,1,1]*.7,'linewidth',1.2)
plot(hca.XLim,-90*[1,1],'--','color',[1,1,1]*.7,'linewidth',1.2)

grid(hca,'off')
hca.YLim = [-1,1]*110;

xlabel(hca,'Time','Fontsize',15,'interpreter','latex')
ylabel(hca,'$\alpha$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')

%% Histogram of thetaBn

bincol = [157,214,166]/255;

fig = figure;
hca = axes(fig);

% thBnEdges is defined above
histogram(hca,thBnV,thBinEdges,'FaceColor',bincol,'edgecolor',textcol);

xlabel(hca,'$\theta_{Bn}$','Fontsize',15,'interpreter','latex')
ylabel(hca,'Number of events','Fontsize',15,'interpreter','latex')

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

hca.LineWidth = 1.2;
hca.FontSize = 14;

%% Histogram of Ma

bincol = [157,214,166]/255;

fig = figure;
hca = axes(fig);

dMaBin = 2;
MaBinEdges = 0:dMaBin:60;
histogram(hca,MaV.*cosd(thVnV),MaBinEdges,'FaceColor',bincol,'edgecolor',textcol);

xlabel(hca,'$M_A$','Fontsize',15,'interpreter','latex')
ylabel(hca,'Number of events','Fontsize',15,'interpreter','latex')

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

hca.LineWidth = 1.2;
hca.FontSize = 14;




