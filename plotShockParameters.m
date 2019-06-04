

instrumentMode = irf_ask('Use EIS or only FPI? (1:FPI, 2:EIS): [%]>','instrumentMode',1);

EmaxFac = irf_ask('Emax must be at least this much greater than Esw: [%]>','EmaxFac',20);

% PLACEHOLDER, SHOULD BE INCLUDED IN SAVED FILE!
shModel = {'farris','slho','per','fa4o','fan4o','foun'};

fprintf('\n')
for jj = 1:length(shModel); fprintf([num2str(jj),': ',shModel{jj},'\n']); end
modelInd = irf_ask('Which shock model: [%]>','modelInd',1);

%% some bin edges
dthBin = 10;
thBinEdges = 0:dthBin:90;
deltathBin = median(diff(thBinEdges))*pi/180;

% mach number
% dMBin = 4;
% MBinEdges = 1:dMBin:16;

MBinEdges = [1,6,10,100];

%% approved data indices

% now theres an error in EfpiMaxV (should be corrected when data is rerun)
dind = find((EfpiMaxV)./(.5*u.mp*(VuV*1e3).^2/u.e)>EmaxFac);

% all
%dind = 1:N;

% redefine N
N = numel(find(dind));

%% Clean data arrays
% single values
TV = TV(dind);
dTV = dTV(dind);
VuV = VuV(dind);
BuV = BuV(dind,:);
NuV = NuV(dind,:);
thBrV = thBrV(dind);
betaiV = betaiV(dind);
accEffV = accEffV(dind);
accEffAltV = accEffAltV(dind);
accEffFpiV = accEffFpiV(dind);
EmaxV = EmaxV(dind);
EfpiMaxV = EfpiMaxV(dind);
RV = RV(dind,:);
hasEISV = hasEISV(dind);
% convert to boolean
hasEISV = (hasEISV==1);
dstV = dstV(dind);
kpV = kpV(dind);
ssnV = ssnV(dind);
s107V = s107V(dind);
aeV = aeV(dind);
lineNumV = lineNumV(dind);
swIdV = swIdV(dind);
imfIdV = imfIdV(dind);

for jj = 1:length(shModel)
    nvecV.(shModel{jj}) = nvecV.(shModel{jj})(dind,:);
    MaV.(shModel{jj}) = MaV.(shModel{jj})(dind);
    MfV.(shModel{jj}) = MfV.(shModel{jj})(dind);
    thBnV.(shModel{jj}) = thBnV.(shModel{jj})(dind);
    thVnV.(shModel{jj}) = thVnV.(shModel{jj})(dind);
    sigV.(shModel{jj}) = sigV.(shModel{jj})(dind);
end

Nevents = numel(find(dind));

% angle between earth-sun line and sc position in xy plane
[alphaV,~] = cart2pol(RV(:,1),RV(:,2),RV(:,3));
alphaV = alphaV*180/pi; % degrees

% angle between earth-sun line and sc position
[phiV,~] = cart2pol(RV(:,1),sqrt(RV(:,2).^2+RV(:,3).^2));
phiV = phiV*180/pi; % degrees

%% Which acceleration efficiency to be used?

accEffV1 = accEffFpiV;

%% Which shock model to be used?

thBnV1 = thBnV.(shModel{modelInd});
thVnV1 = thVnV.(shModel{modelInd});
MaV1 = MaV.(shModel{modelInd});
MfV1 = MfV.(shModel{modelInd});
nvecV1 = nvecV.(shModel{modelInd});
sigV1 = sigV.(shModel{modelInd});


%% Plot simple position
plotShockPos

%% compare shock models
plotModelComp

%% plot parameter space
fig = figure;
hca = axes(fig);

scatter(hca,thBnV1,MaV1,400,col1,'.')

hca.XLim = [0,90];
hca.YLim(1) = 0;

hca.Box = 'on';

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

ylabel(hca,'$M_A$','fontsize',15,'interpreter','latex')
xlabel(hca,'$\theta_{Bn}$ [$^{\circ}$]','fontsize',15,'interpreter','latex')

hca.LineWidth = 1.2;
hca.FontSize = 14;
fig.InvertHardcopy = 'off';


%% Plot acceleration efficiency as a function of shock angle
fig = figure;
hca = axes(fig);

% plot events with Ma as color
hold(hca,'on')
%scatter(hca,thBnV1(hasEISV),accEffV1(hasEISV)*100,60,MaV1(hasEISV),'o','MarkerFaceColor','flat')
%scatter(hca,thBnV1(~hasEISV),accEffV1(~hasEISV)*100,60,MaV1(~hasEISV),'x','linewidth',3);

scatter(hca,thBnV1,accEffV1*100,60,MaV1,'o','MarkerFaceColor','flat')


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
hca.CLim = [5,20];

hca.Box = 'on';

ylabel(hca,'Acceleration efficiency [$\%$]','Fontsize',15,'interpreter','latex')
xlabel(hca,'$\theta_{Bn}$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')
ylabel(hcb,'$M_A$','Fontsize',15,'interpreter','latex')

%title(hca,'Energy flux of ions with $E>10E_{sw}$ measured by MMS-FPI','Fontsize',15,'interpreter','latex','color',textcol)
hleg = irf_legend(hca,['$N = ',num2str(Nevents),'$'],[0.98,0.98],'Fontsize',15,'interpreter','latex','color',textcol);
hleg.BackgroundColor = hca.Color;

hca.LineWidth = 1.2;
hca.FontSize = 14;
hcb.LineWidth = 1.2;
fig.InvertHardcopy = 'off';

%% same as above but with avg with proper Baysian limits
fig = figure;
hca = axes(fig);
hold(hca,'on')
hsc = scatter(hca,thBnV1,accEffV1*100,200,col2,'.');

% set significance
beta = .90;

idTh = discretize(thBnV1,thBinEdges);
accEffAvg = zeros(1,length(thBinEdges)-1);
accEffStd = zeros(1,length(thBinEdges)-1);

% upper and lower limits (difference to average, not actual limits)
accEffErrUp = zeros(1,length(thBinEdges)-1);
accEffErrDown = zeros(1,length(thBinEdges)-1);

for ii = 1:length(thBinEdges)-1
    % mean and std
    mu =  nanmean(accEffV1(idTh==ii));
    sig = nanstd(accEffV1(idTh==ii));
    
    % after a lot(!) of math
    accEffErrUp(ii) = norminv(beta*normcdf(mu/sig)-normcdf(mu/sig)+1)*sig;
    accEffErrDown(ii) = norminv(1-beta*normcdf(mu/sig))*sig;
    
    accEffAvg(ii) = mu;
    accEffStd(ii) = sig;
    
    
end

% clean arrays
accEffErrUp(isnan(accEffErrUp)) = 0;
accEffErrDown(isnan(accEffErrDown)) = 0;

errorbar(hca,thBinEdges(1:end-1)+dthBin/2,accEffAvg*100,accEffErrDown*100,accEffErrUp*100,'-o','color',col1,'linewidth',2.2)
accEffAvg(isnan(accEffAvg)) = 0;
xErr = thBinEdges(1:end-1)+dthBin/2;

% make smooth lines of errors and average
smErrX = linspace(0,90,1e2);
smUpY = interp1(thBinEdges(1:end-1)+dthBin/2,accEffErrUp,smErrX,'pchip');
smDownY = interp1(thBinEdges(1:end-1)+dthBin/2,accEffErrDown,smErrX,'pchip');
smAvgY = interp1(thBinEdges(1:end-1)+dthBin/2,accEffAvg,smErrX,'pchip');

hfill = fill(hca,[xErr,fliplr(xErr)],([accEffAvg,fliplr(accEffAvg)]+[accEffErrUp,fliplr(accEffErrDown)])*100,col1);
% hfill = fill(hca,[smErrX,fliplr(smErrX)],([smAvgY,fliplr(smAvgY)]+[smUpY,fliplr(smDownY)])*100,col1);
hfill.EdgeColor = 'none';
hfill.FaceAlpha = .5;
%uistack(hfill,'bottom');

hca.XLim = [0,90];
%hca.YLim = [0,15];
hca.YLim(1) = 0;

uistack(hsc,'top');

plot(hca,45*[1,1],hca.YLim,'--','color',textcol,'linewidth',1.2)

hca.Box = 'on';

ylabel(hca,'Acceleration efficiency [$\%$]','Fontsize',15,'interpreter','latex')
xlabel(hca,'$\theta_{Bn}$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

%title(hca,'Energy flux of ions with $E>10E_{sw}$ measured by MMS-FPI','Fontsize',15,'interpreter','latex','color',textcol)
hleg = irf_legend(hca,['$N = ',num2str(Nevents),'$'],[0.98,0.98],'Fontsize',15,'interpreter','latex','color',textcol);
irf_legend(hca,[num2str(beta*100),'\% significance'],[0.98,0.92],'Fontsize',15,'interpreter','latex','color',textcol);
hleg.BackgroundColor = hca.Color;

hca.LineWidth = 1.2;
hca.FontSize = 14;
fig.InvertHardcopy = 'off';

%% Mach number dependence between 20 and 40 deg

fig = figure;
hca = axes(fig);
hold(hca,'on')
hsc = scatter(hca,MaV1(thBnV1>20 & thBnV1<40),accEffV1(thBnV1>20 & thBnV1<40)*100,200,col2,'.');
fig.InvertHardcopy = 'off';


%% scatter with errorbars of acceff of thBn-various parameter space

% lower resolution of thBn
dthBin2 = 18;
thBinEdges2 = 0:dthBin2:90;
% Ma
MaBinEdges2 = [0,6,10,100];
% Mf
MfBinEdges2 = [0,5,7,100];
% phi
phiBinEdges2 = [0,30,60,120];


% Ma
fig = figure;
hca = axes(fig);
hold(hca,'on')
% do the plot
plot_acceff_dep(hca,thBnV1,thBinEdges2,accEffV1,MaV1,MaBinEdges2,colmat,'M_A','','line')
% fix stuff
hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
fig.InvertHardcopy = 'off';

% Mf
fig = figure;
hca = axes(fig);
hold(hca,'on')
% do the plot
plot_acceff_dep(hca,thBnV1,thBinEdges2,accEffV1,MfV1,MfBinEdges2,colmat,'M_f','','line')
% fix stuff
hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
fig.InvertHardcopy = 'off';

% phi
fig = figure;
hca = axes(fig);
hold(hca,'on')
% do the plot
plot_acceff_dep(hca,thBnV1,thBinEdges2,accEffV1,phiV,[0,60,110],colmat,'\phi','$^{\circ}$','line')
% fix stuff
hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
fig.InvertHardcopy = 'off';


% thBr
fig = figure;
hca = axes(fig);
hold(hca,'on')
% do the plot
plot_acceff_dep(hca,thBnV1,thBinEdges2,accEffV1,thBrV,[0,45,90],colmat,'\theta_{Br}','$^{\circ}$','line')
% fix stuff
hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
fig.InvertHardcopy = 'off';

%%

% lower resolution of thBn
dthBin2 = 18;
thBinEdges2 = 0:dthBin2:90;

idTh = discretize(thBnV1,thBinEdges2);
idM = discretize(MaV1,MBinEdges);

accEffAvg2D = zeros(length(thBinEdges2)-1,length(MBinEdges)-1);
accEffStd2D = zeros(length(thBinEdges2)-1,length(MBinEdges)-1);
accEffErrUp2D = zeros(length(thBinEdges2)-1,length(MBinEdges)-1);
accEffErrDown2D = zeros(length(thBinEdges2)-1,length(MBinEdges)-1);

for ii = 1:length(thBinEdges2)-1
    for jj = 1:length(MBinEdges)-1
        
        % mean and std
        mu =  nanmean(accEffV1(idTh==ii & idM==jj));
        sig = nanstd(accEffV1(idTh==ii & idM==jj));
        
        % after a lot(!) of math
        accEffErrUp2D(ii,jj) = norminv(beta*normcdf(mu/sig)-normcdf(mu/sig)+1)*sig;
        accEffErrDown2D(ii,jj) = norminv(1-beta*normcdf(mu/sig))*sig;
        
        accEffAvg2D(ii,jj) = mu;
        accEffStd2D(ii,jj) = sig;
       
    end
end

% clean arrays
accEffErrUp(isnan(accEffErrUp)) = 0;
accEffErrDown(isnan(accEffErrDown)) = 0;

xoffs = [-1,0,1,2];

%
fig = figure;
hca = axes(fig);
hold(hca,'on')



for jj = 1:size(accEffAvg2D,2)
    hsc = scatter(hca,thBnV1(idM==jj),accEffV1(idM==jj)*100,60,colmat(jj,:),'o','MarkerFaceColor','flat');
end


for jj = 1:size(accEffAvg2D,2)
    errorbar(hca,thBinEdges2(1:end-1)+dthBin2/2+xoffs(jj),accEffAvg2D(:,jj)*100,accEffErrDown2D(:,jj)*100,accEffErrUp2D(:,jj)*100,'linewidth',3,'color',colmat(jj,:))
end


legStr = cell(1,size(accEffAvg2D,2));
for jj = 1:size(accEffAvg2D,2)
    if jj == 1
        legStr{jj} = ['$M_A{<}',num2str(MBinEdges(jj+1)),'$, \hspace{.5cm} $N = ',num2str(numel(find(idM==jj))),'$'];
    elseif jj == size(accEffAvg2D,2)
        legStr{jj} = ['$M_A{>}',num2str(MBinEdges(jj)),'$, \hspace{.5cm} $N = ',num2str(numel(find(idM==jj))),'$'];
    else
        legStr{jj} = ['$',num2str(MBinEdges(jj)),'{<}M_A{<}',num2str(MBinEdges(jj+1)),'$, \hspace{.5cm} $N = ',num2str(numel(find(idM==jj))),'$'];
        
    end
end

hl = legend(hca,legStr);
hl.Interpreter = 'latex';
hl.FontSize = 15;
hl.Color = figcol;
hl.TextColor = textcol;

hca.Box = 'on';

hca.XLim = [0,90];
%hca.YLim = [0,15];
hca.YLim(1) = 0;

ylabel(hca,'Acceleration efficiency [$\%$]','Fontsize',15,'interpreter','latex')
xlabel(hca,'$\theta_{Bn}$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

%irf_legend(hca,[num2str(beta*100),'\% significance'],[0.98,0.92],'Fontsize',15,'interpreter','latex','color',textcol);
%hl.BackgroundColor = hca.Color;

hca.LineWidth = 1.2;
hca.FontSize = 14;
fig.InvertHardcopy = 'off';

%% Plot acceleration efficiency as a function of angle to sun-earth line
fig = figure;
hca = axes(fig);

% plot events with Ma as color
scatter(hca,phiV,accEffV1*100,400,thBnV1,'.')
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
hca.CLim = [0,90];

hca.Box = 'on';

ylabel(hca,'Acceleration efficiency [$\%$]','Fontsize',15,'interpreter','latex')
xlabel(hca,'$\alpha$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')
ylabel(hcb,'$\theta_{Bn}$','Fontsize',15,'interpreter','latex')

title(hca,'Energy flux of ions with $E>10E_{sw}$ measured by MMS-FPI','Fontsize',15,'interpreter','latex','color',textcol)
hleg = irf_legend(hca,['$N = ',num2str(Nevents),'$'],[0.98,0.98],'Fontsize',15,'interpreter','latex','color',textcol);
hleg.BackgroundColor = hca.Color;

hca.LineWidth = 1.2;
hca.FontSize = 14;
hcb.LineWidth = 1.2;
fig.InvertHardcopy = 'off';


%% plot time

fig = figure;
hca = axes(fig);

irf_plot([TV,alphaV],'.','color',col1,'markersize',20)
hold(hca,'on')
plot(hca.XLim,[0,0],'k--','linewidth',1.8)
plot(hca.XLim,45*[1,1],'--','color',axcol*.3,'linewidth',1.2)
plot(hca.XLim,-45*[1,1],'--','color',axcol*.3,'linewidth',1.2)
plot(hca.XLim,90*[1,1],'--','color',axcol*.7,'linewidth',1.2)
plot(hca.XLim,-90*[1,1],'--','color',axcol*.7,'linewidth',1.2)

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

grid(hca,'off')
hca.YLim = [-1,1]*110;

xlabel(hca,'Time','Fontsize',15,'interpreter','latex')
ylabel(hca,'$\alpha$ [$^{\circ}$]','Fontsize',15,'interpreter','latex')
hca.FontSize = 14;
fig.InvertHardcopy = 'off';


%% Histogram of thetaBn

fig = figure;
hca = axes(fig);

% thBnEdges is defined above
thhist = histogram(hca,thBnV1,thBinEdges,'FaceColor',bincol,'edgecolor',textcol,'linewidth',1.3);

% "analytical"
than = linspace(0,90,1e3);
hold(hca,'on')
%plot(than,sind(than)*N*deltathBin,'linewidth',2,'color',col2)
plot(sort([thBinEdges,thBinEdges(1:end-1)]),[0,sind(sort([thBinEdges,thBinEdges(1:end-2)]+deltathBin*180/pi/2))]*N*deltathBin,'linewidth',2,'color',col2)

hca.XLim = [0,90];
hca.YLim = [0,1.1*max(thhist.Values)];

xlabel(hca,'$\theta_{Bn}$','Fontsize',15,'interpreter','latex')
ylabel(hca,'Number of events','Fontsize',15,'interpreter','latex')

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

hca.LineWidth = 1.2;
hca.FontSize = 14;
fig.InvertHardcopy = 'off';

%% Histogram of Ma & Mf

fig = figure;
hca = axes(fig);

dMaBin = 2;
MaBinEdges = 0:dMaBin:60;
Mahist = histogram(hca,MaV1,MaBinEdges,'FaceColor',bincol,'edgecolor',textcol,'linewidth',1.3);
hold(hca,'on')
Mfhist = histogram(hca,MfV1,MaBinEdges,'FaceColor',col2,'edgecolor',textcol,'linewidth',1.3);

xlabel(hca,'$M_A$','Fontsize',15,'interpreter','latex')
ylabel(hca,'Number of events','Fontsize',15,'interpreter','latex')

hca.XLim = [0,20];
hca.YLim = [0,1.1*max(Mfhist.Values)];


hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

hleg = legend(hca,'Alfv\''en Mach','Magnetosonic Mach');
hleg.Interpreter = 'latex';
hleg.FontSize = 15;

hca.LineWidth = 1.3;
hca.FontSize = 15;
fig.InvertHardcopy = 'off';
hca.Position = [.1,.1,.85,.85];


