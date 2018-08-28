
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


