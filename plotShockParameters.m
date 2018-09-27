
EmaxLim = irf_ask('Upper energy limit of FPI in units of solar wind energy (0 will ignore): [%]>','EmaxLim',0);

colorMode = irf_ask('Color mode (1:Fancy, 2:Boring): [%]>','colorMode',1);

%% colors for plots

switch colorMode
    case 1 % darker, cooler
        cmap = 'strangeways';
        axcol = [1,1,1]*.4;
        figcol = [1,1,1]*.2;
        textcol = [1,1,1]*.95;
        daycol = [1,1,1]*.95;
        col1 = [253,232,159]/255;
        col2 = [211,64,82]/255;
        col3 = [53,151,103]/255;
        bincol = [157,214,166]/255;
    case 2 % lighter, boring
        cmap = 'cool';
        axcol = [1,1,1];
        figcol = [1,1,1];
        textcol = [1,1,1]*.05;
        daycol = [1,1,1];
        col1 = [1,1,1]*.5;
        col2 = [211,64,82]/255;
        col3 = [53,151,103]/255;
        bincol = [157,214,166]/255;
end

colmat = [col1;col2;col3];



%% some bin edges
dthBin = 10;
thBinEdges = 0:dthBin:90;

% mach number
% dMBin = 4;
% MBinEdges = 1:dMBin:16;

MBinEdges = [1,6,10,100];

%% approved data indices

EswV = 1/2*u.mp*VuV.^2;
if EmaxLim ~= 0
    dind = (~isnan(MaV) & EmaxV>EswV*EmaxLim);
else
    dind = ~isnan(MaV);
end

%% Clean data arrays
dTV = dTV(dind);
MaV = MaV(dind);
MfV = MfV(dind);
VuV = VuV(dind);
thBnV = thBnV(dind);
thBrV = thBrV(dind);
thVnV = thVnV(dind);
betaiV = betaiV(dind);
accEffV = accEffV(dind);
accEffAltV = accEffAltV(dind);
accEffFpiV = accEffFpiV(dind);
EmaxV = EmaxV(dind);
RV = RV(dind,:);
sigV = sigV(dind);
hasEISV = hasEISV(dind);
% convert to boolean
hasEISV = (hasEISV==1);

dstV = dstV(dind);
kpV = kpV(dind);
ssnV = ssnV(dind);
s107V = s107V(dind);
aeV = aeV(dind);

lineNumV = lineNumV(dind);

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

scatter(hca,thBnV,MaV,400,col1,'.')

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


%% Plot acceleration efficiency as a function of shock angle
fig = figure;
hca = axes(fig);

% plot events with Ma as color
hold(hca,'on')
scatter(hca,thBnV(hasEISV),accEffV(hasEISV)*100,60,MaV(hasEISV),'o','MarkerFaceColor','flat')
scatter(hca,thBnV(~hasEISV),accEffV(~hasEISV)*100,60,MaV(~hasEISV),'o','MarkerFaceColor','flat')
%foo = scatter(hca,thBnV(~hasEISV),accEffV(~hasEISV)*100,60,MaV(~hasEISV),'x','linewidth',3);

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

%title(hca,'Energy flux of ions with $E>10E_{sw}$ measured by MMS-FPI','Fontsize',15,'interpreter','latex','color',textcol)
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

% upper and lower limits (difference to average, not actual limits)
accEffErrUp = zeros(1,length(thBinEdges)-1);
accEffErrDown = zeros(1,length(thBinEdges)-1);

for ii = 1:length(thBinEdges)-1
    % mean and std
    mu =  nanmean(accEffV(idTh==ii));
    sig = nanstd(accEffV(idTh==ii));
    
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
xErr = thBinEdges(1:end-1)+dthBin/2;

% make smooth lines of errors and average
smErrX = linspace(0,90,1e2);
smUpY = interp1(thBinEdges(1:end-1)+dthBin/2,accEffErrUp,smErrX,'pchip');
smDownY = interp1(thBinEdges(1:end-1)+dthBin/2,accEffErrDown,smErrX,'pchip');
smAvgY = interp1(thBinEdges(1:end-1)+dthBin/2,accEffAvg,smErrX,'pchip');

%hfill = fill(hca,[xErr,fliplr(xErr)],([accEffAvg,fliplr(accEffAvg)]+[accEffErrUp,fliplr(accEffErrDown)])*100,col1);
hfill = fill(hca,[smErrX,fliplr(smErrX)],([smAvgY,fliplr(smAvgY)]+[smUpY,fliplr(smDownY)])*100,col1);
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


%% 2D histogram of acceff of thBn-Ma/Mf parameter space
% sort of hopeful should fix

% lower resolution of thBn
dthBin2 = 18;
thBinEdges2 = 0:dthBin2:90;

idTh = discretize(thBnV,thBinEdges2);
idM = discretize(MaV,MBinEdges);

accEffAvg2D = zeros(length(thBinEdges2)-1,length(MBinEdges)-1);
accEffStd2D = zeros(length(thBinEdges2)-1,length(MBinEdges)-1);
accEffErrUp2D = zeros(length(thBinEdges2)-1,length(MBinEdges)-1);
accEffErrDown2D = zeros(length(thBinEdges2)-1,length(MBinEdges)-1);

for ii = 1:length(thBinEdges2)-1
    for jj = 1:length(MBinEdges)-1
        
        % mean and std
        mu =  nanmean(accEffV(idTh==ii & idM==jj));
        sig = nanstd(accEffV(idTh==ii & idM==jj));
        
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
    hsc = scatter(hca,thBnV(idM==jj),accEffV(idM==jj)*100,60,colmat(jj,:),'o','MarkerFaceColor','flat');
end


for jj = 1:size(accEffAvg2D,2)
    errorbar(hca,thBinEdges2(1:end-1)+dthBin2/2+xoffs(jj),accEffAvg2D(:,jj)*100,accEffErrDown2D(:,jj)*100,accEffErrUp2D(:,jj)*100,'linewidth',3,'color',colmat(jj,:))
end




%plot(hca,45*[1,1],hca.YLim,'--','color',textcol,'linewidth',1.2)


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

%% Plot acceleration efficiency as a function of angle to sun-earth line
fig = figure;
hca = axes(fig);

% plot events with Ma as color
scatter(hca,phiV,accEffV*100,400,MaV,'.')
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


%% Histogram of thetaBn

fig = figure;
hca = axes(fig);

% thBnEdges is defined above
histogram(hca,thBnV,thBinEdges,'FaceColor',bincol,'edgecolor',textcol,'linewidth',1.3);

xlabel(hca,'$\theta_{Bn}$','Fontsize',15,'interpreter','latex')
ylabel(hca,'Number of events','Fontsize',15,'interpreter','latex')

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

hca.LineWidth = 1.2;
hca.FontSize = 14;

%% Histogram of Ma

fig = figure;
hca = axes(fig);

dMaBin = 2;
MaBinEdges = 0:dMaBin:60;
histogram(hca,MaV,MaBinEdges,'FaceColor',bincol,'edgecolor',textcol,'linewidth',1.3);

xlabel(hca,'$M_A$','Fontsize',15,'interpreter','latex')
ylabel(hca,'Number of events','Fontsize',15,'interpreter','latex')

hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

hca.LineWidth = 1.2;
hca.FontSize = 14;


