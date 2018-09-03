% Plot MMS position at shock crossings


%% Initiate figure
fig = figure;
h = gobjects(1,3);

%% xy
hca = subplot(1,3,1:2);

scatter(hca,RV(~isnan(MaV),1),RV(~isnan(MaV),2),200,thBnV(~isnan(MaV)),'.')
axis(hca,'equal')
hold(hca,'on')

hcb = colorbar(hca,'location','northoutside');

% Plot Earth
theta=0:pi/20:pi;
xEarth=sin(theta);yEarth=cos(theta);
patch(-xEarth,yEarth,'k','edgecolor','k','Parent',hca)
patch(xEarth,yEarth,daycol,'edgecolor','k','Parent',hca)
plot(hca,[-xEarth,xEarth],[yEarth,yEarth],'k','linewidth',1.5)

hca.XDir = 'reverse';
hca.YDir = 'reverse';

xlabel(hca,'$X_{GSE}$ [$R_E$]','fontsize',15,'interpreter','latex')
ylabel(hca,'$Y_{GSE}$ [$R_E$]','fontsize',15,'interpreter','latex')
ylabel(hcb,'$\theta_{Bn}$ [$^{\circ}$]','fontsize',15,'interpreter','latex')

hcb.LineWidth = 1.2;

h(1) = hca;

%% xz
hca = subplot(2,3,3);
scatter(hca,RV(~isnan(MaV),1),RV(~isnan(MaV),3),200,thBnV(~isnan(MaV)),'.')
axis(hca,'equal')
hold(hca,'on')

% Plot Earth
theta=0:pi/20:pi;
xEarth=sin(theta);yEarth=cos(theta);
patch(-xEarth,yEarth,'k','edgecolor','k','Parent',hca)
patch(xEarth,yEarth,textcol,'edgecolor','k','Parent',hca)
plot(hca,[-xEarth,xEarth],[yEarth,yEarth],'k','linewidth',1.5)

hca.XDir = 'reverse';

xlabel(hca,'$X_{GSE}$ [$R_E$]','fontsize',15,'interpreter','latex')
ylabel(hca,'$Z_{GSE}$ [$R_E$]','fontsize',15,'interpreter','latex')
ylabel(hcb,'$\theta_{Bn}$ [$^{\circ}$]','fontsize',15,'interpreter','latex')

h(2) = hca;

%
%% yz
hca = subplot(2,3,6);
scatter(hca,RV(~isnan(MaV),2),RV(~isnan(MaV),3),200,thBnV(~isnan(MaV)),'.')
axis(hca,'equal')
hold(hca,'on')

% Plot Earth
theta=0:pi/20:pi;
xEarth=sin(theta);yEarth=cos(theta);
patch([-xEarth,xEarth],[yEarth,fliplr(yEarth)],textcol,'edgecolor','k','Parent',hca)
plot(hca,[-xEarth,xEarth],[yEarth,fliplr(yEarth)],'k','linewidth',1.5)

xlabel(hca,'$Y_{GSE}$ [$R_E$]','fontsize',15,'interpreter','latex')
ylabel(hca,'$Z_{GSE}$ [$R_E$]','fontsize',15,'interpreter','latex')
ylabel(hcb,'$\theta_{Bn}$ [$^{\circ}$]','fontsize',15,'interpreter','latex')

h(3) = hca;


%% plot stuff
for ii = 1:length(h)
    h(ii).Box = 'on';
    h(ii).LineWidth = 1.2;
    h(ii).FontSize = 14;
    
    h(ii).CLim = [0,90];
    sh_cmap(h(ii),cmap)
    h(ii).Color=axcol;
    
    h(ii).XAxis.Color = textcol;
    h(ii).YAxis.Color = textcol;
    h(ii).FontSize = 15;
end
fig.Color=figcol;
hcb.Color = textcol;
hcb.FontSize = 15;
hcb.Ticks = 0:15:90;


%% plot bow shock sigma

fig = figure;
hca = axes(fig);

scatter(alphaV,MaV,300,sigV,'.')
%hca.YLim = [0,90];

hca.Color = [1,1,1]*.5;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;

hcb = colorbar(hca);
sh_cmap(hca,'bluered')
hcb.Color = textcol;

title(hca,'Bow shock model scaling factor','Fontsize',15,'interpreter','latex','color',textcol)

xlabel(hca,'$\alpha$ [$^{\circ}$]','fontsize',15,'interpreter','latex')
ylabel(hca,'$\theta_{Bn}$ [$^{\circ}$]','fontsize',15,'interpreter','latex')
ylabel(hcb,'$\sigma$','fontsize',17,'interpreter','latex')

hca.Box = 'on';
hca.LineWidth = 1.2;
hcb.LineWidth = 1.2;
hca.FontSize = 14;
hca.CLim = 1+[-1,1]*diff(hca.CLim)/2;
