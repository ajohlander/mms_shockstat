% Plot MMS position at shock crossings


%% Initiate figure
figure;
h = gobjects(1,3);

%% xy
hca = subplot(1,3,1:2);

scatter(hca,RV(:,1),RV(:,2),200,thBnV,'.')
axis(hca,'equal')
hold(hca,'on')

hcb = colorbar(hca,'location','northoutside');

% Plot Earth
theta=0:pi/20:pi;
xEarth=sin(theta);yEarth=cos(theta);
patch(-xEarth,yEarth,'k','edgecolor','k','Parent',hca)
patch(xEarth,yEarth,'w','edgecolor','k','Parent',hca)
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
scatter(hca,RV(:,1),RV(:,3),200,thBnV,'.')
axis(hca,'equal')
hold(hca,'on')

% Plot Earth
theta=0:pi/20:pi;
xEarth=sin(theta);yEarth=cos(theta);
patch(-xEarth,yEarth,'k','edgecolor','k','Parent',hca)
patch(xEarth,yEarth,'w','edgecolor','k','Parent',hca)
plot(hca,[-xEarth,xEarth],[yEarth,yEarth],'k','linewidth',1.5)

hca.XDir = 'reverse';

xlabel(hca,'$X_{GSE}$ [$R_E$]','fontsize',15,'interpreter','latex')
ylabel(hca,'$Z_{GSE}$ [$R_E$]','fontsize',15,'interpreter','latex')
ylabel(hcb,'$\theta_{Bn}$ [$^{\circ}$]','fontsize',15,'interpreter','latex')

h(2) = hca;

%
%% yz
hca = subplot(2,3,6);
scatter(hca,RV(:,2),RV(:,3),200,thBnV,'.')
axis(hca,'equal')
hold(hca,'on')

% Plot Earth
theta=0:pi/20:pi;
xEarth=sin(theta);yEarth=cos(theta);
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
end



%% plot angles

fig = figure;
hca = axes(fig);

scatter(atand(RV(:,2)./RV(:,1)),thBnV,200,MaV.*cosd(thVnV),'.')
hcb = colorbar(hca);

xlabel(hca,'$\alpha$ [$^{\circ}$]','fontsize',15,'interpreter','latex')
ylabel(hca,'$\theta_{Bn}$ [$^{\circ}$]','fontsize',15,'interpreter','latex')
ylabel(hcb,'$M_A$','fontsize',15,'interpreter','latex')

hca.Box = 'on';
hca.LineWidth = 1.2;
hcb.LineWidth = 1.2;
hca.FontSize = 14;
hca.CLim(1) = 0;
