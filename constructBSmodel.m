% try to make supercool unofficial MMS SITL bow shock model


global RSC
RSC = RV;
% assume r0=0, alpha=3.5 deg

% alpha = 10*pi/180;
% RVabd = zeros(size(RV));
% 
% Q = [cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
% for ii = 1:size(RV,1)
%     RVabd(ii,:) = Q*RV(ii,:)';
% end
% 
% [theta,r] = cart2pol(RVabd(:,1),sqrt(RVabd(:,2).^2+RVabd(:,3).^2));

%Rfun = @(x)sum((x(1)./r-x(2)*cos(theta)-1).^2);

options = optimset('PlotFcns',@optimplotfval);
pars = fminsearch(@Rfun,[20,1,-0.01],options);


%%
L = pars(1);
eps = pars(2);
alpha = pars(3);


%% construct Rabd

Rabd = zeros(size(RV));
Q = [cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
for ii = 1:size(RV,1)
    Rabd(ii,:) = Q*RV(ii,:)';%-[x0;y0;z0];
end

%%  plot Rabd
fig = figure;
hca = axes(fig);


% plot model
N = 1000;
thModel = linspace(-3*pi/4,3*pi/4,N);
rModel = L./(1+eps*cos(thModel));
[xModel,yModel] = pol2cart(thModel,rModel);

plot(hca,xModel,yModel,'linewidth',4,'color',col2)
hold(hca,'on')

scatter(hca,Rabd(:,1),sqrt(Rabd(:,2).^2+Rabd(:,3).^2),200,col1,'.')
axis(hca,'equal')
%hca.YLim = [0,max(Rabd(:,2))+3];
%hca.XLim = [min(Rabd(:,1))-6,max(Rabd(:,1))+6];

hca.YLim = [0,25];
hca.XLim = [-10,25];


% Plot Earth
theta=0:pi/20:pi;
xEarth=sin(theta);yEarth=cos(theta);
patch(-xEarth,yEarth,'k','edgecolor','k','Parent',hca)
patch(xEarth,yEarth,daysidecol,'edgecolor','k','Parent',hca)
plot(hca,[-xEarth,xEarth],[yEarth,yEarth],'k','linewidth',1.5)

hca.XDir = 'reverse';


xlabel(hca,'$X''$ [$R_E$]','fontsize',15,'interpreter','latex')
ylabel(hca,'$(Y''^2+Z''^2)^{1/2}$ [$R_E$]','fontsize',15,'interpreter','latex')

legstr = {['$L = ',num2str(round(L,1)),'\,R_E$'];...
    ['$\epsilon = ',num2str(round(eps,2)),'$'];...
    ['$\alpha = ',num2str(round(alpha*180/pi,1)),'^{\circ}$']};
irf_legend(hca,legstr,[0.02,0.98],'fontsize',15,'interpreter','latex','color',textcol)

h(1) = hca;

hca.Box = 'on';
hca.LineWidth = 1.2;
hca.Color = axcol;
fig.Color = figcol;
hca.XAxis.Color = textcol;
hca.YAxis.Color = textcol;
hca.FontSize = 15;

%%
function Rval = Rfun(x)
global RSC

L = x(1);
eps = x(2);
alpha = x(3);
% x0 = x(4);
% y0 = x(5);
% z0 = 0;

Rabd = zeros(size(RSC));

Q = [cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
for ii = 1:size(RSC,1)
    Rabd(ii,:) = Q*RSC(ii,:)';%-[x0;y0;z0];
end

[theta,r] = cart2pol(Rabd(:,1),sqrt(Rabd(:,2).^2+Rabd(:,3).^2));

Rval = sum((L./r-eps*cos(theta)-1).^2);

end