% ion composition plot
% this script is called from within plotShockOverview


%% create omni spectrogram structures
hpomnispec = [];
hpomnispec.f = hpca_energy;
hpomnispec.t = hpDist.time.epochUnix;
hpomnispec.p = squeeze(mean(hpDist.data*1e12,2));
hpomnispec.f_label = 'E_i (eV)';

he2pomnispec = hpomnispec;
he2pomnispec.p = squeeze(mean(he2pDist.data*1e12,2));

opomnispec = hpomnispec;
opomnispec.p = squeeze(mean(opDist.data*1e12,2));


%% create reduced velocity spectrogram structures
hpspecx = [];
hpspecx.t = hpred.t.epochUnix+5;
hpspecx.f = hpred.v;
hpspecx.p = hpred.fx;
hpspecx.f_label = 'V_n (km/s)';

he2pspecx = [];
he2pspecx.t = he2pred.t.epochUnix+5;
he2pspecx.f = he2pred.v;
he2pspecx.p = he2pred.fx;
hpspecx.f_label = 'V_n (km/s)';



%% Plot

h = sh_figure(9,[12,16]);

% B 1sc
hca = irf_panel(h,'bxyz');
irf_plot(hca,B)
hold(hca,'on')
irf_plot(hca,B.abs)
ylabel(hca,'$\mathbf{B}$ [nT]','fontsize',15,'interpreter','latex')
irf_legend(hca,{'$B_x$';'$B_y$';'$B_z$';'$|B|$'},[1.02,0.98],'Fontsize',15,'interpreter','latex')
hleg = irf_legend(hca,['MMS',num2str(ic)],[0.02,0.95],'Fontsize',15,'interpreter','latex');
hleg.BackgroundColor = 'w';

% n
hca = irf_panel(h,'n');
irf_plot(hca,ni)
hold(hca,'on')
irf_plot(hca,ne)
ylabel(hca,'$N$ [cm$^{-3}$]','fontsize',15,'interpreter','latex')
hca.YLim(1) = 0;
irf_legend(hca,{'$N_i$';'$N_e$'},[1.02,0.9],'Fontsize',15,'interpreter','latex')

% fi fpi
hca = irf_panel(h,'fi');
iPDistSI = iPDist;
iPDistSI.data = iPDist.data*1e12;
irf_spectrogram(hca,iPDistSI.omni.specrec,'donotshowcolorbar')
hca.YScale = 'log';
sh_cmap(hca,'irf')
ylabel(hca,'E [eV]','fontsize',15,'interpreter','latex')
hca.YTick = 10.^[1,2,3,4];

% fi hp
hca = irf_panel(h,'hpdist');
irf_spectrogram(hca,hpomnispec,'donotshowcolorbar')
hca.YScale = 'log';
sh_cmap(hca,'irf')
ylabel(hca,'E [eV]','fontsize',15,'interpreter','latex')
hca.YTick = 10.^[1,2,3,4];
irf_legend(hca,'$\mathrm{H}^{+}$',[0.02,0.95],'Fontsize',15,'interpreter','latex');

% fi he2p
hca = irf_panel(h,'he2pdist');
irf_spectrogram(hca,he2pomnispec,'donotshowcolorbar')
hca.YScale = 'log';
sh_cmap(hca,'irf')
ylabel(hca,'E [eV]','fontsize',15,'interpreter','latex')
hca.YTick = 10.^[1,2,3,4];
irf_legend(hca,'$\mathrm{He}^{2+}$',[0.02,0.95],'Fontsize',15,'interpreter','latex');

% fi op
hca = irf_panel(h,'opdist');
irf_spectrogram(hca,opomnispec,'donotshowcolorbar')
hca.YScale = 'log';
sh_cmap(hca,'irf')
ylabel(hca,'E [eV]','fontsize',15,'interpreter','latex')
hca.YTick = 10.^[1,2,3,4];
irf_legend(hca,'$\mathrm{O}^{+}$',[0.02,0.95],'Fontsize',15,'interpreter','latex');

% reduced fi
hca = irf_panel(h,'redi');
irf_spectrogram(hca,f1Dn.specrec('1D_velocity'),'donotshowcolorbar')
sh_cmap(hca,'irf')
ylabel(hca,'$V_n$ [km/s]','fontsize',15,'interpreter','latex')

% reduced hp dist
hca = irf_panel(h,'hpx');
irf_spectrogram(hca,hpspecx);
ylabel(hca,'$V_{n}$ [km/s]','interpreter','latex')
sh_cmap(hca,'irf')
irf_legend(hca,'$\mathrm{H}^{+}$',[0.02,0.95],'Fontsize',15,'interpreter','latex');

% reduced he2p dist
hca = irf_panel(h,'he2px');
irf_spectrogram(hca,he2pspecx);
ylabel(hca,'$V_{n}$ [km/s]','interpreter','latex')
sh_cmap(hca,'irf')
irf_legend(hca,'$\mathrm{He}^{2+}$',[0.02,0.95],'Fontsize',15,'interpreter','latex');




%% more stuff

irf_zoom(h,'x',tint)
irf_plot_axis_align(h)
sh_panel_span(h,[axu,axu+axh])
pause(0.01)
for jj = 1:length(h)
    irf_zoom(h(jj),'y',h(jj).YLim)
    h(jj).Position(1) = axl;
    h(jj).Position(3) = axw;
    
    h(jj).Layer = 'top';
    h(jj).FontSize = 15;
    h(jj).LineWidth = 1.3;
    h(jj).YLabel.Position(1) = -0.12;
end


%hca = irf_panel(h,'fi');
% hcb1.Position([2,4]) = hca.Position([2,4]);
% hcb1.Position([1,3]) = [cbl,cbw];
% hca = irf_panel(h,'redi');
% hcb2.Position([2,4]) = hca.Position([2,4]);
% hcb2.Position([1,3]) = [cbl,cbw];


%% colorbars
hcb1 = sh_cbar(h(3:6));
ylabel(hcb1,{'$\log{f_i}$ ';'[s$^3$\,m$^{-6}$]'},'fontsize',14,'interpreter','latex')

hcb2 = sh_cbar(h(7:9));
ylabel(hcb2,{'$\log{F_i}$ ';'[s$^2$\,m$^{-5}$]'},'fontsize',15,'interpreter','latex')



hcb1.LineWidth = 1.3; hcb2.LineWidth = 1.3;

hcb1.Position(1) = 0.85;
hcb2.Position(1) = 0.85;

