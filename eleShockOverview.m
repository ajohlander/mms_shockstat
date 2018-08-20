
% Plot the electron overview
% this script is called from within plotShockOverview

h = sh_figure(7,[12,16]);

% babs 4sc
hca = irf_panel(h,'babs');
mcol = [[0 0 0]; [213,94,0]/255;[0,158,115]/255;[86,180,233]/255;[0 1 1]];
hca.ColorOrder = mcol;
hold(hca,'on')
c_eval('irf_plot(hca,B?.abs,''color'',mcol(?,:));')
ylabel(hca,'$B$ [nT]','fontsize',15,'interpreter','latex')
irf_legend(hca,{'MMS1';'MMS2';'MMS3';'MMS4'},[1.02,0.9],'color','mms','Fontsize',15,'interpreter','latex')
hca.YLimMode ='auto';
hca.YLim(1) = 0;

% B 1sc
hca = irf_panel(h,'bxyz');
irf_plot(hca,B)
ylabel(hca,'$\mathbf{B}$ [nT]','fontsize',15,'interpreter','latex')
irf_legend(hca,{'$B_x$';'$B_y$';'$B_z$'},[1.02,0.9],'Fontsize',15,'interpreter','latex')
hleg = irf_legend(hca,['MMS',num2str(ic)],[0.02,0.95],'Fontsize',15,'interpreter','latex');
hleg.BackgroundColor = 'w';

% E
hca = irf_panel(h,'E');
hl = irf_plot(hca,E);
uistack(hl(2),'bottom')
uistack(hl(3),'bottom')
ylabel(hca,'$\mathbf{E}$ [mV/m]','fontsize',15,'interpreter','latex')
irf_legend(hca,{'$E_x$';'$E_y$';'$E_z$'},[1.02,0.9],'Fontsize',15,'interpreter','latex')

% n
hca = irf_panel(h,'n');
irf_plot(hca,ni)
hold(hca,'on')
irf_plot(hca,ne)
ylabel(hca,'$N$ [cm$^{-3}$]','fontsize',15,'interpreter','latex')
hca.YLim(1) = 0;
irf_legend(hca,{'$N_i$';'$N_e$'},[1.02,0.9],'Fontsize',15,'interpreter','latex')

% Ve
hca = irf_panel(h,'Ve');
irf_plot(hca,Ve)
ylabel(hca,'$\mathbf{V}_e$ [km/s]','fontsize',15,'interpreter','latex')
irf_legend(hca,{'$V_x$';'$V_y$';'$V_z$'},[1.02,0.9],'Fontsize',15,'interpreter','latex')

% Te
hca = irf_panel(h,'Te');
irf_plot(hca,TePar)
hold(hca,'on')
irf_plot(hca,TePerp)
ylabel(hca,'$T_e$ [eV]','fontsize',15,'interpreter','latex')
hca.YLim(1) = 0;
irf_legend(hca,{'$T_{e,\parallel}$';'$T_{e,\perp}$'},[1.02,0.9],'Fontsize',15,'interpreter','latex')

% fi
hca = irf_panel(h,'fe');
ePDistSI = ePDist;
ePDistSI.data = ePDist.data*1e12;
irf_spectrogram(hca,ePDistSI.omni.specrec,'donotshowcolorbar')
hca.YScale = 'log';
sh_cmap(hca,'irf')
ylabel(hca,'Energy [eV]','fontsize',15,'interpreter','latex')
hca.YTick = 10.^[1,2,3,4];
hcb1 = colorbar(hca);
ylabel(hcb1,{'$\log{f_i}$ ';'[s$^3$\,m$^{-6}$]'},'fontsize',14,'interpreter','latex')


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

hcb1.LineWidth = 1.3; hcb2.LineWidth = 1.3;
hca = irf_panel(h,'fi');
hcb1.Position([2,4]) = hca.Position([2,4]);
hcb1.Position([1,3]) = [cbl,cbw];
hca = irf_panel(h,'redi');
hcb2.Position([2,4]) = hca.Position([2,4]);
hcb2.Position([1,3]) = [cbl,cbw];