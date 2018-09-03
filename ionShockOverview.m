
% Plot the ion overview
% this script is called from within plotShockOverview


%% EIS average spectrogram structure
if plotEIS 
    % time step of eis detectors
    dTEis = median(diff(EISdpf0.time.epochUnix));
    
    EISspecrec = [];
    EISspecrec.p = EISpsd.data;
    EISspecrec.t = EISpsd.time.epochUnix+dTEis/2;
    EISspecrec.f = Eeis;
    EISspecrec.f_label = 'E (eV)';
end

%% Figure
h = sh_figure(6,[12,16]);

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

% n
hca = irf_panel(h,'n');
irf_plot(hca,ni)
hold(hca,'on')
irf_plot(hca,ne)
ylabel(hca,'$N$ [cm$^{-3}$]','fontsize',15,'interpreter','latex')
hca.YLim(1) = 0;
irf_legend(hca,{'$N_i$';'$N_e$'},[1.02,0.9],'Fontsize',15,'interpreter','latex')

% Vi
hca = irf_panel(h,'Vi');
irf_plot(hca,Vi)
ylabel(hca,'$\mathbf{V}_i$ [km/s]','fontsize',15,'interpreter','latex')
irf_legend(hca,{'$V_x$';'$V_y$';'$V_z$'},[1.02,0.9],'Fontsize',15,'interpreter','latex')

% fi (FPI and EIS)
hca = irf_panel(h,'fi');
% plot EIS data if read
if plotEIS
    hold(hca,'on')
    irf_spectrogram(hca,EISspecrec,'donotshowcolorbar');
    %hca.YLim(2) = max(EISspecrec.f);
    hca.YLimMode = 'auto';
end
% plot FPI data
iPDistSI = iPDist;
iPDistSI.data = iPDist.data*1e12;
irf_spectrogram(hca,iPDistSI.omni.specrec,'donotshowcolorbar')
hca.YScale = 'log';
sh_cmap(hca,'irf')
ylabel(hca,'Energy [eV]','fontsize',15,'interpreter','latex')
hca.YTick = 10.^[1,2,3,4];
idPanel = find(hca==h);
hcb1 = sh_cbar(h(idPanel));
ylabel(hcb1,{'$\log{f_i}$ ';'[s$^3$\,m$^{-6}$]'},'fontsize',14,'interpreter','latex')
% dividing line
maxETS = irf.ts_scalar(iPDist.time,iPDist.depend{1}(:,end));
irf_plot(hca,maxETS,'k-','linewidth',0.7)
grid(hca,'off')
hca.YLim(1) = 10; % 10 eV limit in plot

% reduced fi
hca = irf_panel(h,'redi');
irf_spectrogram(hca,f1Dn.specrec('1D_velocity'),'donotshowcolorbar')
sh_cmap(hca,'irf')
ylabel(hca,'$V_n$ [km/s]','fontsize',15,'interpreter','latex')
hcb2 = colorbar(hca);
ylabel(hcb2,{'$\log{F_i}$ ';'[s\,m$^{-4}$]'},'fontsize',15,'interpreter','latex')


irf_zoom(h,'x',tint)
irf_plot_axis_align(h)
sh_panel_span(h,[axu,axu+axh])
pause(0.01)

if plotEIS
    % make ion energy panel large
    dhf = 0.3; % how much larger?
    dhh = dhf*h(idPanel).Position(4);
    h(idPanel).Position(4) = h(idPanel).Position(4)+dhh;
    sh_panel_span(h(1:idPanel-1),[h(idPanel-1).Position(2)+dhh*(idPanel-1)/(length(h)-1),axu+axh])
    h(idPanel).Position(2) = h(idPanel-1).Position(2)-h(idPanel).Position(4);
    sh_panel_span(h(idPanel+1:end),[axu,h(idPanel).Position(2)])
end

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
%hca = irf_panel(h,'fi');
    hcb1.Position([2,4]) = h(idPanel).Position([2,4]);
hcb1.Position([1,3]) = [cbl,cbw];
hca = irf_panel(h,'redi');
hcb2.Position([2,4]) = hca.Position([2,4]);
hcb2.Position([1,3]) = [cbl,cbw];