%% read time intervals

rePlotAll = irf_ask('Redo all plots? (0:no, 1:yes) [%]>','rePlotAll',0);

filePath = irf_ask('Choose path for plots [%]>','filePath','../mms_shock_plots/ions/');

u = irf_units;

% normal bs model to use (cannot be slho)
shModel = 'farris';

% for reduced distributions
nMC = 2e2;
vlim = 800;
vg = linspace(-vlim,vlim,100);

% plot position parameters
axl = 0.125; % axis left
axu = 0.06;
axw = .71;
axh = .87;
cbl = axl+axw+.01;
cbw = .03;


%% Make plot
tline = 1;
lineNum = 0;

fid = fopen(listFileName);

% loop through skipped lines
for ii = 1:startLine-1
    lineNum = lineNum+1;
    tline = fgets(fid);
end

% super trooper looper
while tline ~= -1
    lineNum = lineNum+1;
    disp(['Current line number: ',num2str(lineNum)])
    
    %% end loop if requested
    if lineNum >= stopLine; disp('Reached stop line, exiting...'); break; end
    
    %% read line from file
    tline = fgets(fid);
    
    if tline(1) == -1 || ~strcmp(tline(1),'2')
        disp('done!')
        break;
    end
    tintStr = [tline(1:10),'T',tline(12:19),'/',tline(23:32),'T',tline(34:41)];
    
    % all other strings if needed
    restStr = tline(44:end);
    
    % find commas
    idcomma = strfind(restStr,',');
    % find start paranthesis
    idpar = strfind(restStr,'(');
    
    fomStr = restStr(1:idcomma(1)-1);
    sitlStr = restStr(idcomma(1)+2:idpar(1)-1);
    
    %%
    descStr = restStr(idcomma(2)+2:end);
    if length(descStr)>50
        descStr = descStr(1:50);
    end
    
    %% get time interval and check if already done
    tint = irf.tint(tintStr);
    
    tstr = tint(1).toUtc;
    tstr([5,8,14,17]) = '';
    tstr = tstr(1:15);
    
    listing = dir(filePath);
    fnCell = {listing.name};
    fileName = ['sh_ions_',tstr,'.png'];
    
    if ismember(fileName,fnCell) && ~rePlotAll
        disp(['Skipping ',tstr])
        continue;
    end
    
    
    
    %% get data
    c_eval('B? = mms.get_data(''B_gse_fgm_brst_l2'',tint,?);')
    c_eval('B = B?;',ic)
    
    if isempty(B)
        disp(['no FGM data for ',tstr,', skipping...'])
        continue;
    end
    
    ni = mms.get_data('Ni_fpi_brst_l2',tint,ic);
    ne = mms.get_data('Ne_fpi_brst_l2',tint,ic);
    Vi = mms.get_data('Vi_gse_fpi_brst_l2',tint,ic);
    
    try
        R = mms.get_data('R_gse',tint);
    catch
       disp(['error reading sc position for ',tstr,', skipping...'])
       continue;
    end
    

    try
        iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
        iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);
        % ignore flux where count is 1 (also makes function faster)
        iPDist.data(iPDist.data<1.1*iPDistErr.data) = 0;
    catch
       disp(['error reading PDist for ',tstr,', skipping...'])
       continue;
    end
    % check if all data is there
    if isempty(B1) || isempty(B2) || isempty(B3) || isempty(B4) || isempty(B) || ...
            isempty(ni) || isempty(ne) || isempty(Vi) || isempty(R) || ... % || isempty(E)
            isempty(iPDist)
        disp(['no data for ',tstr,', skipping...'])
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
    
    ff= irf_get_data(tint,'bx,by,bz,Ma,v,n','omni_min');
    
    if ~isempty(ff)
        Bu = mean(ff(:,2:4),1);
        if isnan(Bu); Bu = nan(1,3); end
        Ma = mean(ff(:,5));
        Vu = mean(ff(:,6));
        Nu = mean(ff(:,7));
        thBn = acosd(dot(Bu,nvec)/(norm(Bu)));
        if thBn>90 
            thBn = 180-thBn;
        end
    else
        disp('failed to read omni data, moving on with plot...')
        thBn = nan;
        Ma = nan;
        Bu = nan(1,3);
        Nu = nan;
        Vu = nan;
    end
    
    
    %% set title strings
    meanRx = mean(R.gseR1(:,1)+R.gseR2(:,1)+R.gseR3(:,1)+R.gseR4(:,1))/4/u.RE*1e3;
    meanRy = mean(R.gseR1(:,2)+R.gseR2(:,2)+R.gseR3(:,2)+R.gseR4(:,2))/4/u.RE*1e3;
    meanRz = mean(R.gseR1(:,3)+R.gseR2(:,3)+R.gseR3(:,3)+R.gseR4(:,3))/4/u.RE*1e3;
    Rstr = [num2str(round(meanRx,1)),',',num2str(round(meanRy,1)),',',num2str(round(meanRz,1))];
    nstr = [num2str(round(nvec(1),2)),',',num2str(round(nvec(2),2)),',',num2str(round(nvec(3),2))];
    
    Bustr = [num2str(round(Bu(1),2)),',',num2str(round(Bu(2),2)),',',num2str(round(Bu(3),2))];
    Vustr = num2str(round(Vu));
    Nustr = num2str(round(Nu));
    thstr = num2str(round(thBn));
    Mastr = num2str(round(Ma));
    
    if isnan(Bu(1)); Bustr = '-'; end
    if isnan(Vu); Vustr = '-'; end
    if isnan(Nu); Nustr = '-'; end
    if isnan(thBn); thstr = '-'; end
    if isnan(Ma); Mastr = '-'; end
    
    %% reduce fi
    f1Dn = iPDist.reduce('1D',nvec,'vg',vg,'nMC',nMC);
    
    %% plot
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

    % fi
	hca = irf_panel(h,'fi');
    iPDistSI = iPDist;
    iPDistSI.data = iPDist.data*1e12;
    irf_spectrogram(hca,iPDistSI.omni.specrec,'donotshowcolorbar')
    hca.YScale = 'log';
    sh_cmap(hca,'irf')
    ylabel(hca,'Energy [eV]','fontsize',15,'interpreter','latex')
    hca.YTick = 10.^[1,2,3,4];
    hcb1 = colorbar(hca);
    ylabel(hcb1,{'$\log{f_i}$ ';'[s$^3$\,m$^{-6}$]'},'fontsize',14,'interpreter','latex')

    
    % reduced fi
	hca = irf_panel(h,'redi');
    irf_spectrogram(hca,f1Dn.specrec('1D_velocity'),'donotshowcolorbar')
    sh_cmap(hca,'irf')
    ylabel(hca,'$V_n$ [km/s]','fontsize',15,'interpreter','latex')
    hcb2 = colorbar(hca);
    ylabel(hcb2,{'$\log{F_i}$ ';'[s$^2$\,m$^{-5}$]'},'fontsize',15,'interpreter','latex')

    

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
    
    %% titles
%     irf_legend(h(1),['$\mathbf{R} = [',Rstr,']\,R_E$'],[0.02,1.045],'Fontsize',15,'interpreter','latex')
%     irf_legend(h(1),['$\mathbf{\hat{n}} = [',nstr,']$'],[0.98,1.045],'Fontsize',15,'interpreter','latex')
    
    %% new title
    
    ll = -0.14; % legend left
    lr = 1.14; % right
    
    %descStr
    
    irf_legend(h(1),descStr(1:end-1),[ll,1.35],'Fontsize',15,'interpreter','latex','horizontalalignment','left')
    irf_legend(h(1),['FOM: ',fomStr],[lr,1.35],'Fontsize',15,'interpreter','latex','horizontalalignment','right')
    
    irf_legend(h(1),['$\mathbf{R} = [',Rstr,']\,R_E$'],[ll,1.2],'Fontsize',15,'interpreter','latex','horizontalalignment','left')
    irf_legend(h(1),['$\mathbf{\hat{n}} = [',nstr,']$'],[lr,1.2],'Fontsize',15,'interpreter','latex','horizontalalignment','right')
   
    irf_legend(h(1),['$\mathbf{B}_u = [',Bustr,']$\,nT,'],[ll,1.045],'Fontsize',14,'interpreter','latex','horizontalalignment','left')
    irf_legend(h(1),['$N_u = ',Nustr,'$\,cm$^{-3}$, $V_u = ',Vustr,'$\,km/s, $M_A = ',Mastr,'$, $\theta_{Bn} = ',thstr,'^{\circ}$'],[lr,1.045],'Fontsize',14,'interpreter','latex','horizontalalignment','right')
    %% save figure
    % give it the time stamp of the start time
    irf_print_fig(h(1),[filePath,fileName(1:end-4)],'png')
    
    %% close figure
    close(h(1).Parent)
    
end


