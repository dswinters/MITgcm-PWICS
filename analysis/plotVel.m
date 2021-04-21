function plotVel(theta,kTopo)

om = 1.36*1e-4; % forcing frequency
fs = 8; fn = 'times';
thetaPrefix = sprintf('theta%3.2f_',theta); % File prefix for theta
kTopoPrefix = sprintf('kTopo%.8f_',kTopo); % File prefix for kTopo
rname = sprintf('run_%s%s',thetaPrefix,kTopoPrefix(1:end-1));
froot = fullfile('..','runs',rname);
fig_dir = sprintf('plotVel/theta%3.2f_kTopo%.8f',theta,kTopo);
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

% Load grid and time info
gridm = rdmnc(fullfile(froot, 'grid*'));
datt = rdmnc(fullfile(froot,'outs_sn.*'),'T','iter');
files = dir(fullfile(froot,'outs_sn.*.nc'));         % all files
fids = extractBetween({files.name},'outs_sn.','.t'); % time identifiers (ignore tile suffixes)
fids = unique(fids);
yf = gridm.Y(find(gridm.Depth(1,:)==0,1,'last')) + [0 100e3];

% Load corrugation parameters (contains location of flux line)
load(fullfile(froot,'corrugation_params.mat'),'xSin1');
[~,nx] = min(abs(gridm.Xp1-xSin1));

figure('position',[131 68 1020 560])
cm = cmocean('balance',101);

v0 = 0.1; % forcing amplitude
uclim = 0.1*v0;
vclim = 0.1*v0;


for i = 1:length(datt.iter)
    dat = load_data(froot,fids,datt.iter(i),{'UVEL','VVEL'});


    % First plot: u plan-view
    clf
    subplot(221)
    pcolor(gridm.Xp1,gridm.Y,squeeze(dat.UVEL(:,:,1))'); hold on
    title('u (m/s)')
    caxis([-1 1]*uclim);
    %
    colorbar
    shading flat
    colormap(gca,cm)
    plot(xSin1*[1 1]+70e3,yf,'k--')
    xlabel('x (km)')
    ylabel('y (km)')
    ylim([0 1e6]);
    xticklabels(xticks/1e3);
    yticklabels(yticks/1e3);
    contour(gridm.XC',gridm.YC',gridm.Depth',[0:500:4e3],'k');
    contour(gridm.XC',gridm.YC',gridm.Depth',[0 0],'color','k','linewidth',2);

    % Second plot: v plan-view
    subplot(222)
    pcolor(gridm.X,gridm.Yp1,squeeze(dat.VVEL(:,:,1))'); hold on
    title('v (m/s)')
    caxis([-1 1]*vclim);
    %
    colorbar
    shading flat
    colormap(gca,cm)
    plot(xSin1*[1 1]+70e3,yf,'k--')
    xlabel('x (km)')
    ylabel('y (km)')
    ylim([0 1e6]);
    xticklabels(xticks/1e3);
    yticklabels(yticks/1e3);
    contour(gridm.XC',gridm.YC',gridm.Depth',[0:500:4e3],'k');
    contour(gridm.XC',gridm.YC',gridm.Depth',[0 0],'color','k','linewidth',2);

    % Third plot: u slice at flux calculation line
    subplot(223)
    pcolor(gridm.Y,gridm.Z,squeeze(dat.UVEL(nx,:,:))'); hold on
    title('u (m/s)')
    caxis([-1 1]*uclim);
    %
    plot(yf(end)*[1 1],ylim,'k--')
    colorbar
    shading flat
    colormap(gca,cm)
    plot(gridm.Y,-gridm.Depth(end,:),'k-');
    xlim([0 500e3]);
    xlabel('y (km)')
    ylabel('z (km)')
    xticklabels(xticks/1000)
    yticklabels(yticks/1000)

    % Fourth plot: v slice at flux calculation zone
    subplot(224)
    pcolor(gridm.Yp1,gridm.Z,squeeze(dat.VVEL(nx,:,:))'); hold on
    title('v (m/s)')
    caxis([-1 1]*vclim);
    %
    plot(yf(end)*[1 1],ylim,'k--')
    colorbar
    shading flat
    colormap(gca,cm)
    plot(gridm.Y,-gridm.Depth(end,:),'k-');
    xlim([0 500e3]);
    xlabel('y (km)')
    ylabel('z (km)')
    xticklabels(xticks/1000)
    yticklabels(yticks/1000)

    ttxt = sprintf('\\theta=%.1f^\\circ, kTopo=%.2f | T=%.1f cycles | [nx ny nz]=[%d %d %d]',...
                   theta,kTopo,datt.T(i)/(2*pi/om),...
                   length(gridm.X),length(gridm.Y),length(gridm.Z));
    hax = axes('visible','off','position',[0 0 1 1]);
    xlim(hax,[-1 1]); ylim(hax,[-1 1]);
    text(hax,0,1,ttxt,...
         'fontweight','bold',...
         'verticalalignment','top',...
         'horizontalalignment','center');

    fout=fullfile(fig_dir,sprintf('vel_%04d',i));
    print('-djpeg90','-r300',fout)
    fprintf('\rSaved %s [%0.2f%%]',fout,100 * i / length(datt.iter));

end




function [dat] = load_data(froot,fids,iter,vars)

    % If the simulation is too large for a single .nc file, rdmnc fails when
    % trying to load "outs_sn.*.nc". This workaround keeps trying each file
    % until the one containing the desired iteration is found.
    % TODO: There's probably a better way to do this.
    fileFound = false;
    nf = 0;
    while ~fileFound & nf < length(fids)
        nf = nf + 1; % try next file
        try
            dat = rdmnc(fullfile(froot,['outs_sn.' fids{nf},'.t*.nc']),vars{:},iter);
            fileFound = true;
        catch err
        end
    end
    if ~fileFound
        error('Failed loading data from iter %d',iter)
    end
