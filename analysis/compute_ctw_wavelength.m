function compute_ctw_wavelength(theta,kTopo)

om = 1.36*1e-4; % forcing frequency
fs = 8; fn = 'times';
thetaPrefix = sprintf('theta%3.2f_',theta); % File prefix for theta
kTopoPrefix = sprintf('kTopo%.8f_',kTopo); % File prefix for kTopo
rname = sprintf('run_%s%s',thetaPrefix,kTopoPrefix(1:end-1));
froot = fullfile('..','runs',rname);

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

figure('position',[245 774 990 920])
cm = cmocean('balance',101);

v0 = 0.1; % forcing amplitude
uclim = 0.1*v0;
vclim = 0.1*v0;

dat = load_data(froot,fids,datt.iter(end),{'UVEL','VVEL'});

yidx = find(gridm.Depth(1,:)==0,1,'last')+1;
xidx = find(gridm.Xp1>xSin1,1,'first'):find(abs(dat.UVEL(:,yidx,1))>1e-5,1,'last');

clf
subplot(211)
pcolor(gridm.Xp1,gridm.Y,squeeze(dat.UVEL(:,:,1))'); hold on
title('u (m/s)')
caxis([-1 1]*uclim);
%
colorbar
shading flat
colormap(gca,cm)
plot(gridm.Xp1(xidx([1 end])),gridm.Y(yidx)*[1 1],'r-','linewidth',2)
xlabel('x (km)')
ylabel('y (km)')
ylim([0 1e6]);
xticklabels(xticks/1e3);
yticklabels(yticks/1e3);
contour(gridm.XC',gridm.YC',gridm.Depth',[0:500:4e3],'k');
contour(gridm.XC',gridm.YC',gridm.Depth',[0 0],'color','k','linewidth',2);

subplot(212); cla
x = gridm.Xp1(xidx);
u = dat.UVEL(xidx,yidx,1);
x2 = x(1):10:x(end);
u2 = interp1(x,u,x2);

plot(x,u,'k-'); hold on
plot(xlim,[0 0],'k--')

zidx = find([false, diff(sign(u2))~=0]);
plot(x2(zidx),u2(zidx),'ko','markerfacecolor','k')

B = [1+0*zidx',[0:length(zidx)-1]']\(x2(zidx)');
plot(B(1)+B(2)*[0:length(zidx)-1]',0*zidx,'r.');

xlabel('x (km)')
xticklabels(xticks/1e3);
ylabel('u (m/s)')
grid on
title(sprintf('Best-fit Wavelength: %.2fkm',B(2)/1e3))

ttxt = sprintf('\\theta=%.1f^\\circ, kTopo=%.2f | T=%.1f cycles | [nx ny nz]=[%d %d %d]',...
               theta,kTopo,datt.T(end)/(2*pi/om),...
               length(gridm.X),length(gridm.Y),length(gridm.Z));
hax = axes('visible','off','position',[0 0 1 1]);
xlim(hax,[-1 1]); ylim(hax,[-1 1]);
text(hax,0,1,ttxt,...
     'fontweight','bold',...
     'verticalalignment','top',...
     'horizontalalignment','center');

dir_out = 'ctw_wavelength';
if ~exist(dir_out,'dir'); mkdir(dir_out); end
f_out = sprintf('ctw_wavelength_theta%3.2f_kTopo%.8f.jpg',theta,kTopo);
print('-djpeg90','-r300',fullfile(dir_out,f_out));
disp(['Saved ' f_out])

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
