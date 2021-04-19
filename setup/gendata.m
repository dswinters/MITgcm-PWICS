function params = gendata(theta, kTopo, rdir, flags)

lprof = 20e3; % depth of corrugations
shelf_offset = 2.1*abs(lprof);

%% setup for wave generation along shelf topography at small promontory.
%% x is along shore
%% y is off shore
%% RCM Sept 2018

rname = sprintf('run_theta%3.2f_kTopo%.8f',theta,kTopo);
prec='real*4';
ieee='b';
fs = 6; fn = 'times';

% % use the setup from Jim Lerczak's model to set bathymetry & stratification
% jl = load('linearModel/subcritTopo_modes.mat'); % warning: in sigma coordinates!

g = 9.81; 
f = 1e-4; % gsw_f(43.29);
om = 1.36*f;
np = [4 2]; % processors [nx, ny]

% % this was a crude estimate neglecting the shelf
% cp = sqrt(g*max(abs(jl.h)));
% k = om/cp;
% lam = 2*pi/k;

% Ld = cp/f; % offshore length scale

modeno = 1; % which vertical wave mode to force with

% sponge cells on north, east and west boundaries
nsponge = 15;

% vertical grid parameters
Lz = 4000; % [m]

% vertical grid, exponential high res near surface
nzc = 30;
mindz = 25; maxdz = 500;
dz = smooth([ones(1,floor(700/mindz))*mindz logspace(log10(mindz),log10(maxdz),nzc)],5)';
dz = smooth(dz,5)';
tmp = [cumsum(dz)];
ind = find(tmp<Lz);
dze = Lz-max(abs(tmp(ind)));  % make sure that depth is Lz (I am not totally sure this is ok)
dz = [dz(ind) dze];
zf = -[0 cumsum(dz)]; % this is RF
z = 0.5*(zf(1:end-1)+zf(2:end));
nzc = length(dz);

%% now stratification - from Jim's linear code
rho0 = 999.8; g = 9.81; alpha = 2e-4;

r1 = 992; r2 = 995;
r0 = (r1+r2)/2; dr = r2-r1;
N2back = (2*pi/(0.5*3600))^2;
mupyc = 400;
Zpyc = -400;
r = r2 - 0.5*dr*(1+tanh((z-Zpyc)/mupyc)) - z*N2back*r0/g;
n2 = 0.5*(dr/r0)*(g/mupyc)*sech((z-Zpyc)/mupyc).^2 + N2back;

t = (1-r/rho0)/alpha+5;

% get vertical modes
n2tmp = fliplr(n2);
zrg = linspace(z(end),z(1),100); % need a regular grid for eigenmode code
n2rg = interp1(fliplr(z),n2tmp,zrg);
[wmodes,kpw,wt] = vmodes_w(zrg,n2rg,om,f);
wV = fliplr(interp1(zrg,wmodes(:,modeno),fliplr(z)));
wV = wV/max(abs(wV));
uV = fliplr(interp1(zrg,ddz(zrg,1)*wmodes(:,modeno),fliplr(z)));
uV = uV/max(abs(uV));
kV = kpw(modeno);
lamV = 2*pi/kV; % wavelength of mode

[wmom,kpp,wt] = vmodes_w(zrg,n2rg,om*1.01,f);
[wmop,kpm,wt] = vmodes_w(zrg,n2rg,om*0.99,f);
cg = (om*1.01 - om*0.99)/(kpp(1)-kpm(1));  % estimate group velocity

% cg = (Lz/pi)*(sqrt(om^2 -f^2)*(mean(n2)-om^2)^(3/2))/(om*(mean(n2)-f^2)); % an inferior estimate
% lamV = 188e3;
% mode 1 has lam ~ 335km
% mode 2 has lam ~ 188km

lTopo = 2*pi/kTopo;
LF = 1200e3; % length of forcing line
lamX = lamV/sind(theta); % projected wavelength of incident wave at coast
yF0 = 2*lamV+shelf_offset; % north position of forcing line

%horizontal grid parameters
dx_outer = lamV/10; % set the outer grid resolution to resolve the forcing wave
dy_outer = dx_outer;
dx_inner = 5e3; %2.e3;
dy_inner = 5e3; %2.e3;

% Expand the domain to accommodate large incident angles + some minimum sponge
% region.
high_res_pad = 300e3;
Lx = max(12*lamV, 2*(dx_outer*nsponge + (yF0+shelf_offset)*abs(tand(theta)) + LF/2 + high_res_pad));
Ly = 3*lamV; % in m

Tend = 10*(2*pi/om); % run for 10 forcing cycles
% disp(['Theta = ' num2str(theta) ': Simulation should end at ' num2str(Tend/3600,4) ' hr.'])

x0inner = Lx/2 - LF/2 - high_res_pad; % let's start the high res region just before the sine topo
x1inner = Lx/2 + LF/2 + high_res_pad; % how long after sine topo should we have high res?

xSin1 = x1inner - 500e3; % end corrugations at some constant offset from high-res region

xF0 = xSin1 - LF - (yF0-shelf_offset)*tand(theta); % west position of forcing line
xF1 = xF0 + LF;

% now make horizontal grids
% x
dx = [dx_outer*ones(1,ceil(x0inner/dx_outer)) ...
      dx_inner*ones(1,ceil((x1inner-x0inner)/dx_inner))...
      dx_outer*ones(1,ceil((Lx-x1inner)/dx_outer))];
% add cells on right to make domain divisible by # of processors
n_add = np(1) - mod(length(dx),np(1));
dx = [dx, dx_outer*ones(1,n_add)];

% smooth dx
dx = smooth(smooth(dx,5),5)';
xg = [0 cumsum(dx)];
xc = 0.5*(xg(2:end)+xg(1:end-1));

% y
dy = [dy_inner*ones(1,floor((yF0 + 0.25*lamV)/dy_inner)) ...
      dy_outer*ones(1,floor((Ly-0.75*lamV)/dy_outer))];
% add cells offshore to make domain divisible by # of processors
n_add = np(2) - mod(length(dy),np(2));
dy = [dy, dy_outer*ones(1,n_add)];

% smooth dy
dy = smooth(smooth(dy,5),5)';
yg = [0 cumsum(dy)];
yc = 0.5*(yg(2:end)+yg(1:end-1));

nxc = length(xc); nyc = length(yc);
Lx = max(xc); Ly = max(yc);

% rbcs mask: forcing will be applied within this region
[~,xi0] = min(abs(xc-xF0)); [~,yi0] = min(abs(yc-yF0));
[~,xi1] = min(abs(xc-xF1));
xind = xi0:xi1;
maskxy = zeros(nxc,nyc);
maskxy(xi0:xi1,yi0) = 1;
mask = repmat(maskxy,[1 1 nzc]);

% make topo, again from  Jim's linear code
hsh = 250; ysh = 75e3; dysl = 25e3;
prof = -hsh -0.5*(Lz-hsh)*(1+tanh((yc-ysh)/dysl));
%prof(1) = 0; % vertical wall
PROF = NaN(nxc,nyc);

% xSin0 = Lx/2 - LF/2; % start sine wave around here
% xSin1 = Lx/2 + LF/2; % stop sine wave after this

xSin1 = x1inner - 500e3; % stop sine wave after this
nLambdas = floor(LF/lTopo); % this must be odd, otherwise there will be a sharp cutoff at the ends
xSin0 = xSin1 - lTopo*nLambdas; % start sine wave around here

cutoff = 0*xc;
cutoff(xc >= xSin0 & xc <= xSin1) = 1;
for ii = 1:nxc
    ycshift = yc + shelf_offset + lprof*(cutoff(ii)* (0.5 - 0.5*cos((xc(ii)-xSin1)*kTopo)));
    tmp = interp1(ycshift,prof,yc);
    PROF(ii,:) = tmp;
end

bads = find(isnan(PROF)); % these are where min(ycshift) was > 0
PROF(bads) = 0;

%% output some parameters
params.t_end = Tend;
params.nxc = nxc;
params.nyc = nyc;
params.nzc = nzc;
params.np = np;

%% strat done, grids done, topo done

% initial fields
T = permute(repmat(t',[1 nxc nyc]),[2 3 1]);

U = 0*T;
V = 0*T;

% boundary sponge region fields

Uzonal = zeros(nyc,nzc,2); % [ny nz nt]
Vzonal = zeros(nyc,nzc,2);
Tzonal = repmat(squeeze(T(1,:,:)),[1 1 2]);

Umerid = zeros(nxc,nzc,2); % [nz nx nt]
Vmerid = zeros(nxc,nzc,2);
Tmerid = repmat(squeeze(T(:,end,:)),[1 1 2]);

% rbcs forcing: modal structures to be read into rbcs_fields_load
UMODE = permute(repmat(uV',[1 nxc nyc]),[2 3 1]);
UMODE = cat(4,UMODE,UMODE);

%% write some files
if flags.write_kTopo_dependent && flags.write_theta_dependent
    % These only need to be written once
    openfile =@(name) fopen(sprintf('../input/shared/%s.bin',name),'w',ieee);

    % dz and dy are independent of theta & kTopo
    fid=openfile('delY');
    fwrite(fid,dy,prec);
    fclose(fid);

    fid=openfile('delZ');
    fwrite(fid,dz,prec);
    fclose(fid);

    % zonal files only depend on Y
    fid=openfile('Uzonal');
    fwrite(fid,Uzonal,prec);
    fclose(fid);

    fid=openfile('Vzonal');
    fwrite(fid,Vzonal,prec);
    fclose(fid);

    fid=openfile('Tzonal');
    fwrite(fid,Tzonal,prec);
    fclose(fid);
end

if flags.write_kTopo_dependent
    % These only need to be written once per kTopo
    openfile =@(name) fopen(sprintf('../input/generated/kTopo%.8f_%s.bin',kTopo,name),'w',ieee);
end

if flags.write_theta_dependent
    % Written once per theta
    openfile =@(name) fopen(sprintf('../input/generated/theta%3.2f_%s.bin',theta,name),'w',ieee);

    % Anything depending on X depends on theta (I think)
    fid=openfile('Uinit');
    fwrite(fid,U,prec);
    fclose(fid);

    fid=openfile('Vinit');
    fwrite(fid,V,prec);
    fclose(fid);

    fid=openfile('Tinit');
    fwrite(fid,T,prec);
    fclose(fid);

    fid=openfile('delX');
    fwrite(fid,dx,prec);
    fclose(fid);

    fid=openfile('Umerid');
    fwrite(fid,Umerid,prec);
    fclose(fid);

    fid=openfile('Vmerid');
    fwrite(fid,Vmerid,prec);
    fclose(fid);

    fid=openfile('Tmerid');
    fwrite(fid,Tmerid,prec);
    fclose(fid);

    % RBCS forcing files depend on theta
    fid=openfile('umode');
    fwrite(fid,UMODE,prec);
    fclose(fid);

    fid=openfile('mask');
    fwrite(fid,mask,prec);
    fclose(fid);
end

% Everything else is specific to each run
openfile =@(name) fopen(fullfile(rdir,[name '.bin']),'w',ieee);
fid=openfile('topog');
fwrite(fid,PROF,prec);
fclose(fid);

% Output corrugation parameters to .mat file in run directory so we don't need
% to copy/paste all of gendata when calculating fluxes etc.
save(fullfile(rdir, 'corrugation_params.mat'),'lprof','xSin0','xSin1');
disp(['Saved ' fullfile(rdir, 'corrugation_params.mat')])

% Make setup figure
figure('visible','on')
set(gcf,'position',[0 0 1000 300],'color','w')
subplot('position',[0.08 0.1 0.15 0.8])
plot(yc/1e3,prof)
title('Bathymetry')
set(gca,'fontsize',fs,'fontname',fn)

subplot('position',[0.3 0.1 0.15 0.8])
plot(t,z)
grid on; hold on
plot(t,z,'.')
title('Temperature')
set(gca,'fontsize',fs,'fontname',fn)

subplot('position',[0.48 0.3 0.05 0.6])
plot(dy/1e3,yc/1e3)
grid on; hold on
plot(dy/1e3,yc/1e3,'.')
ylim([0 Ly/1e3])
set(gca,'fontsize',fs,'fontname',fn,'layer','top')
ylabel('y [km]')

subplot('position',[0.55 0.3 0.15 0.6])
set(gcf,'paperpositionmode','auto','color','w')
pcolor(xc,yc,mask(:,:,1)'); shading flat;
cb = colorbar;
set(cb,'position',[0.71 0.3 0.015 0.6])
hold on
contour(xc,yc,PROF',[-4000:500:0],'k')
plot([xc(nsponge) xc(nsponge)],[yc(1) yc(end)],'k')
plot([xc(nxc-nsponge) xc(nxc-nsponge)],[yc(1) yc(end)],'k')
plot([xc(1) xc(end)],[yc(nyc-nsponge) yc(nyc-nsponge)],'k')
contour(xc,yc,PROF',[0 0],'linewidth',2,'color','k');
yf = shelf_offset + [0 100e3]; % y-coords of flux calc region
plot((xSin1*[1 1]+70e3),yf,'r-','linewidth',1);
patch([xF0 xF1, xSin1 xSin1-LF],...
      [yF0*[1 1], shelf_offset*[1 1]],...
      [1 1 1],'edgecolor','none','facealpha',0.5)
caxis([-1 1])
grid on
xlim([0 Lx])
ylim([0 Ly])
xticklabels(xticks/1e3)
yticklabels([])
title('MASK')
set(gca,'fontsize',fs,'fontname',fn,'layer','top','yticklabel',{})

subplot('position',[0.55 0.1 0.15 0.1])
plot(xc/1e3,dx/1e3)
hold on
plot(xc/1e3,dx/1e3,'.')
plot(xSin0/1e3*[1 1],ylim,'k--')
plot(xSin1/1e3*[1 1],ylim,'k--')
xlim([0 Lx/1e3])
grid on
set(gca,'fontsize',fs,'fontname',fn)
xlabel('x [km]')
ylabel('dx [km]')

subplot('position',[0.78 0.1 0.17 0.8])
set(gcf,'paperpositionmode','auto','color','w')
plot(uV,z)
grid on
% pcolor(yc/1e3,z,squeeze(UMODE(1,:,:,1))'); shading flat; colorbar
% hold on
% plot(yc/1e3,prof,'k')
% caxis([-1 1])
title({['modal structure, \lambda = ' num2str(lamV/1e3,4) ' km,'],['\omega = ' num2str(om,3) 'rad/s']})
set(gca,'fontsize',fs,'fontname',fn)
xlabel('U ')

fout = fullfile('figures',sprintf('setup_%s.jpg',rname));
print('-djpeg90','-r300',fout)
disp(['Saved ' fout])
close all
