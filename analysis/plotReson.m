function [aflxBC aflx] = plotReson(theta,kTopo)

fs = 14; fn = 'times';
lam = 2*pi./kTopo;
dx_inner = 5e3; %2.e3;
dy_inner = 5e3; %2.e3;
om = 1.36*1e-4;

datt = [];
aflxBC = {[]};
aflx = {[]};
yf = {};

for i = 1:length(theta)
    for j = 1:length(kTopo)
        % Load run data
        thetaPrefix = sprintf('theta%3.2f_',theta(i)); % File prefix for theta
        kTopoPrefix = sprintf('kTopo%.8f_',kTopo(j)); % File prefix for kTopo
        rdir = '..';
        rname = sprintf('run_%s%s',thetaPrefix,kTopoPrefix(1:end-1));
        load(fullfile(rdir,'runs',rname,'corrugation_params.mat'),'xSin1')
        datt = rdmnc(fullfile(rdir,'runs',rname,'outs_sn.*'),'T','iter');
        fluxFile = fullfile(rdir,'analysis','calcFluxBC',[rname '_flux_wave_avg_end.mat']);

        if exist(fluxFile,'file')
            fl = load(fluxFile);
            if j==1
                gridm = rdmnc(fullfile(rdir,'runs',rname,'grid*'));
                yf{i} = fl.yc(find(gridm.Depth(1,:)==0,1,'last')) + [0 100e3];
            end
            indy = fl.yc >= yf{i}(1) & fl.yc <= yf{i}(2);


            % X location where we computed fluxes (70km after end of sine topog)
            [~,indx] = min(abs(fl.xc-xSin1+70e3));

            aflxBC{i}(j) = squeeze(nansum(fl.upbcNz(indx,indy,:)*dy_inner,2));
            aflx{i}(j) = squeeze(nansum(fl.upz(indx,indy,:)*dy_inner,2));
        else
            aflxBC{i}(j) = nan;
            aflx{i}(j) = nan;
        end
        fprintf('\r%s [%d of %d]',fluxFile,j,length(kTopo))
    end
end
fprintf('\n')
save plotReson.mat aflxBC aflx indx datt yf

figure('position',[0 0 1065 400],'color','w','paperpositionmode','auto');

% subplot(2,2,3)
pl1 = [];
pl2 = [];
cols = get(groot,'DefaultAxesColorOrder');
clf

for i = 1:length(aflx)
    subplot(1,2,1); hold on
    pl1(i) = plot(lam,aflxBC{i}(end,:),'-','linewidth',1,'color',cols(i,:)); grid on
    xlabel('\lambda Topo')
    ylabel('Baroclinic Flux (Along-shelf/incident)')
    set(gca,'fontsize',8)

    % subplot(2,2,4)
    subplot(1,2,2); hold on
    pl2(i) = plot(lam,aflx{i}(end,:),'-','linewidth',1,'color',cols(i,:)); grid on
    xlabel('\lambda Topo')
    ylabel('Total Flux (Along-shelf/Incident)')
    set(gca,'fontsize',8)
end

legstr = arrayfun(@(x) sprintf('\\theta=%.1f^\\circ',x), theta, 'uni', false);

hl = legend(pl1,legstr,'location','northeast');
hl = legend(pl2,legstr,'location','northeast');

print('-dpng','pwics_resonance.png')
save plotReson.mat aflxBC aflx theta kTopo
