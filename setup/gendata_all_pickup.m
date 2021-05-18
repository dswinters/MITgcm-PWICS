clear all, close all

deltaT = 500;
theta = [-60:20:60]; % angle of PW wave vector wrt line perpendicular to coast [deg]
lTopo = linspace(20e3,700e3,30);
kTopo = 2*pi./lTopo;
params = gendata_params();

% Some flags so we can generate files on e.g. the first theta of every kTopo,
% only once per group of runs, etc...
flags.write_theta_dependent = true;
flags.write_kTopo_dependent = true;
output_freq = -3600; % write every 3600 seconds


for i = 1:length(theta) % theta index
    thetaPrefix = sprintf('theta%3.2f_',theta(i)); % File prefix for theta
    flags.write_theta_dependent = true;

    for j = 1:length(kTopo) % kTopo index
        kTopoPrefix = sprintf('kTopo%.8f_',kTopo(j)); % File prefix for kTopo

        % Set up run directories
        rname = sprintf('run_%s%s',thetaPrefix,kTopoPrefix(1:end-1));
        rdir = fullfile('..','runs',rname);
        if ~exist(rdir,'dir'); mkdir(rdir); end

        % Generate run data
        params = gendata(theta(i), kTopo(j), rdir, flags);
        niter0 = floor(params.t_chkpt/deltaT);

        % Template files have been set up with FIELD = PLACEHOLDER for fields that need
        % to be written. Make a cell array of substitutions with entries that
        % look like: {filename, {field1, string; ... fieldN, string}, prefix}
        other_subs = {};
        thsubs = {};
        ktoposubs = {};

        % Substitutions depending only on theta
        if flags.write_theta_dependent % (i.e. do these on the first kTopo of every theta)
            thsubs = {
                '../code/templates/rbcs_fields_load.F',...
                    {'th', sprintf('%.1f * pi/180.',theta(i));
                     'om', sprintf('%.3e * pi/180.',params.om)},...
                     thetaPrefix;
                '../code/templates/SIZE.h',...
                    {'sNx', sprintf('%d', params.nxc/params.np(1));
                     'sNy', sprintf('%d', params.nyc/params.np(2));
                     'nPx', sprintf('%d', params.np(1));
                     'nPy', sprintf('%d', params.np(2));
                     'Nr', sprintf('%d', params.nzc)},...
                     thetaPrefix;
                '../input/templates/data',...
                    {'deltaT', sprintf('%.1f',deltaT);
                     'pChkptFreq', sprintf('%.1f',params.t_chkpt);
                     'niter0', sprintf('%d',niter0);
                     'nTimeSteps',  sprintf('%d',ceil(params.t_end/deltaT)-niter0)},...
                     thetaPrefix;
                '../input/templates/data.diagnostics',...
                     {'frequency(3)', sprintf('%.1f',output_freq);
                      'timePhase(3)', sprintf('%.1f',output_freq);
                      'frequency(4)', sprintf('%.1f',output_freq);
                      'timePhase(4)', sprintf('%.1f',output_freq)},...
                     thetaPrefix;
                '../input/templates/data.obcs',...
                    {'OB_Ieast', sprintf('%d*-1', params.nyc);
                     'OB_Iwest', sprintf('%d*1', params.nyc);
                     'OB_Jnorth', sprintf('%d*-1', params.nxc);
                     'OB_Jsouth', sprintf('%d*1', params.nxc)},...
                     thetaPrefix;
                     };
        end

        if flags.write_kTopo_dependent % Substitutions depending only on kTopo
            % I don't think we have any of these
        end

        substitutions = cat(1,thsubs,ktoposubs);
        for nf = 1:size(substitutions,1)
            ftxt = fileread(substitutions{nf,1}); % load template file text
            for ns = 1:size(substitutions{nf,2},1)
                expr = sprintf('%s *= *PLACEHOLDER',substitutions{nf,2}{ns,1}); % create regexp
                expr = strrep(expr,'(','\('); % Replace parens for regexp function
                expr = strrep(expr,')','\)');

                [i1,i2] = regexp(ftxt,expr,'start','end'); % find placeholder
                % update placeholder text with substitution text
                ftxt = cat(2, ...
                           ftxt(1:i1-1), ...
                           strrep(ftxt(i1:i2),'PLACEHOLDER',substitutions{nf,2}{ns,2}),...
                           ftxt(i2+1:end));
            end
            % Write new text to file named with prefix
            [fdir,fname,fext] = fileparts(substitutions{nf,1});
            dir_out = strrep(fdir,'templates','generated');
            newfile = fullfile(dir_out,sprintf('%s%s%s',substitutions{nf,3},fname,fext));
            fid = fopen(newfile,'w');
            fprintf(fid,'%s\n',ftxt);
            fclose(fid);
        end % doing substitutions in template files

        % run shell script to link files
        if flags.write_theta_dependent; build_flag='--build'; else build_flag=''; end
        cmd=sprintf('./pwics_setup.sh %s %s %s',thetaPrefix,kTopoPrefix,build_flag);
        system(cmd);

        flags.write_theta_dependent = false; % We've finished (at least) one kTopo
    end % loop over kTopos

    flags.write_kTopo_dependent = false; % We've looped over all kTopos
end % loop over thetas

!ls ../runs > run_list.txt
