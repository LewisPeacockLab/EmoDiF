function [] = imdif_create_pcit_jobs(submit)
%does not require any arguments


SUBMIT_JOBS = submit;

REUSE_ANALYSIS_ID = false;
analysis_id_touse = '130417_1225_17_030';

HICONF_HITS = false;
if HICONF_HITS; pcit_dist = 'bernoulli';
else pcit_dist = 'normal';
end
  
% create job .m file

fprintf('+ Creating job ==> %s\n',filename);

fid = fopen(filename, 'w');
fprintf(fid, 'function[] = importance_sampler_%s_r%d()\n', analysis_id, run_number);
fprintf(fid, '[analysis_settings raw_data] = repref_configure_importance_sampler(''%s'', ''%s'');\n',analysis_id, data_id);

fprintf(fid, 'analysis_settings.distribution = ''%s'';\n', pcit_dist);

if DO_PCIT_SCRAMBLE
    fprintf(fid, 'analysis_settings.scramble = true;\n');
    fprintf(fid, 'analysis_settings.scramble_run = %d;\n', SCRAMBLED_ITERS);
    fprintf(fid, 'analysis_settings.scramble_style = ''across_subjects_across_categories'';\n');
end

fprintf(fid, 'imdif_importance_sampler(raw_data, analysis_settings, %d)\n', run_number);
fclose(fid);


% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

% submit job to rondo cluster @ Princeton
if SUBMIT_JOBS
    
    cmd = sprintf('submit_long -l vf=14G %s',filename);
    [s r] = system(cmd);
    fprintf('+ Submitting job ==> %s',r);
    if ~REUSE_ANALYSIS_ID; fprintf('  Data: %s\n',data_id); end
    fprintf('\n');
    pause(0.25); % pause for 0.25 sec
    
end

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


end