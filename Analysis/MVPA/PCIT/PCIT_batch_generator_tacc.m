for x = 1:200
fname = fopen('run_commands.txt','a')
fprintf(fname, 'batch_matlab.sh "imdif_bootstrap_run_importance_sampler_tacc(1,%d)"\n',x)
fclose(fname)
end