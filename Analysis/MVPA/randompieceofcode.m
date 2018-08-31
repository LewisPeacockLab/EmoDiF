%rsa sanity check

    corr_matrix_previewcompare_full = zeros(60,60);
    for x = 1:60 %trial number
        
        for y = 1:60 %trial number
            
            corr_matrix_previewcompare_full(x, y) = corr(rsa_102.preview.mean.patterns(:,x), rsa_104.preview.mean.patterns(:,y));
       
        end
    end
    
    compare_fig = figure;
    set(compare_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(corr_matrix_previewcompare_full); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview 101 raw patterns averaged over TR %d: TR %d run 1',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('Preview 102 raw patterns averaged over TR %d: TR %d run 1',args.DFencode.meanTR_start, args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
    
    saveas(compare_fig, sprintf('emodif%s_compare1v2',subjNum),'png')