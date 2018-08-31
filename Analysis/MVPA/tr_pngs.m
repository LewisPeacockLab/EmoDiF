%load rsa structure


corr_matrix_match=rsa.DFencode.bytrialTR.smatrixbytr.corr_matrix_match

tr1 = corr_matrix_match(1).corr_matrix_match;
tr1_fig = figure;
    set(tr1_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(tr1); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d tr 1',3, 5),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d tr 1',1, 3),'FontSize',15,'FontWeight','bold');
    
    saveas(tr1_fig, sprintf('emodif%s_tr1','103'),'png')
    
    tr2 = corr_matrix_match(2).corr_matrix_match;
    tr2_fig = figure;
    set(tr2_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(tr2); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d tr 1',3, 5),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d tr 1',1, 3),'FontSize',15,'FontWeight','bold');
    
    saveas(tr1_fig, sprintf('emodif%s_tr2','103'),'png')
    
        tr3 = corr_matrix_match(3).corr_matrix_match;
        tr3_fig = figure;
    set(tr3_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(tr3); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d tr 1',3, 5),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d tr 1',1, 3),'FontSize',15,'FontWeight','bold');
    
    saveas(tr3_fig, sprintf('emodif%s_tr3','103'),'png')
    
            tr4 = corr_matrix_match(4).corr_matrix_match;
        tr4_fig = figure;
    set(tr4_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(tr4); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d tr 1',3, 5),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d tr 1',1, 3),'FontSize',15,'FontWeight','bold');
    
    saveas(tr4_fig, sprintf('emodif%s_tr4','103'),'png')
    
                
    
    
    tr5 = corr_matrix_match(5).corr_matrix_match;
                
                
        tr5_fig = figure;
    set(tr5_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(tr5); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d tr 1',3, 5),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d tr 1',1, 3),'FontSize',15,'FontWeight','bold');
    
    saveas(tr5_fig, sprintf('emodif%s_tr5','103'),'png')
    
        tr6 = corr_matrix_match(6).corr_matrix_match;
                
                
        tr6_fig = figure;
    set(tr6_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(tr6); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d tr 1',3, 6),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d tr 1',1, 5),'FontSize',15,'FontWeight','bold');
    
    saveas(tr6_fig, sprintf('emodif%s_tr6','103'),'png')
    
            tr7 = corr_matrix_match(7).corr_matrix_match;
                
                
        tr7_fig = figure;
    set(tr7_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(tr7); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d tr 1',3, 5),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d tr 1',1, 5),'FontSize',15,'FontWeight','bold');
    
    saveas(tr7_fig, sprintf('emodif%s_tr7','103'),'png')