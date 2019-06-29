function emodif_mvpa_decoding_bar_compare_fsor(subjNum,test_phase, maskName, test_date)
%emodif_mvpa_decoding_bar_compare_fsor('101','DFencode','tempoccfusi_pHg_LOC_combined_epi_space','25-Jun-2018')
%emodif_mvpa_decoding_bar_compare('102','DFencode','03-May-2018')
%emodif_mvpa_decoding_bar_compare('103','DFencode','07-May-2018')
%emodif_mvpa_decoding_bar_compare('104','DFencode','13-Jun-2018')
%emodif creation of bar graphs from result file - DFencode

 version = '2018Jul12';
 
  args.subjNum = subjNum;
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  %for astoria
   args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
% for tigger
%    args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);

  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.output_dir = sprintf('%s/results/%s/%s/%s',args.subj_dir, test_phase, maskName, test_date);
  args.script_dir = pwd;
  
  cd(args.output_dir)
if strcmp(test_phase, 'DFencode') == true
    load('results_DFencode.mat');    %this loads the result file 
    
    %figure out how many categories this parsed file has
    
    categories = results.remember.all.acts;
    
    
    %Remember - PRECUE
    
    results_rem_precue_acts_face = nanmean(results.remember.all.acts.face(:,1:3),2); %merging 1-3TR
    results.remember.all.bar.precue.face.bytrial = results_rem_precue_acts_face;
    results.remember.all.bar.precue.face.nanmean = nanmean(results.remember.all.bar.precue.face.bytrial);
    results.remember.all.bar.precue.face.std = std(results.remember.all.bar.precue.face.bytrial);
    
    results_rem_precue_acts_scene = nanmean(results.remember.all.acts.scene(:,1:3),2); %merging 1-3TR
    results.remember.all.bar.precue.scene.bytrial = results_rem_precue_acts_scene;
    results.remember.all.bar.precue.scene.nanmean = nanmean(results.remember.all.bar.precue.scene.bytrial);
    results.remember.all.bar.precue.scene.std = std(results.remember.all.bar.precue.scene.bytrial);
    
    results_rem_precue_acts_object = nanmean(results.remember.all.acts.object(:,1:3),2); %merging 1-3TR
    results.remember.all.bar.precue.object.bytrial = results_rem_precue_acts_object;
    results.remember.all.bar.precue.object.nanmean = nanmean(results.remember.all.bar.precue.object.bytrial);
    results.remember.all.bar.precue.object.std = std(results.remember.all.bar.precue.object.bytrial);
    
    results_rem_precue_acts_rest = nanmean(results.remember.all.acts.rest(:,1:3),2); %merging 1-3TR
    results.remember.all.bar.precue.rest.bytrial = results_rem_precue_acts_rest;
    results.remember.all.bar.precue.rest.nanmean = nanmean(results.remember.all.bar.precue.rest.bytrial);
    results.remember.all.bar.precue.rest.std = std(results.remember.all.bar.precue.rest.bytrial);
    
    results.remember.all.bar.precue.scene_facenorm.bytrial = results.remember.all.bar.precue.scene.bytrial-results.remember.all.bar.precue.face.bytrial;
    results.remember.all.bar.precue.scene_facenorm.nanmean = nanmean(results.remember.all.bar.precue.scene_facenorm.bytrial);
    results.remember.all.bar.precue.scene_facenorm.std = std(results.remember.all.bar.precue.scene_facenorm.bytrial);
    
    
    %Remember - POSTCUE
    
    results_rem_postcue_acts_face = nanmean(results.remember.all.acts.face(:,5:7),2); %merging 1-3TR
    results.remember.all.bar.postcue.face.bytrial = results_rem_postcue_acts_face;
    results.remember.all.bar.postcue.face.nanmean = nanmean(results.remember.all.bar.postcue.face.bytrial );
    results.remember.all.bar.postcue.face.std = std(results.remember.all.bar.postcue.face.bytrial );
    
    results_rem_postcue_acts_scene = nanmean(results.remember.all.acts.scene(:,5:7),2); %merging 1-3TR
    results.remember.all.bar.postcue.scene.bytrial = results_rem_postcue_acts_scene;
    results.remember.all.bar.postcue.scene.nanmean = nanmean(results.remember.all.bar.postcue.scene.bytrial );
    results.remember.all.bar.postcue.scene.std = std(results.remember.all.bar.postcue.scene.bytrial );
    
    results_rem_postcue_acts_object = nanmean(results.remember.all.acts.object(:,5:7),2); %merging 1-3TR
    results.remember.all.bar.postcue.object.bytrial = results_rem_postcue_acts_object;
    results.remember.all.bar.postcue.object.nanmean = nanmean(results.remember.all.bar.postcue.object.bytrial );
    results.remember.all.bar.postcue.object.std = std(results.remember.all.bar.postcue.object.bytrial );
    
    results_rem_postcue_acts_rest = nanmean(results.remember.all.acts.rest(:,5:7),2); %merging 1-3TR
    results.remember.all.bar.postcue.rest.bytrial = results_rem_postcue_acts_rest;
    results.remember.all.bar.postcue.rest.nanmean = nanmean(results.remember.all.bar.postcue.rest.bytrial );
    results.remember.all.bar.postcue.rest.std = std(results.remember.all.bar.postcue.rest.bytrial);
    
    results.remember.all.bar.postcue.scene_facenorm.bytrial = results.remember.all.bar.postcue.scene.bytrial-results.remember.all.bar.postcue.face.bytrial;
    results.remember.all.bar.postcue.scene_facenorm.nanmean = nanmean(results.remember.all.bar.postcue.scene_facenorm.bytrial);
    results.remember.all.bar.postcue.scene_facenorm.std = std(results.remember.all.bar.postcue.scene_facenorm.bytrial);
    
    results.remember.all.bar.postminuspre.facenorm.scene.bytrial= results.remember.all.bar.postcue.scene_facenorm.bytrial - results.remember.all.bar.precue.scene_facenorm.bytrial;
    results.remember.all.bar.postminuspre.facenorm.scene.nanmean = nanmean(results.remember.all.bar.postminuspre.facenorm.scene.bytrial);
    results.remember.all.bar.postminuspre.facenorm.scene.std = std(results.remember.all.bar.postminuspre.facenorm.scene.nanmean);
    
    results.remember.all.bar.postminuspre.nofacenorm.scene.bytrial = results.remember.all.bar.postcue.scene.bytrial -results.remember.all.bar.precue.scene.bytrial;
    results.remember.all.bar.postminuspre.nofacenorm.scene.nanmean = nanmean(results.remember.all.bar.postminuspre.nofacenorm.scene.bytrial);
    results.remember.all.bar.postminuspre.nofacenorm.scene.std = std(results.remember.all.bar.postminuspre.nofacenorm.scene.nanmean);
    
    %Forget  - PRECUE
    
    results_for_precue_acts_face = nanmean(results.forget.all.acts.face(:,1:3),2); %merging 1-3TR
    results.forget.all.bar.precue.face.bytrial = results_for_precue_acts_face;
    results.forget.all.bar.precue.face.nanmean = nanmean(results.forget.all.bar.precue.face.bytrial);
    results.forget.all.bar.precue.face.std = std(results.forget.all.bar.precue.face.bytrial);
    
    results_for_precue_acts_scene = nanmean(results.forget.all.acts.scene(:,1:3),2); %merging 1-3TR
    results.forget.all.bar.precue.scene.bytrial = results_for_precue_acts_scene;
    results.forget.all.bar.precue.scene.nanmean = nanmean(results.forget.all.bar.precue.scene.bytrial);
    results.forget.all.bar.precue.scene.std = std(results.forget.all.bar.precue.scene.bytrial);
    
    results_for_precue_acts_object = nanmean(results.forget.all.acts.object(:,1:3),2); %merging 1-3TR
    results.forget.all.bar.precue.object.bytrial = results_for_precue_acts_object;
    results.forget.all.bar.precue.object.nanmean = nanmean(results.forget.all.bar.precue.object.bytrial);
    results.forget.all.bar.precue.object.std = std(results.forget.all.bar.precue.object.bytrial);
   
    results_for_precue_acts_rest = nanmean(results.forget.all.acts.rest(:,1:3),2); %merging 1-3TR
    results.forget.all.bar.precue.rest.bytrial = results_for_precue_acts_rest;
    results.forget.all.bar.precue.rest.nanmean = nanmean(results.forget.all.bar.precue.rest.bytrial);
    results.forget.all.bar.precue.rest.std = std(results.forget.all.bar.precue.rest.bytrial);
    
    results.forget.all.bar.precue.scene_facenorm.bytrial = results.forget.all.bar.precue.scene.bytrial-results.forget.all.bar.precue.face.bytrial;
    results.forget.all.bar.precue.scene_facenorm.nanmean = nanmean(results.forget.all.bar.precue.scene_facenorm.bytrial);
    results.forget.all.bar.precue.scene_facenorm.std = std(results.forget.all.bar.precue.scene_facenorm.bytrial);
    
    %forget - POSTCUE
    
    results_for_postcue_acts_face = nanmean(results.forget.all.acts.face(:,5:7),2); %merging 1-3TR
    results.forget.all.bar.postcue.face.bytrial = results_for_postcue_acts_face;
    results.forget.all.bar.postcue.face.nanmean = nanmean(results.forget.all.bar.postcue.face.bytrial);
    results.forget.all.bar.postcue.face.std = std(results.forget.all.bar.postcue.face.bytrial);
    
    results_for_postcue_acts_scene = nanmean(results.forget.all.acts.scene(:,5:7),2); %merging 1-3TR
    results.forget.all.bar.postcue.scene.bytrial = results_for_postcue_acts_scene;
    results.forget.all.bar.postcue.scene.nanmean = nanmean(results.forget.all.bar.postcue.scene.bytrial);
    results.forget.all.bar.postcue.scene.std = std(results.forget.all.bar.postcue.scene.bytrial);
    
    results_for_postcue_acts_object = nanmean(results.forget.all.acts.object(:,5:7),2); %merging 1-3TR
    results.forget.all.bar.postcue.object.bytrial = results_for_postcue_acts_object;
    results.forget.all.bar.postcue.object.nanmean = nanmean(results.forget.all.bar.postcue.object.bytrial);
    results.forget.all.bar.postcue.object.std = std(results.forget.all.bar.postcue.object.bytrial);
      
    results_for_postcue_acts_rest = nanmean(results.forget.all.acts.rest(:,5:7),2); %merging 1-3TR
    results.forget.all.bar.postcue.rest.bytrial = results_for_postcue_acts_rest;
    results.forget.all.bar.postcue.rest.nanmean = nanmean(results.forget.all.bar.postcue.rest.bytrial);
    results.forget.all.bar.postcue.rest.std = std(results.forget.all.bar.postcue.rest.bytrial);

    results.forget.all.bar.postcue.scene_facenorm.bytrial = results.forget.all.bar.postcue.scene.bytrial-results.forget.all.bar.postcue.face.bytrial;
    
    results.forget.all.bar.postcue.scene_facenorm.bytrial = results.forget.all.bar.postcue.scene.bytrial-results.forget.all.bar.postcue.face.bytrial;
    results.forget.all.bar.postcue.scene_facenorm.nanmean = nanmean(results.forget.all.bar.postcue.scene_facenorm.bytrial);
    results.forget.all.bar.postcue.scene_facenorm.std = std(results.forget.all.bar.postcue.scene_facenorm.bytrial);
    
    results.forget.all.bar.postminuspre.facenorm.scene.bytrial= results.forget.all.bar.postcue.scene_facenorm.bytrial - results.forget.all.bar.precue.scene_facenorm.bytrial;
    results.forget.all.bar.postminuspre.facenorm.scene.nanmean = nanmean(results.forget.all.bar.postminuspre.facenorm.scene.bytrial);
    results.forget.all.bar.postminuspre.facenorm.scene.std = std(results.forget.all.bar.postminuspre.facenorm.scene.nanmean);
        
    results.forget.all.bar.postminuspre.nofacenorm.scene.bytrial = results.forget.all.bar.postcue.scene.bytrial -results.forget.all.bar.precue.scene.bytrial;
    results.forget.all.bar.postminuspre.nofacenorm.scene.nanmean = nanmean(results.forget.all.bar.postminuspre.nofacenorm.scene.bytrial);
    results.forget.all.bar.postminuspre.nofacenorm.scene.std = std(results.forget.all.bar.postminuspre.nofacenorm.scene.nanmean);
    
    %remember - MERGED - OLD - PRECUE
    results_rem_merged_old_precue_acts_face = nanmean(results.remember.merged.acts.OLD.face(:,1:3),2); %merging 1-3TR
    results.remember.merged.bar.OLD.precue.face.bytrial =  results_rem_merged_old_precue_acts_face;
    results.remember.merged.bar.OLD.precue.face.nanmean = nanmean(results.remember.merged.bar.OLD.precue.face.bytrial);
    results.remember.merged.bar.OLD.precue.face.std = std(results.remember.merged.bar.OLD.precue.face.bytrial);
    
    results_rem_merged_old_precue_acts_scene = nanmean(results.remember.merged.acts.OLD.scene(:,1:3),2); %merging 1-3TR
    results.remember.merged.bar.OLD.precue.scene.bytrial =  results_rem_merged_old_precue_acts_scene;
    results.remember.merged.bar.OLD.precue.scene.nanmean = nanmean(results.remember.merged.bar.OLD.precue.scene.bytrial);
    results.remember.merged.bar.OLD.precue.scene.std = std(results.remember.merged.bar.OLD.precue.scene.bytrial);
    
    results_rem_merged_old_precue_acts_object = nanmean(results.remember.merged.acts.OLD.object(:,1:3),2); %merging 1-3TR
    results.remember.merged.bar.OLD.precue.object.bytrial =  results_rem_merged_old_precue_acts_object;
    results.remember.merged.bar.OLD.precue.object.nanmean = nanmean(results.remember.merged.bar.OLD.precue.object.bytrial);
    results.remember.merged.bar.OLD.precue.object.std = std(results.remember.merged.bar.OLD.precue.object.bytrial);
    
    results_rem_merged_old_precue_acts_rest = nanmean(results.remember.merged.acts.OLD.rest(:,1:3),2); %merging 1-3TR
    results.remember.merged.bar.OLD.precue.rest.bytrial =  results_rem_merged_old_precue_acts_rest;
    results.remember.merged.bar.OLD.precue.rest.nanmean = nanmean(results.remember.merged.bar.OLD.precue.rest.bytrial);
    results.remember.merged.bar.OLD.precue.rest.std = std(results.remember.merged.bar.OLD.precue.rest.bytrial);
  % 
  
  %Remember - MERGED - OLD - POSTCUE
  
    results_rem_merged_old_postcue_acts_face = nanmean(results.remember.merged.acts.OLD.face(:,5:7),2); %merging 1-3TR
    results.remember.merged.bar.OLD.postcue.face.bytrial =  results_rem_merged_old_postcue_acts_face;
    results.remember.merged.bar.OLD.postcue.face.nanmean = nanmean(results.remember.merged.bar.OLD.postcue.face.bytrial);
    results.remember.merged.bar.OLD.postcue.face.std = std(results.remember.merged.bar.OLD.postcue.face.bytrial);
    
    results_rem_merged_old_postcue_acts_scene = nanmean(results.remember.merged.acts.OLD.scene(:,5:7),2); %merging 1-3TR
    results.remember.merged.bar.OLD.postcue.scene.bytrial =  results_rem_merged_old_postcue_acts_scene;
    results.remember.merged.bar.OLD.postcue.scene.nanmean = nanmean(results.remember.merged.bar.OLD.postcue.scene.bytrial);
    results.remember.merged.bar.OLD.postcue.scene.std = std(results.remember.merged.bar.OLD.postcue.scene.bytrial);
    
    results_rem_merged_old_postcue_acts_object = nanmean(results.remember.merged.acts.OLD.object(:,5:7),2); %merging 1-3TR
    results.remember.merged.bar.OLD.postcue.object.bytrial =  results_rem_merged_old_postcue_acts_object;
    results.remember.merged.bar.OLD.postcue.object.nanmean = nanmean(results.remember.merged.bar.OLD.postcue.object.bytrial);
    results.remember.merged.bar.OLD.postcue.object.std = std(results.remember.merged.bar.OLD.postcue.object.bytrial);
       
    results_rem_merged_old_postcue_acts_rest = nanmean(results.remember.merged.acts.OLD.rest(:,5:7),2); %merging 1-3TR
    results.remember.merged.bar.OLD.postcue.rest.bytrial =  results_rem_merged_old_postcue_acts_rest;
    results.remember.merged.bar.OLD.postcue.rest.nanmean = nanmean(results.remember.merged.bar.OLD.postcue.rest.bytrial);
    results.remember.merged.bar.OLD.postcue.rest.std = std(results.remember.merged.bar.OLD.postcue.rest.bytrial);
    
    %Remember - MERGED - NEW- PRECUE
    if subjNum == '103' | subjNum == '113' | subjNum == '115' | subjNum == '114' | subjNum == '119';% has 0 here 
        
        results.remember.merged.bar.NEW.precue.face.bytrial =  NaN;
        results.remember.merged.bar.NEW.precue.face.nanmean = NaN;
        results.remember.merged.bar.NEW.precue.face.std = NaN;
        
        
        results.remember.merged.bar.NEW.precue.scene.bytrial =  NaN;
        results.remember.merged.bar.NEW.precue.scene.nanmean = NaN;
        results.remember.merged.bar.NEW.precue.scene.std = NaN;
        
        
        results.remember.merged.bar.NEW.precue.object.bytrial =  NaN;
        results.remember.merged.bar.NEW.precue.object.nanmean = NaN;
        results.remember.merged.bar.NEW.precue.object.std = NaN;
              
        
        results.remember.merged.bar.NEW.precue.rest.bytrial =  NaN;
        results.remember.merged.bar.NEW.precue.rest.nanmean = NaN;
        results.remember.merged.bar.NEW.precue.rest.std = NaN;
        
        %
        
        %Remember - MERGED - NEW- POSTCUE
        
        
        results.remember.merged.bar.NEW.postcue.face.bytrial =  NaN;
        results.remember.merged.bar.NEW.postcue.face.nanmean = NaN;
        results.remember.merged.bar.NEW.postcue.face.std = NaN;
        
        
        results.remember.merged.bar.NEW.postcue.scene.bytrial =  NaN;
        results.remember.merged.bar.NEW.postcue.scene.nanmean = NaN;
        results.remember.merged.bar.NEW.postcue.scene.std = NaN;
        
        
        results.remember.merged.bar.NEW.postcue.object.bytrial =  NaN;
        results.remember.merged.bar.NEW.postcue.object.nanmean = NaN;
        results.remember.merged.bar.NEW.postcue.object.std = NaN;
            
        
        results.remember.merged.bar.NEW.postcue.rest.bytrial =  NaN;
        results.remember.merged.bar.NEW.postcue.rest.nanmean = NaN;
        results.remember.merged.bar.NEW.postcue.rest.std = NaN;
        
    else
        
        % Remember - Merged - New - Precue
        results_rem_merged_new_precue_acts_face = nanmean(results.remember.merged.acts.NEW.face(:,1:3),2); %merging 1-3TR
        results.remember.merged.bar.NEW.precue.face.bytrial =  results_rem_merged_new_precue_acts_face;
        results.remember.merged.bar.NEW.precue.face.nanmean = nanmean(results.remember.merged.bar.NEW.precue.face.bytrial);
        results.remember.merged.bar.NEW.precue.face.std = std(results.remember.merged.bar.NEW.precue.face.bytrial);
        
        results_rem_merged_new_precue_acts_scene = nanmean(results.remember.merged.acts.NEW.scene(:,1:3),2); %merging 1-3TR
        results.remember.merged.bar.NEW.precue.scene.bytrial =  results_rem_merged_new_precue_acts_scene;
        results.remember.merged.bar.NEW.precue.scene.nanmean = nanmean(results.remember.merged.bar.NEW.precue.scene.bytrial);
        results.remember.merged.bar.NEW.precue.scene.std = std(results.remember.merged.bar.NEW.precue.scene.bytrial);
        
        results_rem_merged_new_precue_acts_object = nanmean(results.remember.merged.acts.NEW.object(:,1:3),2); %merging 1-3TR
        results.remember.merged.bar.NEW.precue.object.bytrial =  results_rem_merged_new_precue_acts_object;
        results.remember.merged.bar.NEW.precue.object.nanmean = nanmean(results.remember.merged.bar.NEW.precue.object.bytrial);
        results.remember.merged.bar.NEW.precue.object.std = std(results.remember.merged.bar.NEW.precue.object.bytrial);
               
        results_rem_merged_new_precue_acts_rest = nanmean(results.remember.merged.acts.NEW.rest(:,1:3),2); %merging 1-3TR
        results.remember.merged.bar.NEW.precue.rest.bytrial =  results_rem_merged_new_precue_acts_rest;
        results.remember.merged.bar.NEW.precue.rest.nanmean = nanmean(results.remember.merged.bar.NEW.precue.rest.bytrial);
        results.remember.merged.bar.NEW.precue.rest.std = std(results.remember.merged.bar.NEW.precue.rest.bytrial);
        
        %
        
        %Remember - MERGED - NEW- POSTCUE
        
        results_rem_merged_new_postcue_acts_face = nanmean(results.remember.merged.acts.NEW.face(:,5:7),2); %merging 1-3TR
        results.remember.merged.bar.NEW.postcue.face.bytrial =  results_rem_merged_new_postcue_acts_face;
        results.remember.merged.bar.NEW.postcue.face.nanmean = nanmean(results.remember.merged.bar.NEW.postcue.face.bytrial);
        results.remember.merged.bar.NEW.postcue.face.std = std(results.remember.merged.bar.NEW.postcue.face.bytrial);
        
        results_rem_merged_new_postcue_acts_scene = nanmean(results.remember.merged.acts.NEW.scene(:,5:7),2); %merging 1-3TR
        results.remember.merged.bar.NEW.postcue.scene.bytrial =  results_rem_merged_new_postcue_acts_scene;
        results.remember.merged.bar.NEW.postcue.scene.nanmean = nanmean(results.remember.merged.bar.NEW.postcue.scene.bytrial);
        results.remember.merged.bar.NEW.postcue.scene.std = std(results.remember.merged.bar.NEW.postcue.scene.bytrial);
        
        results_rem_merged_new_postcue_acts_object = nanmean(results.remember.merged.acts.NEW.object(:,5:7),2); %merging 1-3TR
        results.remember.merged.bar.NEW.postcue.object.bytrial =  results_rem_merged_new_postcue_acts_object;
        results.remember.merged.bar.NEW.postcue.object.nanmean = nanmean(results.remember.merged.bar.NEW.postcue.object.bytrial);
        results.remember.merged.bar.NEW.postcue.object.std = std(results.remember.merged.bar.NEW.postcue.object.bytrial);
               
        results_rem_merged_new_postcue_acts_rest = nanmean(results.remember.merged.acts.NEW.rest(:,5:7),2); %merging 1-3TR
        results.remember.merged.bar.NEW.postcue.rest.bytrial =  results_rem_merged_new_postcue_acts_rest;
        results.remember.merged.bar.NEW.postcue.rest.nanmean = nanmean(results.remember.merged.bar.NEW.postcue.rest.bytrial);
        results.remember.merged.bar.NEW.postcue.rest.std = std(results.remember.merged.bar.NEW.postcue.rest.bytrial);
    end
    
    %%%%%%%%% HIGH CONF %%%%%%%%%%%%%
  
          %Remember - HICONF - OLD - PRECUE
    results_rem_highconf_old_precue_acts_face = nanmean(results.remember.highconf.acts.OLD.face(:,1:3),2); %merging 1-3TR
    results.remember.highconf.bar.OLD.precue.face.bytrial =  results_rem_highconf_old_precue_acts_face;
    results.remember.highconf.bar.OLD.precue.face.nanmean = nanmean(results.remember.highconf.bar.OLD.precue.face.bytrial);
    results.remember.highconf.bar.OLD.precue.face.std = std(results.remember.highconf.bar.OLD.precue.face.bytrial);
    
    results_rem_highconf_old_precue_acts_scene = nanmean(results.remember.highconf.acts.OLD.scene(:,1:3),2); %merging 1-3TR
    results.remember.highconf.bar.OLD.precue.scene.bytrial =  results_rem_highconf_old_precue_acts_scene;
    results.remember.highconf.bar.OLD.precue.scene.nanmean = nanmean(results.remember.highconf.bar.OLD.precue.scene.bytrial);
    results.remember.highconf.bar.OLD.precue.scene.std = std(results.remember.highconf.bar.OLD.precue.scene.bytrial);
    
    results_rem_highconf_old_precue_acts_object = nanmean(results.remember.highconf.acts.OLD.object(:,1:3),2); %merging 1-3TR
    results.remember.highconf.bar.OLD.precue.object.bytrial =  results_rem_highconf_old_precue_acts_object;
    results.remember.highconf.bar.OLD.precue.object.nanmean = nanmean(results.remember.highconf.bar.OLD.precue.object.bytrial);
    results.remember.highconf.bar.OLD.precue.object.std = std(results.remember.highconf.bar.OLD.precue.object.bytrial);
     
    results_rem_highconf_old_precue_acts_rest = nanmean(results.remember.highconf.acts.OLD.rest(:,1:3),2); %merging 1-3TR
    results.remember.highconf.bar.OLD.precue.rest.bytrial =  results_rem_highconf_old_precue_acts_rest;
    results.remember.highconf.bar.OLD.precue.rest.nanmean = nanmean(results.remember.highconf.bar.OLD.precue.rest.bytrial);
    results.remember.highconf.bar.OLD.precue.rest.std = std(results.remember.highconf.bar.OLD.precue.rest.bytrial);
  % 
  
  %Remember - highconf - OLD - POSTCUE
  
  results_rem_highconf_old_postcue_acts_face = nanmean(results.remember.highconf.acts.OLD.face(:,5:7),2); %merging 1-3TR
  results.remember.highconf.bar.OLD.postcue.face.bytrial =  results_rem_highconf_old_postcue_acts_face;
  results.remember.highconf.bar.OLD.postcue.face.nanmean = nanmean(results.remember.highconf.bar.OLD.postcue.face.bytrial);
  results.remember.highconf.bar.OLD.postcue.face.std = std(results.remember.highconf.bar.OLD.postcue.face.bytrial);
  
  results_rem_highconf_old_postcue_acts_scene = nanmean(results.remember.highconf.acts.OLD.scene(:,5:7),2); %merging 1-3TR
  results.remember.highconf.bar.OLD.postcue.scene.bytrial =  results_rem_highconf_old_postcue_acts_scene;
  results.remember.highconf.bar.OLD.postcue.scene.nanmean = nanmean(results.remember.highconf.bar.OLD.postcue.scene.bytrial);
  results.remember.highconf.bar.OLD.postcue.scene.std = std(results.remember.highconf.bar.OLD.postcue.scene.bytrial);
  
  results_rem_highconf_old_postcue_acts_object = nanmean(results.remember.highconf.acts.OLD.object(:,5:7),2); %merging 1-3TR
  results.remember.highconf.bar.OLD.postcue.object.bytrial =  results_rem_highconf_old_postcue_acts_object;
  results.remember.highconf.bar.OLD.postcue.object.nanmean = nanmean(results.remember.highconf.bar.OLD.postcue.object.bytrial);
  results.remember.highconf.bar.OLD.postcue.object.std = std(results.remember.highconf.bar.OLD.postcue.object.bytrial);
    
  results_rem_highconf_old_postcue_acts_rest = nanmean(results.remember.highconf.acts.OLD.rest(:,5:7),2); %merging 1-3TR
  results.remember.highconf.bar.OLD.postcue.rest.bytrial =  results_rem_highconf_old_postcue_acts_rest;
  results.remember.highconf.bar.OLD.postcue.rest.nanmean = nanmean(results.remember.highconf.bar.OLD.postcue.rest.bytrial);
  results.remember.highconf.bar.OLD.postcue.rest.std = std(results.remember.highconf.bar.OLD.postcue.rest.bytrial);
  
  if subjNum == '103' | subjNum == '113'
      %Remember - highconf - NEW- PRECUE
      
      results.remember.highconf.bar.NEW.precue.face.bytrial =  NaN;
      results.remember.highconf.bar.NEW.precue.face.nanmean = NaN;
      results.remember.highconf.bar.NEW.precue.face.std = NaN;
      
      
      results.remember.highconf.bar.NEW.precue.scene.bytrial = NaN;
      results.remember.highconf.bar.NEW.precue.scene.nanmean = NaN;
      results.remember.highconf.bar.NEW.precue.scene.std = NaN;
      
      
      results.remember.highconf.bar.NEW.precue.object.bytrial =  NaN;
      results.remember.highconf.bar.NEW.precue.object.nanmean = NaN;
      results.remember.highconf.bar.NEW.precue.object.std = NaN;
        
      
      results.remember.highconf.bar.NEW.precue.rest.bytrial =  NaN;
      results.remember.highconf.bar.NEW.precue.rest.nanmean = NaN;
      results.remember.highconf.bar.NEW.precue.rest.std = NaN;
      %
      
      %Remember - highconf - NEW- POSTCUE
      
      
      results.remember.highconf.bar.NEW.postcue.face.bytrial =  NaN;
      results.remember.highconf.bar.NEW.postcue.face.nanmean = NaN;
      results.remember.highconf.bar.NEW.postcue.face.std = NaN;
      
      
      results.remember.highconf.bar.NEW.postcue.scene.bytrial =  NaN;
      results.remember.highconf.bar.NEW.postcue.scene.nanmean = NaN;
      results.remember.highconf.bar.NEW.postcue.scene.std = NaN;
      
      
      results.remember.highconf.bar.NEW.postcue.object.bytrial =  NaN;
      results.remember.highconf.bar.NEW.postcue.object.nanmean = NaN;
      results.remember.highconf.bar.NEW.postcue.object.std = NaN;
         
      
      results.remember.highconf.bar.NEW.postcue.rest.bytrial = NaN;
      results.remember.highconf.bar.NEW.postcue.rest.nanmean = NaN;
      results.remember.highconf.bar.NEW.postcue.rest.std = NaN;
  else
      
      %Remember - highconf - NEW- PRECUE
      results_rem_highconf_new_precue_acts_face = nanmean(results.remember.highconf.acts.NEW.face(:,1:3),2); %merging 1-3TR
      results.remember.highconf.bar.NEW.precue.face.bytrial =  results_rem_highconf_new_precue_acts_face;
      results.remember.highconf.bar.NEW.precue.face.nanmean = nanmean(results.remember.highconf.bar.NEW.precue.face.bytrial);
      results.remember.highconf.bar.NEW.precue.face.std = std(results.remember.highconf.bar.NEW.precue.face.bytrial);
      
      results_rem_highconf_new_precue_acts_scene = nanmean(results.remember.highconf.acts.NEW.scene(:,1:3),2); %merging 1-3TR
      results.remember.highconf.bar.NEW.precue.scene.bytrial =  results_rem_highconf_new_precue_acts_scene;
      results.remember.highconf.bar.NEW.precue.scene.nanmean = nanmean(results.remember.highconf.bar.NEW.precue.scene.bytrial);
      results.remember.highconf.bar.NEW.precue.scene.std = std(results.remember.highconf.bar.NEW.precue.scene.bytrial);
      
      results_rem_highconf_new_precue_acts_object = nanmean(results.remember.highconf.acts.NEW.object(:,1:3),2); %merging 1-3TR
      results.remember.highconf.bar.NEW.precue.object.bytrial =  results_rem_highconf_new_precue_acts_object;
      results.remember.highconf.bar.NEW.precue.object.nanmean = nanmean(results.remember.highconf.bar.NEW.precue.object.bytrial);
      results.remember.highconf.bar.NEW.precue.object.std = std(results.remember.highconf.bar.NEW.precue.object.bytrial);
        
      results_rem_highconf_new_precue_acts_rest = nanmean(results.remember.highconf.acts.NEW.rest(:,1:3),2); %merging 1-3TR
      results.remember.highconf.bar.NEW.precue.rest.bytrial =  results_rem_highconf_new_precue_acts_rest;
      results.remember.highconf.bar.NEW.precue.rest.nanmean = nanmean(results.remember.highconf.bar.NEW.precue.rest.bytrial);
      results.remember.highconf.bar.NEW.precue.rest.std = std(results.remember.highconf.bar.NEW.precue.rest.bytrial);
      %
      
      %Remember - highconf - NEW- POSTCUE
      
      results_rem_highconf_new_postcue_acts_face = nanmean(results.remember.highconf.acts.NEW.face(:,5:7),2); %merging 1-3TR
      results.remember.highconf.bar.NEW.postcue.face.bytrial =  results_rem_highconf_new_postcue_acts_face;
      results.remember.highconf.bar.NEW.postcue.face.nanmean = nanmean(results.remember.highconf.bar.NEW.postcue.face.bytrial);
      results.remember.highconf.bar.NEW.postcue.face.std = std(results.remember.highconf.bar.NEW.postcue.face.bytrial);
      
      results_rem_highconf_new_postcue_acts_scene = nanmean(results.remember.highconf.acts.NEW.scene(:,5:7),2); %merging 1-3TR
      results.remember.highconf.bar.NEW.postcue.scene.bytrial =  results_rem_highconf_new_postcue_acts_scene;
      results.remember.highconf.bar.NEW.postcue.scene.nanmean = nanmean(results.remember.highconf.bar.NEW.postcue.scene.bytrial);
      results.remember.highconf.bar.NEW.postcue.scene.std = std(results.remember.highconf.bar.NEW.postcue.scene.bytrial);
      
      results_rem_highconf_new_postcue_acts_object = nanmean(results.remember.highconf.acts.NEW.object(:,5:7),2); %merging 1-3TR
      results.remember.highconf.bar.NEW.postcue.object.bytrial =  results_rem_highconf_new_postcue_acts_object;
      results.remember.highconf.bar.NEW.postcue.object.nanmean = nanmean(results.remember.highconf.bar.NEW.postcue.object.bytrial);
      results.remember.highconf.bar.NEW.postcue.object.std = std(results.remember.highconf.bar.NEW.postcue.object.bytrial);
         
      results_rem_highconf_new_postcue_acts_rest = nanmean(results.remember.highconf.acts.NEW.rest(:,5:7),2); %merging 1-3TR
      results.remember.highconf.bar.NEW.postcue.rest.bytrial =  results_rem_highconf_new_postcue_acts_rest;
      results.remember.highconf.bar.NEW.postcue.rest.nanmean = nanmean(results.remember.highconf.bar.NEW.postcue.rest.bytrial);
      results.remember.highconf.bar.NEW.postcue.rest.std = std(results.remember.highconf.bar.NEW.postcue.rest.bytrial);
  end
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% FORGET
    
        %forget - MERGED - OLD - PRECUE
    results_for_merged_old_precue_acts_face = nanmean(results.forget.merged.acts.OLD.face(:,1:3),2); %merging 1-3TR
    results.forget.merged.bar.OLD.precue.face.bytrial =  results_for_merged_old_precue_acts_face;
    results.forget.merged.bar.OLD.precue.face.nanmean = nanmean(results.forget.merged.bar.OLD.precue.face.bytrial);
    results.forget.merged.bar.OLD.precue.face.std = std(results.forget.merged.bar.OLD.precue.face.bytrial);
    
    results_for_merged_old_precue_acts_scene = nanmean(results.forget.merged.acts.OLD.scene(:,1:3),2); %merging 1-3TR
    results.forget.merged.bar.OLD.precue.scene.bytrial =  results_for_merged_old_precue_acts_scene;
    results.forget.merged.bar.OLD.precue.scene.nanmean = nanmean(results.forget.merged.bar.OLD.precue.scene.bytrial);
    results.forget.merged.bar.OLD.precue.scene.std = std(results.forget.merged.bar.OLD.precue.scene.bytrial);
    
    results_for_merged_old_precue_acts_object = nanmean(results.forget.merged.acts.OLD.object(:,1:3),2); %merging 1-3TR
    results.forget.merged.bar.OLD.precue.object.bytrial =  results_for_merged_old_precue_acts_object;
    results.forget.merged.bar.OLD.precue.object.nanmean = nanmean(results.forget.merged.bar.OLD.precue.object.bytrial);
    results.forget.merged.bar.OLD.precue.object.std = std(results.forget.merged.bar.OLD.precue.object.bytrial);
     
    results_for_merged_old_precue_acts_rest = nanmean(results.forget.merged.acts.OLD.rest(:,1:3),2); %merging 1-3TR
    results.forget.merged.bar.OLD.precue.rest.bytrial =  results_for_merged_old_precue_acts_rest;
    results.forget.merged.bar.OLD.precue.rest.nanmean = nanmean(results.forget.merged.bar.OLD.precue.rest.bytrial);
    results.forget.merged.bar.OLD.precue.rest.std = std(results.forget.merged.bar.OLD.precue.rest.bytrial);
  % 
  
  %forget - MERGED - OLD - POSTCUE
  
    results_for_merged_old_postcue_acts_face = nanmean(results.forget.merged.acts.OLD.face(:,5:7),2); %merging 1-3TR
    results.forget.merged.bar.OLD.postcue.face.bytrial =  results_for_merged_old_postcue_acts_face;
    results.forget.merged.bar.OLD.postcue.face.nanmean = nanmean(results.forget.merged.bar.OLD.postcue.face.bytrial);
    results.forget.merged.bar.OLD.postcue.face.std = std(results.forget.merged.bar.OLD.postcue.face.bytrial);
    
    results_for_merged_old_postcue_acts_scene = nanmean(results.forget.merged.acts.OLD.scene(:,5:7),2); %merging 1-3TR
    results.forget.merged.bar.OLD.postcue.scene.bytrial =  results_for_merged_old_postcue_acts_scene;
    results.forget.merged.bar.OLD.postcue.scene.nanmean = nanmean(results.forget.merged.bar.OLD.postcue.scene.bytrial);
    results.forget.merged.bar.OLD.postcue.scene.std = std(results.forget.merged.bar.OLD.postcue.scene.bytrial);
    
    results_for_merged_old_postcue_acts_object = nanmean(results.forget.merged.acts.OLD.object(:,5:7),2); %merging 1-3TR
    results.forget.merged.bar.OLD.postcue.object.bytrial =  results_for_merged_old_postcue_acts_object;
    results.forget.merged.bar.OLD.postcue.object.nanmean = nanmean(results.forget.merged.bar.OLD.postcue.object.bytrial);
    results.forget.merged.bar.OLD.postcue.object.std = std(results.forget.merged.bar.OLD.postcue.object.bytrial);
   
    results_for_merged_old_postcue_acts_rest = nanmean(results.forget.merged.acts.OLD.rest(:,5:7),2); %merging 1-3TR
    results.forget.merged.bar.OLD.postcue.rest.bytrial =  results_for_merged_old_postcue_acts_rest;
    results.forget.merged.bar.OLD.postcue.rest.nanmean = nanmean(results.forget.merged.bar.OLD.postcue.rest.bytrial);
    results.forget.merged.bar.OLD.postcue.rest.std = std(results.forget.merged.bar.OLD.postcue.rest.bytrial);
    
        %forget - MERGED - NEW- PRECUE
        
        if subjNum == '105' | subjNum == '113'
            
    results_for_merged_new_precue_acts_face = NaN;
    results.forget.merged.bar.NEW.precue.face.bytrial =  NaN;
    results.forget.merged.bar.NEW.precue.face.nanmean = NaN;
    results.forget.merged.bar.NEW.precue.face.std = NaN;
    
    results_for_merged_new_precue_acts_scene =NaN;
    results.forget.merged.bar.NEW.precue.scene.bytrial = NaN;
    results.forget.merged.bar.NEW.precue.scene.nanmean = NaN;
    results.forget.merged.bar.NEW.precue.scene.std = NaN;
    
    results_for_merged_new_precue_acts_object = NaN;
    results.forget.merged.bar.NEW.precue.object.bytrial =  NaN;
    results.forget.merged.bar.NEW.precue.object.nanmean = NaN;
    results.forget.merged.bar.NEW.precue.object.std = NaN;
      
    results_for_merged_new_precue_acts_rest =NaN;
    results.forget.merged.bar.NEW.precue.rest.bytrial =  NaN;
    results.forget.merged.bar.NEW.precue.rest.nanmean = NaN;
    results.forget.merged.bar.NEW.precue.rest.std = NaN;
  % 
  
  %forget - MERGED - NEW- POSTCUE
  
    results_for_merged_new_postcue_acts_face = NaN;
    results.forget.merged.bar.NEW.postcue.face.bytrial = NaN;
    results.forget.merged.bar.NEW.postcue.face.nanmean = NaN;
    results.forget.merged.bar.NEW.postcue.face.std = NaN;
    
    results_for_merged_new_postcue_acts_scene= NaN;
    results.forget.merged.bar.NEW.postcue.scene.bytrial = NaN;
    results.forget.merged.bar.NEW.postcue.scene.nanmean = NaN;
    results.forget.merged.bar.NEW.postcue.scene.std = NaN;
    
    results_for_merged_new_postcue_acts_object = NaN;
    results.forget.merged.bar.NEW.postcue.object.bytrial = NaN;
    results.forget.merged.bar.NEW.postcue.object.nanmean = NaN;
    results.forget.merged.bar.NEW.postcue.object.std = NaN;
        
    results_for_merged_new_postcue_acts_rest= NaN;
    results.forget.merged.bar.NEW.postcue.rest.bytrial = NaN;
    results.forget.merged.bar.NEW.postcue.rest.nanmean = NaN;
    results.forget.merged.bar.NEW.postcue.rest.std = NaN;
        
        else
    results_for_merged_new_precue_acts_face = nanmean(results.forget.merged.acts.NEW.face(:,1:3),2); %merging 1-3TR
    results.forget.merged.bar.NEW.precue.face.bytrial =  results_for_merged_new_precue_acts_face;
    results.forget.merged.bar.NEW.precue.face.nanmean = nanmean(results.forget.merged.bar.NEW.precue.face.bytrial);
    results.forget.merged.bar.NEW.precue.face.std = std(results.forget.merged.bar.NEW.precue.face.bytrial);
    
    results_for_merged_new_precue_acts_scene = nanmean(results.forget.merged.acts.NEW.scene(:,1:3),2); %merging 1-3TR
    results.forget.merged.bar.NEW.precue.scene.bytrial =  results_for_merged_new_precue_acts_scene;
    results.forget.merged.bar.NEW.precue.scene.nanmean = nanmean(results.forget.merged.bar.NEW.precue.scene.bytrial);
    results.forget.merged.bar.NEW.precue.scene.std = std(results.forget.merged.bar.NEW.precue.scene.bytrial);
    
    results_for_merged_new_precue_acts_object = nanmean(results.forget.merged.acts.NEW.object(:,1:3),2); %merging 1-3TR
    results.forget.merged.bar.NEW.precue.object.bytrial =  results_for_merged_new_precue_acts_object;
    results.forget.merged.bar.NEW.precue.object.nanmean = nanmean(results.forget.merged.bar.NEW.precue.object.bytrial);
    results.forget.merged.bar.NEW.precue.object.std = std(results.forget.merged.bar.NEW.precue.object.bytrial);
       
    results_for_merged_new_precue_acts_rest = nanmean(results.forget.merged.acts.NEW.rest(:,1:3),2); %merging 1-3TR
    results.forget.merged.bar.NEW.precue.rest.bytrial =  results_for_merged_new_precue_acts_rest;
    results.forget.merged.bar.NEW.precue.rest.nanmean = nanmean(results.forget.merged.bar.NEW.precue.rest.bytrial);
    results.forget.merged.bar.NEW.precue.rest.std = std(results.forget.merged.bar.NEW.precue.rest.bytrial);
  % 
  
  %forget - MERGED - NEW- POSTCUE
  
    results_for_merged_new_postcue_acts_face = nanmean(results.forget.merged.acts.NEW.face(:,5:7),2); %merging 5-7
    results.forget.merged.bar.NEW.postcue.face.bytrial =  results_for_merged_new_postcue_acts_face;
    results.forget.merged.bar.NEW.postcue.face.nanmean = nanmean(results.forget.merged.bar.NEW.postcue.face.bytrial);
    results.forget.merged.bar.NEW.postcue.face.std = std(results.forget.merged.bar.NEW.postcue.face.bytrial);
    
    results_for_merged_new_postcue_acts_scene = nanmean(results.forget.merged.acts.NEW.scene(:,5:7),2); %merging 1-3TR
    results.forget.merged.bar.NEW.postcue.scene.bytrial =  results_for_merged_new_postcue_acts_scene;
    results.forget.merged.bar.NEW.postcue.scene.nanmean = nanmean(results.forget.merged.bar.NEW.postcue.scene.bytrial);
    results.forget.merged.bar.NEW.postcue.scene.std = std(results.forget.merged.bar.NEW.postcue.scene.bytrial);
    
    results_for_merged_new_postcue_acts_object = nanmean(results.forget.merged.acts.NEW.object(:,5:7),2); %merging 1-3TR
    results.forget.merged.bar.NEW.postcue.object.bytrial =  results_for_merged_new_postcue_acts_object;
    results.forget.merged.bar.NEW.postcue.object.nanmean = nanmean(results.forget.merged.bar.NEW.postcue.object.bytrial);
    results.forget.merged.bar.NEW.postcue.object.std = std(results.forget.merged.bar.NEW.postcue.object.bytrial);
       
    results_for_merged_new_postcue_acts_rest = nanmean(results.forget.merged.acts.NEW.rest(:,5:7),2); %merging 1-3TR
    results.forget.merged.bar.NEW.postcue.rest.bytrial =  results_for_merged_new_postcue_acts_rest;
    results.forget.merged.bar.NEW.postcue.rest.nanmean = nanmean(results.forget.merged.bar.NEW.postcue.rest.bytrial);
    results.forget.merged.bar.NEW.postcue.rest.std = std(results.forget.merged.bar.NEW.postcue.rest.bytrial);
        end
    
    %%%%%%%%% HIGH CONF %%%%%%%%%%%%%
  
          %forget - HICONF - OLD - PRECUE
    results_for_highconf_old_precue_acts_face = nanmean(results.forget.highconf.acts.OLD.face(:,1:3),2); %merging 1-3TR
    results.forget.highconf.bar.OLD.precue.face.bytrial =  results_for_highconf_old_precue_acts_face;
    results.forget.highconf.bar.OLD.precue.face.nanmean = nanmean(results.forget.highconf.bar.OLD.precue.face.bytrial);
    results.forget.highconf.bar.OLD.precue.face.std = std(results.forget.highconf.bar.OLD.precue.face.bytrial);
    
    results_for_highconf_old_precue_acts_scene = nanmean(results.forget.highconf.acts.OLD.scene(:,1:3),2); %merging 1-3TR
    results.forget.highconf.bar.OLD.precue.scene.bytrial =  results_for_highconf_old_precue_acts_scene;
    results.forget.highconf.bar.OLD.precue.scene.nanmean = nanmean(results.forget.highconf.bar.OLD.precue.scene.bytrial);
    results.forget.highconf.bar.OLD.precue.scene.std = std(results.forget.highconf.bar.OLD.precue.scene.bytrial);
    
    results_for_highconf_old_precue_acts_object = nanmean(results.forget.highconf.acts.OLD.object(:,1:3),2); %merging 1-3TR
    results.forget.highconf.bar.OLD.precue.object.bytrial =  results_for_highconf_old_precue_acts_object;
    results.forget.highconf.bar.OLD.precue.object.nanmean = nanmean(results.forget.highconf.bar.OLD.precue.object.bytrial);
    results.forget.highconf.bar.OLD.precue.object.std = std(results.forget.highconf.bar.OLD.precue.object.bytrial);
       
    results_for_highconf_old_precue_acts_rest = nanmean(results.forget.highconf.acts.OLD.rest(:,1:3),2); %merging 1-3TR
    results.forget.highconf.bar.OLD.precue.rest.bytrial =  results_for_highconf_old_precue_acts_rest;
    results.forget.highconf.bar.OLD.precue.rest.nanmean = nanmean(results.forget.highconf.bar.OLD.precue.rest.bytrial);
    results.forget.highconf.bar.OLD.precue.rest.std = std(results.forget.highconf.bar.OLD.precue.rest.bytrial);
  % 
  
  %forget - highconf - OLD - POSTCUE
  
    results_for_highconf_old_postcue_acts_face = nanmean(results.forget.highconf.acts.OLD.face(:,5:7),2); %merging 1-3TR
    results.forget.highconf.bar.OLD.postcue.face.bytrial =  results_for_highconf_old_postcue_acts_face;
    results.forget.highconf.bar.OLD.postcue.face.nanmean = nanmean(results.forget.highconf.bar.OLD.postcue.face.bytrial);
    results.forget.highconf.bar.OLD.postcue.face.std = std(results.forget.highconf.bar.OLD.postcue.face.bytrial);
    
    results_for_highconf_old_postcue_acts_scene = nanmean(results.forget.highconf.acts.OLD.scene(:,5:7),2); %merging 1-3TR
    results.forget.highconf.bar.OLD.postcue.scene.bytrial =  results_for_highconf_old_postcue_acts_scene;
    results.forget.highconf.bar.OLD.postcue.scene.nanmean = nanmean(results.forget.highconf.bar.OLD.postcue.scene.bytrial);
    results.forget.highconf.bar.OLD.postcue.scene.std = std(results.forget.highconf.bar.OLD.postcue.scene.bytrial);
    
    results_for_highconf_old_postcue_acts_object = nanmean(results.forget.highconf.acts.OLD.object(:,5:7),2); %merging 1-3TR
    results.forget.highconf.bar.OLD.postcue.object.bytrial =  results_for_highconf_old_postcue_acts_object;
    results.forget.highconf.bar.OLD.postcue.object.nanmean = nanmean(results.forget.highconf.bar.OLD.postcue.object.bytrial);
    results.forget.highconf.bar.OLD.postcue.object.std = std(results.forget.highconf.bar.OLD.postcue.object.bytrial);
      
    results_for_highconf_old_postcue_acts_rest = nanmean(results.forget.highconf.acts.OLD.rest(:,5:7),2); %merging 1-3TR
    results.forget.highconf.bar.OLD.postcue.rest.bytrial =  results_for_highconf_old_postcue_acts_rest;
    results.forget.highconf.bar.OLD.postcue.rest.nanmean = nanmean(results.forget.highconf.bar.OLD.postcue.rest.bytrial);
    results.forget.highconf.bar.OLD.postcue.rest.std = std(results.forget.highconf.bar.OLD.postcue.rest.bytrial);
    
        %forget - highconf - NEW- PRECUE
        
         if subjNum == '105'
                
    results.forget.highconf.bar.NEW.precue.face.bytrial =  NaN;
    results.forget.highconf.bar.NEW.precue.face.nanmean = NaN;
    results.forget.highconf.bar.NEW.precue.face.std = NaN;
    
   
    results.forget.highconf.bar.NEW.precue.scene.bytrial =  NaN;
    results.forget.highconf.bar.NEW.precue.scene.nanmean = NaN;
    results.forget.highconf.bar.NEW.precue.scene.std = NaN;
    
   
    results.forget.highconf.bar.NEW.precue.object.bytrial =  NaN;
    results.forget.highconf.bar.NEW.precue.object.nanmean = NaN;
    results.forget.highconf.bar.NEW.precue.object.std = NaN;
    
    results.forget.highconf.bar.NEW.precue.rest.bytrial =  NaN;
    results.forget.highconf.bar.NEW.precue.rest.nanmean = NaN;
    results.forget.highconf.bar.NEW.precue.rest.std = NaN;
  % 
  
  %forget - highconf - NEW- POSTCUE
  

    results.forget.highconf.bar.NEW.postcue.face.bytrial =  NaN;
    results.forget.highconf.bar.NEW.postcue.face.nanmean = NaN;
    results.forget.highconf.bar.NEW.postcue.face.std = NaN;
    

    results.forget.highconf.bar.NEW.postcue.scene.bytrial =  NaN;
    results.forget.highconf.bar.NEW.postcue.scene.nanmean = NaN;
    results.forget.highconf.bar.NEW.postcue.scene.std = NaN;
    

    results.forget.highconf.bar.NEW.postcue.object.bytrial =  NaN;
    results.forget.highconf.bar.NEW.postcue.object.nanmean = NaN;
    results.forget.highconf.bar.NEW.postcue.object.std = NaN;

    results.forget.highconf.bar.NEW.postcue.rest.bytrial =  NaN;
    results.forget.highconf.bar.NEW.postcue.rest.nanmean = NaN;
    results.forget.highconf.bar.NEW.postcue.rest.std = NaN;
         else
             
    results_for_highconf_new_precue_acts_face = nanmean(results.forget.highconf.acts.NEW.face(:,1:3),2); %merging 1-3TR
    results.forget.highconf.bar.NEW.precue.face.bytrial =  results_for_highconf_new_precue_acts_face;
    results.forget.highconf.bar.NEW.precue.face.nanmean = nanmean(results.forget.highconf.bar.NEW.precue.face.bytrial);
    results.forget.highconf.bar.NEW.precue.face.std = std(results.forget.highconf.bar.NEW.precue.face.bytrial);
    
    results_for_highconf_new_precue_acts_scene = nanmean(results.forget.highconf.acts.NEW.scene(:,1:3),2); %merging 1-3TR
    results.forget.highconf.bar.NEW.precue.scene.bytrial =  results_for_highconf_new_precue_acts_scene;
    results.forget.highconf.bar.NEW.precue.scene.nanmean = nanmean(results.forget.highconf.bar.NEW.precue.scene.bytrial);
    results.forget.highconf.bar.NEW.precue.scene.std = std(results.forget.highconf.bar.NEW.precue.scene.bytrial);
    
    results_for_highconf_new_precue_acts_object = nanmean(results.forget.highconf.acts.NEW.object(:,1:3),2); %merging 1-3TR
    results.forget.highconf.bar.NEW.precue.object.bytrial =  results_for_highconf_new_precue_acts_object;
    results.forget.highconf.bar.NEW.precue.object.nanmean = nanmean(results.forget.highconf.bar.NEW.precue.object.bytrial);
    results.forget.highconf.bar.NEW.precue.object.std = std(results.forget.highconf.bar.NEW.precue.object.bytrial);
     
    results_for_highconf_new_precue_acts_rest = nanmean(results.forget.highconf.acts.NEW.rest(:,1:3),2); %merging 1-3TR
    results.forget.highconf.bar.NEW.precue.rest.bytrial =  results_for_highconf_new_precue_acts_rest;
    results.forget.highconf.bar.NEW.precue.rest.nanmean = nanmean(results.forget.highconf.bar.NEW.precue.rest.bytrial);
    results.forget.highconf.bar.NEW.precue.rest.std = std(results.forget.highconf.bar.NEW.precue.rest.bytrial);
  % 
  
  %forget - highconf - NEW- POSTCUE
  
    results_for_highconf_new_postcue_acts_face = nanmean(results.forget.highconf.acts.NEW.face(:,5:7),2); %merging 1-3TR
    results.forget.highconf.bar.NEW.postcue.face.bytrial =  results_for_highconf_new_postcue_acts_face;
    results.forget.highconf.bar.NEW.postcue.face.nanmean = nanmean(results.forget.highconf.bar.NEW.postcue.face.bytrial);
    results.forget.highconf.bar.NEW.postcue.face.std = std(results.forget.highconf.bar.NEW.postcue.face.bytrial);
    
    results_for_highconf_new_postcue_acts_scene = nanmean(results.forget.highconf.acts.NEW.scene(:,5:7),2); %merging 1-3TR
    results.forget.highconf.bar.NEW.postcue.scene.bytrial =  results_for_highconf_new_postcue_acts_scene;
    results.forget.highconf.bar.NEW.postcue.scene.nanmean = nanmean(results.forget.highconf.bar.NEW.postcue.scene.bytrial);
    results.forget.highconf.bar.NEW.postcue.scene.std = std(results.forget.highconf.bar.NEW.postcue.scene.bytrial);
    
    results_for_highconf_new_postcue_acts_object = nanmean(results.forget.highconf.acts.NEW.object(:,5:7),2); %merging 1-3TR
    results.forget.highconf.bar.NEW.postcue.object.bytrial =  results_for_highconf_new_postcue_acts_object;
    results.forget.highconf.bar.NEW.postcue.object.nanmean = nanmean(results.forget.highconf.bar.NEW.postcue.object.bytrial);
    results.forget.highconf.bar.NEW.postcue.object.std = std(results.forget.highconf.bar.NEW.postcue.object.bytrial);
    
    results_for_highconf_new_postcue_acts_rest = nanmean(results.forget.highconf.acts.NEW.rest(:,5:7),2); %merging 1-3TR
    results.forget.highconf.bar.NEW.postcue.rest.bytrial =  results_for_highconf_new_postcue_acts_rest;
    results.forget.highconf.bar.NEW.postcue.rest.nanmean = nanmean(results.forget.highconf.bar.NEW.postcue.rest.bytrial);
    results.forget.highconf.bar.NEW.postcue.rest.std = std(results.forget.highconf.bar.NEW.postcue.rest.bytrial);
         end
    
    %%%%%% BY EMOTION %%%%%
    
        %Remember - neutral- PRECUE
    results_rem_neu_precue_acts_face = nanmean(results.remember.neu.acts.face(:,1:3),2); %merging 1-3TR
    results.remember.neu.bar.precue.face.bytrial =  results_rem_neu_precue_acts_face;
    results.remember.neu.bar.precue.face.nanmean = nanmean(results.remember.neu.bar.precue.face.bytrial);
    results.remember.neu.bar.precue.face.std = std(results.remember.neu.bar.precue.face.bytrial);
    
    results_rem_neu_precue_acts_scene = nanmean(results.remember.neu.acts.scene(:,1:3),2); %merging 1-3TR
    results.remember.neu.bar.precue.scene.bytrial =  results_rem_neu_precue_acts_scene;
    results.remember.neu.bar.precue.scene.nanmean = nanmean(results.remember.neu.bar.precue.scene.bytrial);
    results.remember.neu.bar.precue.scene.std = std(results.remember.neu.bar.precue.scene.bytrial);
    
    results_rem_neu_precue_acts_object = nanmean(results.remember.neu.acts.object(:,1:3),2); %merging 1-3TR
    results.remember.neu.bar.precue.object.bytrial =  results_rem_neu_precue_acts_object;
    results.remember.neu.bar.precue.object.nanmean = nanmean(results.remember.neu.bar.precue.object.bytrial);
    results.remember.neu.bar.precue.object.std = std(results.remember.neu.bar.precue.object.bytrial);
        
    results_rem_neu_precue_acts_rest = nanmean(results.remember.neu.acts.rest(:,1:3),2); %merging 1-3TR
    results.remember.neu.bar.precue.rest.bytrial =  results_rem_neu_precue_acts_rest;
    results.remember.neu.bar.precue.rest.nanmean = nanmean(results.remember.neu.bar.precue.rest.bytrial);
    results.remember.neu.bar.precue.rest.std = std(results.remember.neu.bar.precue.rest.bytrial);
    %Remember - neutral- postcue
    results_rem_neu_postcue_acts_face = nanmean(results.remember.neu.acts.face(:,5:7),2); %merging 1-3TR
    results.remember.neu.bar.postcue.face.bytrial =  results_rem_neu_postcue_acts_face;
    results.remember.neu.bar.postcue.face.nanmean = nanmean(results.remember.neu.bar.postcue.face.bytrial);
    results.remember.neu.bar.postcue.face.std = std(results.remember.neu.bar.postcue.face.bytrial);
    
    results_rem_neu_postcue_acts_scene = nanmean(results.remember.neu.acts.scene(:,5:7),2); %merging 1-3TR
    results.remember.neu.bar.postcue.scene.bytrial =  results_rem_neu_postcue_acts_scene;
    results.remember.neu.bar.postcue.scene.nanmean = nanmean(results.remember.neu.bar.postcue.scene.bytrial);
    results.remember.neu.bar.postcue.scene.std = std(results.remember.neu.bar.postcue.scene.bytrial);
    
    results_rem_neu_postcue_acts_object = nanmean(results.remember.neu.acts.object(:,5:7),2); %merging 1-3TR
    results.remember.neu.bar.postcue.object.bytrial =  results_rem_neu_postcue_acts_object;
    results.remember.neu.bar.postcue.object.nanmean = nanmean(results.remember.neu.bar.postcue.object.bytrial);
    results.remember.neu.bar.postcue.object.std = std(results.remember.neu.bar.postcue.object.bytrial);
      
    results_rem_neu_postcue_acts_rest = nanmean(results.remember.neu.acts.rest(:,5:7),2); %merging 1-3TR
    results.remember.neu.bar.postcue.rest.bytrial =  results_rem_neu_postcue_acts_rest;
    results.remember.neu.bar.postcue.rest.nanmean = nanmean(results.remember.neu.bar.postcue.rest.bytrial);
    results.remember.neu.bar.postcue.rest.std = std(results.remember.neu.bar.postcue.rest.bytrial);
    
            %Remember -negative- PRECUE
    results_rem_neg_precue_acts_face = nanmean(results.remember.neg.acts.face(:,1:3),2); %merging 1-3TR
    results.remember.neg.bar.precue.face.bytrial =  results_rem_neg_precue_acts_face;
    results.remember.neg.bar.precue.face.nanmean = nanmean(results.remember.neg.bar.precue.face.bytrial);
    results.remember.neg.bar.precue.face.std = std(results.remember.neg.bar.precue.face.bytrial);
    
    results_rem_neg_precue_acts_scene = nanmean(results.remember.neg.acts.scene(:,1:3),2); %merging 1-3TR
    results.remember.neg.bar.precue.scene.bytrial =  results_rem_neg_precue_acts_scene;
    results.remember.neg.bar.precue.scene.nanmean = nanmean(results.remember.neg.bar.precue.scene.bytrial);
    results.remember.neg.bar.precue.scene.std = std(results.remember.neg.bar.precue.scene.bytrial);
    
    results_rem_neg_precue_acts_object = nanmean(results.remember.neg.acts.object(:,1:3),2); %merging 1-3TR
    results.remember.neg.bar.precue.object.bytrial =  results_rem_neg_precue_acts_object;
    results.remember.neg.bar.precue.object.nanmean = nanmean(results.remember.neg.bar.precue.object.bytrial);
    results.remember.neg.bar.precue.object.std = std(results.remember.neg.bar.precue.object.bytrial);
        
    results_rem_neg_precue_acts_rest = nanmean(results.remember.neg.acts.rest(:,1:3),2); %merging 1-3TR
    results.remember.neg.bar.precue.rest.bytrial =  results_rem_neg_precue_acts_rest;
    results.remember.neg.bar.precue.rest.nanmean = nanmean(results.remember.neg.bar.precue.rest.bytrial);
    results.remember.neg.bar.precue.rest.std = std(results.remember.neg.bar.precue.rest.bytrial);
    %Remember -negative- postcue
    results_rem_neg_postcue_acts_face = nanmean(results.remember.neg.acts.face(:,5:7),2); %merging 1-3TR
    results.remember.neg.bar.postcue.face.bytrial =  results_rem_neg_postcue_acts_face;
    results.remember.neg.bar.postcue.face.nanmean = nanmean(results.remember.neg.bar.postcue.face.bytrial);
    results.remember.neg.bar.postcue.face.std = std(results.remember.neg.bar.postcue.face.bytrial);
    
    results_rem_neg_postcue_acts_scene = nanmean(results.remember.neg.acts.scene(:,5:7),2); %merging 1-3TR
    results.remember.neg.bar.postcue.scene.bytrial =  results_rem_neg_postcue_acts_scene;
    results.remember.neg.bar.postcue.scene.nanmean = nanmean(results.remember.neg.bar.postcue.scene.bytrial);
    results.remember.neg.bar.postcue.scene.std = std(results.remember.neg.bar.postcue.scene.bytrial);
    
    results_rem_neg_postcue_acts_object = nanmean(results.remember.neg.acts.object(:,5:7),2); %merging 1-3TR
    results.remember.neg.bar.postcue.object.bytrial =  results_rem_neg_postcue_acts_object;
    results.remember.neg.bar.postcue.object.nanmean = nanmean(results.remember.neg.bar.postcue.object.bytrial);
    results.remember.neg.bar.postcue.object.std = std(results.remember.neg.bar.postcue.object.bytrial);
        
    results_rem_neg_postcue_acts_rest = nanmean(results.remember.neg.acts.rest(:,5:7),2); %merging 1-3TR
    results.remember.neg.bar.postcue.rest.bytrial =  results_rem_neg_postcue_acts_rest;
    results.remember.neg.bar.postcue.rest.nanmean = nanmean(results.remember.neg.bar.postcue.rest.bytrial);
    results.remember.neg.bar.postcue.rest.std = std(results.remember.neg.bar.postcue.rest.bytrial);
    
    %%%%% forget %%%%
     %forget - neutral- PRECUE
    results_for_neu_precue_acts_face = nanmean(results.forget.neu.acts.face(:,1:3),2); %merging 1-3TR
    results.forget.neu.bar.precue.face.bytrial =  results_for_neu_precue_acts_face;
    results.forget.neu.bar.precue.face.nanmean = nanmean(results.forget.neu.bar.precue.face.bytrial);
    results.forget.neu.bar.precue.face.std = std(results.forget.neu.bar.precue.face.bytrial);
    
    results_for_neu_precue_acts_scene = nanmean(results.forget.neu.acts.scene(:,1:3),2); %merging 1-3TR
    results.forget.neu.bar.precue.scene.bytrial =  results_for_neu_precue_acts_scene;
    results.forget.neu.bar.precue.scene.nanmean = nanmean(results.forget.neu.bar.precue.scene.bytrial);
    results.forget.neu.bar.precue.scene.std = std(results.forget.neu.bar.precue.scene.bytrial);
    
    results_for_neu_precue_acts_object = nanmean(results.forget.neu.acts.object(:,1:3),2); %merging 1-3TR
    results.forget.neu.bar.precue.object.bytrial =  results_for_neu_precue_acts_object;
    results.forget.neu.bar.precue.object.nanmean = nanmean(results.forget.neu.bar.precue.object.bytrial);
    results.forget.neu.bar.precue.object.std = std(results.forget.neu.bar.precue.object.bytrial);
     
    results_for_neu_precue_acts_rest = nanmean(results.forget.neu.acts.rest(:,1:3),2); %merging 1-3TR
    results.forget.neu.bar.precue.rest.bytrial =  results_for_neu_precue_acts_rest;
    results.forget.neu.bar.precue.rest.nanmean = nanmean(results.forget.neu.bar.precue.rest.bytrial);
    results.forget.neu.bar.precue.rest.std = std(results.forget.neu.bar.precue.rest.bytrial);
    %forget - neutral- postcue
    results_for_neu_postcue_acts_face = nanmean(results.forget.neu.acts.face(:,5:7),2); %merging 1-3TR
    results.forget.neu.bar.postcue.face.bytrial =  results_for_neu_postcue_acts_face;
    results.forget.neu.bar.postcue.face.nanmean = nanmean(results.forget.neu.bar.postcue.face.bytrial);
    results.forget.neu.bar.postcue.face.std = std(results.forget.neu.bar.postcue.face.bytrial);
    
    results_for_neu_postcue_acts_scene = nanmean(results.forget.neu.acts.scene(:,5:7),2); %merging 1-3TR
    results.forget.neu.bar.postcue.scene.bytrial =  results_for_neu_postcue_acts_scene;
    results.forget.neu.bar.postcue.scene.nanmean = nanmean(results.forget.neu.bar.postcue.scene.bytrial);
    results.forget.neu.bar.postcue.scene.std = std(results.forget.neu.bar.postcue.scene.bytrial);
    
    results_for_neu_postcue_acts_object = nanmean(results.forget.neu.acts.object(:,5:7),2); %merging 1-3TR
    results.forget.neu.bar.postcue.object.bytrial =  results_for_neu_postcue_acts_object;
    results.forget.neu.bar.postcue.object.nanmean = nanmean(results.forget.neu.bar.postcue.object.bytrial);
    results.forget.neu.bar.postcue.object.std = std(results.forget.neu.bar.postcue.object.bytrial);
     
    results_for_neu_postcue_acts_rest = nanmean(results.forget.neu.acts.rest(:,5:7),2); %merging 1-3TR
    results.forget.neu.bar.postcue.rest.bytrial =  results_for_neu_postcue_acts_rest;
    results.forget.neu.bar.postcue.rest.nanmean = nanmean(results.forget.neu.bar.postcue.rest.bytrial);
    results.forget.neu.bar.postcue.rest.std = std(results.forget.neu.bar.postcue.rest.bytrial);
    
            %forget -negative- PRECUE
    results_for_neg_precue_acts_face = nanmean(results.forget.neg.acts.face(:,1:3),2); %merging 1-3TR
    results.forget.neg.bar.precue.face.bytrial =  results_for_neg_precue_acts_face;
    results.forget.neg.bar.precue.face.nanmean = nanmean(results.forget.neg.bar.precue.face.bytrial);
    results.forget.neg.bar.precue.face.std = std(results.forget.neg.bar.precue.face.bytrial);
    
    results_for_neg_precue_acts_scene = nanmean(results.forget.neg.acts.scene(:,1:3),2); %merging 1-3TR
    results.forget.neg.bar.precue.scene.bytrial =  results_for_neg_precue_acts_scene;
    results.forget.neg.bar.precue.scene.nanmean = nanmean(results.forget.neg.bar.precue.scene.bytrial);
    results.forget.neg.bar.precue.scene.std = std(results.forget.neg.bar.precue.scene.bytrial);
    
    results_for_neg_precue_acts_object = nanmean(results.forget.neg.acts.object(:,1:3),2); %merging 1-3TR
    results.forget.neg.bar.precue.object.bytrial =  results_for_neg_precue_acts_object;
    results.forget.neg.bar.precue.object.nanmean = nanmean(results.forget.neg.bar.precue.object.bytrial);
    results.forget.neg.bar.precue.object.std = std(results.forget.neg.bar.precue.object.bytrial);
     
    results_for_neg_precue_acts_rest = nanmean(results.forget.neg.acts.rest(:,1:3),2); %merging 1-3TR
    results.forget.neg.bar.precue.rest.bytrial =  results_for_neg_precue_acts_rest;
    results.forget.neg.bar.precue.rest.nanmean = nanmean(results.forget.neg.bar.precue.rest.bytrial);
    results.forget.neg.bar.precue.rest.std = std(results.forget.neg.bar.precue.rest.bytrial);
    %forget -negative- postcue
    results_for_neg_postcue_acts_face = nanmean(results.forget.neg.acts.face(:,5:7),2); %merging 1-3TR
    results.forget.neg.bar.postcue.face.bytrial =  results_for_neg_postcue_acts_face;
    results.forget.neg.bar.postcue.face.nanmean = nanmean(results.forget.neg.bar.postcue.face.bytrial);
    results.forget.neg.bar.postcue.face.std = std(results.forget.neg.bar.postcue.face.bytrial);
    
    results_for_neg_postcue_acts_scene = nanmean(results.forget.neg.acts.scene(:,5:7),2); %merging 1-3TR
    results.forget.neg.bar.postcue.scene.bytrial =  results_for_neg_postcue_acts_scene;
    results.forget.neg.bar.postcue.scene.nanmean = nanmean(results.forget.neg.bar.postcue.scene.bytrial);
    results.forget.neg.bar.postcue.scene.std = std(results.forget.neg.bar.postcue.scene.bytrial);
    
    results_for_neg_postcue_acts_object = nanmean(results.forget.neg.acts.object(:,5:7),2); %merging 1-3TR
    results.forget.neg.bar.postcue.object.bytrial =  results_for_neg_postcue_acts_object;
    results.forget.neg.bar.postcue.object.nanmean = nanmean(results.forget.neg.bar.postcue.object.bytrial);
    results.forget.neg.bar.postcue.object.std = std(results.forget.neg.bar.postcue.object.bytrial);
       
    results_for_neg_postcue_acts_rest = nanmean(results.forget.neg.acts.rest(:,5:7),2); %merging 1-3TR
    results.forget.neg.bar.postcue.rest.bytrial =  results_for_neg_postcue_acts_rest;
    results.forget.neg.bar.postcue.rest.nanmean = nanmean(results.forget.neg.bar.postcue.rest.bytrial);
    results.forget.neg.bar.postcue.rest.std = std(results.forget.neg.bar.postcue.rest.bytrial);
    
  filename = sprintf('results_%s_bar_compare.mat',test_phase);
  save(filename,'results');

  %%%%%%%%%TBD extreme confidence
end
cd(args.script_dir);
end
  
 