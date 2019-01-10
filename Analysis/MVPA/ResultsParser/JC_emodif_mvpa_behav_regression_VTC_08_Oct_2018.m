% "Supersubject" logistics regression (categorical/dummy 1, 0) accuracy vs.
% MVAP PostCue evidence
% 
%
% 
% Input
% Requires "results_DFencode_bar_compare.mat' for each partcipant. 
%
% Output
% Mat file containing the summary tables
% Excel (or CSV on Macs) summary table
% 
%
% Judy Chiu
% Oct. 2018 
%
% corrected "Remember "postcue" mistake" -Dec 13. 2018. 


%% Before running whole script
clear all; close all; clc;

computeruse=2;  %  1 for mac, 2 for PC

%%  Input Arguments %%%%%%%

% code directory
codedirMAC='/Users/jc/Box Sync/CMF_BoxWork/Exp1_JudyMVPA/EmoDF_Scripts/EmoDF_LocGit/MVPA_Results';
codedirPC='C:\Users\CMF_JC\Box Sync\CMF_BoxWork\Exp1_JudyMVPA\EmoDF_Scripts\EmoDF_LocGit ' ;

% slash
slashMAC='/';
slashPC='\';

% Data directory
datadirMAC='/Users/jc/Box Sync/SharedFiles/EmoDF_LewPeaLab/emodif_data/';
datadirPC='C:\Users\CMF_JC\Box Sync\SharedFiles\EmoDF_LewPeaLab\emodif_data\';

% Output directory
outdirMAC='/Users/jc/EmoDFmac/Analysis/MVPA_Results';
outdirPC='C:\Users\CMF_JC\Box Sync\CMF_BoxWork\Exp1_JudyMVPA\Pilot_Programs\JudyfMRI_Pilot\MR_PilotData\Analyzed\MVPA_Results_Excel';

if computeruse==1
    codedir=codedirMAC;
    datadir=datadirMAC;
    slashsym=slashMAC;
    outdir=outdirMAC;
    
    %Subfolder directories for MAC
    phase1dir='localizer/';
    phase2dir='DFencode/';
    phase3dir='Preview/';
    mask1dir='tempoccfusi_pHg_LOC_combined_epi_space/';
    mask2dir='JC_Combine_PHc_epi_space/';
    mask3dir='JC_AMYHP_epi_space/';
    mask4dir='scene-all_005/';
    
    
    % Specify folder info
    datecell={'21-Aug-2018/','29-Aug-2018/','18-Sep-2018/','24-Sep-2018/','25-Sep-2018/'};
    
    
elseif computeruse==2
    codedir=codedirPC;
    datadir=datadirPC;
    slashsym=slashPC;
    outdir=outdirPC;
    %Subfolder directories for DELL
    phase1dir='localizer\';
    phase2dir='DFencode\';
    phase3dir='Preview\';
    mask1dir='tempoccfusi_pHg_LOC_combined_epi_space\';
    mask2dir='JC_Combine_PHc_epi_space\';
    mask3dir='JC_AMYHP_epi_space\';
    mask4dir='scene-all_005\';
    
    % Specify folder info
    % Specify folder info
    datecell={'21-Aug-2018\','29-Aug-2018\','18-Sep-2018\','24-Sep-2018\','25-Sep-2018\'};
    
    
else
    showerr=' Must specify anlaysis computer'
end


% subjects to run through
sublist=[101,102,103,104,105,106,107,108,109,110,111,113,114,115,116,117,118,119,120,121,122,123,124,125];
datedir_ind=[1,1,2,1,4,4,5,4,4,4,5,5,4,5,5,5,5,4,4,4,4,4,5,5];
phasedir=phase2dir;
maskdir=mask1dir;

% output filename
phasedirstr=phasedir(1:end-1);


if size(maskdir,2)==size(mask1dir,2)
    maskdirstr='VTC';
elseif size(maskdir,2)==size(mask2dir,2)
    maskdirstr='PHC';
elseif size(maskdir,2)==size(mask3dir,2)
    maskdirstr='AMY';
elseif size(maskdir,2)==size(mask4dir,2)
    maskdirstr='fPHC';
else
end;

outname=strcat('emodif_',phasedirstr,'_',maskdirstr,'_MVPA_Behav_Regression_',date,'.xlsx');

% output .MAT file name
matname=strcat('EmoDiF_mvpa_behav_regression_',maskdirstr,'_',date,'.mat');


%% Parameters for script

% row and column indices
% % % % % stim type - -allsub_data(:,%,:)
% % % % r_wordcol=1;
% % % % r_scenecol=2;
% % % % r_facecol=3;
% % % % r_objcol=4;
% % % % r_restcol=5;
% % % % f_wordcol=6;
% % % % f_scenecol=7;
% % % % f_facecol=8;
% % % % f_objcol=9;
% % % % f_restcol=10;
% % % % 

%Data indices supersub_RF (:,:,%); thissub_RF (:,%) 
totcols=14;

subcol=1;
trialcol=2;
catcol=3;
RFcol=4;
precuecol=5;
postcuecol=6;
postminprecol=7;
TR1col=8;
TR2col=9;
TR3col=10;
TR4col=11;
TR5col=12;
TR6col=13;
TR7col=14;

% % % % precuestd=14;
% % % % postcuestd=15;
% % % % postminprestd=16;
% % % % TR1std=17;
% % % % TR7std=23;


% preallocate all sub matrix with 1st subjects' dimensions
supersub_RF=nan(5,60,14);  % 5 categories,trials per subject 1, 13 cols of data

for sub=1:size(sublist,2)
    
    clear results;
    
    thissub_RF_word=nan(60,totcols);
    thissub_RF_scene=nan(60,totcols);
    thissub_RF_face=nan(60,totcols);
    thissub_RF_object=nan(60,totcols);
    thissub_RF_rest=nan(60,totcols);
    
    
    subdir=strcat('emodif_', num2str(sublist(sub)),slashsym);
    datedir=datecell{datedir_ind(sub)};
    
    % load results struct from the specified phase
    % specify date by subject
    loaddir=strcat(datadir,subdir,'results',slashsym,phasedir,maskdir,datedir)
    
    load(strcat(loaddir,'results_DFencode_bar_compare.mat'));
    
    % Participant and condition info
    % Word
    thissub_RF_word(:,subcol)=sublist(sub);
    thissub_RF_word(:,catcol)=[4]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
    thissub_RF_word(1:30,RFcol)=1;  %remember
    thissub_RF_word(31:60,RFcol)=2; % forget
    
    % Scene
    thissub_RF_scene(:,subcol)=sublist(sub);
    thissub_RF_scene(:,catcol)=[2]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
    thissub_RF_scene(1:30,RFcol)=1;  %remember
    thissub_RF_scene(31:60,RFcol)=2; % forget
    
    %Face
    thissub_RF_face(:,subcol)=sublist(sub);
    thissub_RF_face(:,catcol)=[1]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
    thissub_RF_face(1:30,RFcol)=1;  %remember
    thissub_RF_face(31:60,RFcol)=2; % forget
    
    %Object
    thissub_RF_object(:,subcol)=sublist(sub);
    thissub_RF_object(:,catcol)=[3]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
    thissub_RF_object(1:30,RFcol)=1;  %remember
    thissub_RF_object(31:60,RFcol)=2; % forget
    
    %Rest
    thissub_RF_rest(:,subcol)=sublist(sub);
    thissub_RF_rest(:,catcol)=[5]; % stim category dummy code 1=f, 2=s, 3=o, 4=w, 5=r;
    thissub_RF_rest(1:30,RFcol)=1;  %remember
    thissub_RF_rest(31:60,RFcol)=2; % forget
    
    for rr=1:60  %trials for sub
        
        
        
        %% By main effect of Remember vs. Forget %%%%%
        if rr<31
            
            %%%%%%%%%  By TR timecourse graphs ( 1 by 7 timepoints) %%%%%%%%%
            
            %trial number
            thissub_RF_word(rr,trialcol)=rr;
            thissub_RF_scene(rr,trialcol)=rr;
            thissub_RF_face(rr,trialcol)=rr;
            thissub_RF_object(rr,trialcol)=rr;
            thissub_RF_rest(rr,trialcol)=rr;
            
            % Remember
               
            thissub_RF_word(rr,TR1col:TR7col)=results.remember.all.acts.word(rr,:);
            thissub_RF_scene(rr,TR1col:TR7col)=results.remember.all.acts.scene(rr,:);
            thissub_RF_face(rr,TR1col:TR7col)=results.remember.all.acts.face(rr,:);
            thissub_RF_object(rr,TR1col:TR7col)=results.remember.all.acts.object(rr,:);
            thissub_RF_rest(rr,TR1col:TR7col)=results.remember.all.acts.rest(rr,:);
            
            
            %%%%%%%%%%  Averaged TR bar graphs ( 1 number per data point) %%%%%%%%%%%%
            
            % Remember Precue
            thissub_RF_word(rr,precuecol)=mean(results.remember.all.acts.word(rr,1:3));
            thissub_RF_scene(rr,precuecol)=mean(results.remember.all.acts.scene(rr,1:3));
            thissub_RF_face(rr,precuecol)=mean(results.remember.all.acts.face(rr,1:3));
            thissub_RF_object(rr,precuecol)=mean(results.remember.all.acts.object(rr,1:3));
            thissub_RF_rest(rr,precuecol)=mean(results.remember.all.acts.rest(rr,1:3));
            
            % Remember Postcue
            thissub_RF_word(rr,postcuecol)=mean(results.remember.all.acts.word(rr,5:7));
            thissub_RF_scene(rr,postcuecol)=mean(results.remember.all.acts.scene(rr,5:7));
            thissub_RF_face(rr,postcuecol)=mean(results.remember.all.acts.face(rr,5:7));
            thissub_RF_object(rr,postcuecol)=mean(results.remember.all.acts.object(rr,5:7));
            thissub_RF_rest(rr,postcuecol)=mean(results.remember.all.acts.rest(rr,5:7));
            
            
        else %if rr>30
            
            % adjust row reference variable
            rrf=rr-30;
            
            %trial number
            thissub_RF_word(rr,trialcol)=rrf;
            thissub_RF_scene(rr,trialcol)=rrf;
            thissub_RF_face(rr,trialcol)=rrf;
            thissub_RF_object(rr,trialcol)=rrf;
            thissub_RF_rest(rr,trialcol)=rrf;
            
            % Forget
            thissub_RF_word(rr,TR1col:TR7col)=results.forget.all.acts.word(rrf,:);
            thissub_RF_scene(rr,TR1col:TR7col)=results.forget.all.acts.scene(rrf,:);
            thissub_RF_face(rr,TR1col:TR7col)=results.forget.all.acts.face(rrf,:);
            thissub_RF_object(rr,TR1col:TR7col)=results.forget.all.acts.object(rrf,:);
            thissub_RF_rest(rr,TR1col:TR7col)=results.forget.all.acts.rest(rrf,:);
            
            
            %%%%%%%%%%  Averaged TR bar graphs ( 1 number per data point) %%%%%%%%%%%%
            
            % Forget Precue
            thissub_RF_word(rr,precuecol)=mean(results.forget.all.acts.word(rrf,1:3));
            thissub_RF_scene(rr,precuecol)=mean(results.forget.all.acts.scene(rrf,1:3));
            thissub_RF_face(rr,precuecol)=mean(results.forget.all.acts.face(rrf,1:3));
            thissub_RF_object(rr,precuecol)=mean(results.forget.all.acts.object(rrf,1:3));
            thissub_RF_rest(rr,precuecol)=mean(results.forget.all.acts.rest(rrf,1:3));
            
            % forget Postcue
            thissub_RF_word(rr,postcuecol)=mean(results.forget.all.acts.word(rrf,5:7));
            thissub_RF_scene(rr,postcuecol)=mean(results.forget.all.acts.scene(rrf,5:7));
            thissub_RF_face(rr,postcuecol)=mean(results.forget.all.acts.face(rrf,5:7));
            thissub_RF_object(rr,postcuecol)=mean(results.forget.all.acts.object(rrf,5:7));
            thissub_RF_rest(rr,postcuecol)=mean(results.forget.all.acts.rest(rrf,5:7));
           
        end
    end
    
    % Within Subject Post-PreCue
    thissub_RF_word(:,postminprecol)=thissub_RF_word(:,postcuecol)-thissub_RF_word(:,precuecol);
    thissub_RF_scene(:,postminprecol)=thissub_RF_scene(:,postcuecol)-thissub_RF_scene(:,precuecol);
    thissub_RF_face(:,postminprecol)=thissub_RF_face(:,postcuecol)-thissub_RF_face(:,precuecol);
    thissub_RF_object(:,postminprecol)=thissub_RF_object(:,postcuecol)-thissub_RF_object(:,precuecol);
    thissub_RF_rest(:,postminprecol)=thissub_RF_rest(:,postcuecol)-thissub_RF_rest(:,precuecol);
    
     
    superrows=[1+(60*(sub-1)):60+(60*(sub-1))];
    % Populate All Sub Matrix
    supersub_RF(1,superrows,:)=thissub_RF_word;
    supersub_RF(2,superrows,:)=thissub_RF_scene;
    supersub_RF(3,superrows,:)=thissub_RF_face;
    supersub_RF(4,superrows,:)=thissub_RF_object;
    supersub_RF(5,superrows,:)=thissub_RF_rest;
  
end


% output to excel
cd(outdir);

xlsheader={'subject','trialcount','Category_fsowr','Rem/For','Precue','Postcue','PostminPre','TR1','TR2','TR3','TR4','TR5','TR6','TR7'};
xlswrite(outname,xlsheader,'wordRF','A1');
xlswrite(outname,squeeze(supersub_RF(1,:,:)),'wordRF','A2');

xlswrite(outname,xlsheader,'sceneRF','A1');
xlswrite(outname,squeeze(supersub_RF(2,:,:)),'sceneRF','A2');


xlswrite(outname,xlsheader,'faceRF','A1');
xlswrite(outname,squeeze(supersub_RF(3,:,:)),'faceRF','A2');


xlswrite(outname,xlsheader,'objRF','A1');
xlswrite(outname,squeeze(supersub_RF(4,:,:)),'objRF','A2');

xlswrite(outname,xlsheader,'restRF','A1');
xlswrite(outname,squeeze(supersub_RF(5,:,:)),'restRF','A2');

% save data in .MAT file
save(matname,'supersub_RF','datadir','codedir','outdir','sublist','datecell','datedir_ind','xlsheader','outname','matname');
cd(codedir);




%%  By Interation of Emotion by RF  %%%%%%%%%%%%%%%%%%%
% % % % % % 
% % % % % % clear allsub_RF_Neg;
% % % % % % clear allsub_RF_Neu;
% % % % % % allsub_RF_Neg=nan(size(sublist,2),10,23);
% % % % % % allsub_RF_Neu=nan(size(sublist,2),10,23);
% % % % % % 
% % % % % % for sub=1:size(sublist,2)
% % % % % %     clear results;
% % % % % %     
% % % % % %     subdir=strcat('emodif_', num2str(sublist(sub)),slashsym);
% % % % % %     datedir=datecell{datedir_ind(sub)};
% % % % % %     
% % % % % %     % load results struct from the specified phase
% % % % % %     % specify date by subject
% % % % % %     loaddir=strcat(datadir,subdir,'results',slashsym,phasedir,maskdir,datedir)
% % % % % %     
% % % % % %     load(strcat(loaddir,'results_DFencode_bar_compare.mat'));
% % % % % %     
% % % % % %     % Participant and condition info
% % % % % %     allsub_RF_Neg(sub,:,1)=sublist(sub);
% % % % % %     allsub_RF_Neg(sub,1:5,2)=[4,2,1,3,5]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
% % % % % %     allsub_RF_Neg(sub,6:10,2)=[4,2,1,3,5];
% % % % % %     
% % % % % %     allsub_RF_Neu(sub,:,1)=sublist(sub);
% % % % % %     allsub_RF_Neu(sub,1:5,2)=[4,2,1,3,5]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
% % % % % %     allsub_RF_Neu(sub,6:10,2)=[4,2,1,3,5];
% % % % % %     
% % % % % %     allsub_RF_Neg(sub,1:5,3)=1;  %remember
% % % % % %     allsub_RF_Neg(sub,6:10,3)=2; % forget
% % % % % %     
% % % % % %     allsub_RF_Neu(sub,1:5,3)=1;  %remember
% % % % % %     allsub_RF_Neu(sub,6:10,3)=2; % forget
% % % % % %     
% % % % % %     
% % % % % %     %%%%%%%%%%  Averaged TR bar graphs ( 1 number per data point) %%%%%%%%%%%%
% % % % % %     
% % % % % %     %%%%%%% Negative %%%%%%%%%%%%%%%
% % % % % %     % Remember Precue
% % % % % %     allsub_RF_Neg(sub,r_wordcol,precuecol)=results.remember.neg.bar.precue.word.nanmean;
% % % % % %     allsub_RF_Neg(sub,r_scenecol,precuecol)=results.remember.neg.bar.precue.scene.nanmean;
% % % % % %     allsub_RF_Neg(sub,r_facecol,precuecol)=results.remember.neg.bar.precue.face.nanmean;
% % % % % %     allsub_RF_Neg(sub,r_objcol,precuecol)=results.remember.neg.bar.precue.object.nanmean;
% % % % % %     allsub_RF_Neg(sub,r_restcol,precuecol)=results.remember.neg.bar.precue.rest.nanmean;
% % % % % %     
% % % % % %     allsub_RF_Neg(sub,r_wordcol,precuestd)=results.remember.neg.bar.precue.word.std;
% % % % % %     allsub_RF_Neg(sub,r_scenecol,precuestd)=results.remember.neg.bar.precue.scene.std;
% % % % % %     allsub_RF_Neg(sub,r_facecol,precuestd)=results.remember.neg.bar.precue.face.std;
% % % % % %     allsub_RF_Neg(sub,r_objcol,precuestd)=results.remember.neg.bar.precue.object.std;
% % % % % %     allsub_RF_Neg(sub,r_restcol,precuestd)=results.remember.neg.bar.precue.rest.std;
% % % % % %     
% % % % % %     
% % % % % %     % Remember Postcue
% % % % % %     allsub_RF_Neg(sub,r_wordcol,postcuecol)=results.remember.neg.bar.postcue.word.nanmean;
% % % % % %     allsub_RF_Neg(sub,r_scenecol,postcuecol)=results.remember.neg.bar.postcue.scene.nanmean;
% % % % % %     allsub_RF_Neg(sub,r_facecol,postcuecol)=results.remember.neg.bar.postcue.face.nanmean;
% % % % % %     allsub_RF_Neg(sub,r_objcol,postcuecol)=results.remember.neg.bar.postcue.object.nanmean;
% % % % % %     allsub_RF_Neg(sub,r_restcol,postcuecol)=results.remember.neg.bar.postcue.rest.nanmean;
% % % % % %     
% % % % % %     allsub_RF_Neg(sub,r_wordcol,postcuestd)=results.remember.neg.bar.postcue.word.std;
% % % % % %     allsub_RF_Neg(sub,r_scenecol,postcuestd)=results.remember.neg.bar.postcue.scene.std;
% % % % % %     allsub_RF_Neg(sub,r_facecol,postcuestd)=results.remember.neg.bar.postcue.face.std;
% % % % % %     allsub_RF_Neg(sub,r_objcol,postcuestd)=results.remember.neg.bar.postcue.object.std;
% % % % % %     allsub_RF_Neg(sub,r_restcol,postcuestd)=results.remember.neg.bar.postcue.rest.std;
% % % % % %     
% % % % % %     
% % % % % %     % % % % %     % Remember Post-Pre
% % % % % %     % % % % %     allsub_RF_Neg(sub,r_wordcol,postminprecol)=results.remember.neg.bar.postminuspre.word.nanmean;
% % % % % %     % % % % %     allsub_RF_Neg(sub,r_scenecol,postminprecol)=results.remember.neg.bar.postminuspre.scene.nanmean;
% % % % % %     % % % % %
% % % % % %     % % % % %
% % % % % %     % % % % %     allsub_RF_Neg(sub,r_wordcol,postminprestd)=results.remember.neg.bar.postminuspre.word.std;
% % % % % %     % % % % %     allsub_RF_Neg(sub,r_scenecol,postminprestd)=results.remember.neg.bar.postminuspre.scene.std;
% % % % % %     % % % % %
% % % % % %     
% % % % % %     
% % % % % %     % Forget Precue
% % % % % %     allsub_RF_Neg(sub,f_wordcol,precuecol)=results.forget.neg.bar.precue.word.nanmean;
% % % % % %     allsub_RF_Neg(sub,f_scenecol,precuecol)=results.forget.neg.bar.precue.scene.nanmean;
% % % % % %     allsub_RF_Neg(sub,f_facecol,precuecol)=results.forget.neg.bar.precue.face.nanmean;
% % % % % %     allsub_RF_Neg(sub,f_objcol,precuecol)=results.forget.neg.bar.precue.object.nanmean;
% % % % % %     allsub_RF_Neg(sub,f_restcol,precuecol)=results.forget.neg.bar.precue.rest.nanmean;
% % % % % %     
% % % % % %     allsub_RF_Neg(sub,f_wordcol,precuestd)=results.forget.neg.bar.precue.word.std;
% % % % % %     allsub_RF_Neg(sub,f_scenecol,precuestd)=results.forget.neg.bar.precue.scene.std;
% % % % % %     allsub_RF_Neg(sub,f_facecol,precuestd)=results.forget.neg.bar.precue.face.std;
% % % % % %     allsub_RF_Neg(sub,f_objcol,precuestd)=results.forget.neg.bar.precue.object.std;
% % % % % %     allsub_RF_Neg(sub,f_restcol,precuestd)=results.forget.neg.bar.precue.rest.std;
% % % % % %     
% % % % % %     % Forget Postcue
% % % % % %     allsub_RF_Neg(sub,f_wordcol,postcuecol)=results.forget.neg.bar.postcue.word.nanmean;
% % % % % %     allsub_RF_Neg(sub,f_scenecol,postcuecol)=results.forget.neg.bar.postcue.scene.nanmean;
% % % % % %     allsub_RF_Neg(sub,f_facecol,postcuecol)=results.forget.neg.bar.postcue.face.nanmean;
% % % % % %     allsub_RF_Neg(sub,f_objcol,postcuecol)=results.forget.neg.bar.postcue.object.nanmean;
% % % % % %     allsub_RF_Neg(sub,f_restcol,postcuecol)=results.forget.neg.bar.postcue.rest.nanmean;
% % % % % %     
% % % % % %     
% % % % % %     allsub_RF_Neg(sub,f_wordcol,postcuestd)=results.forget.neg.bar.postcue.word.std;
% % % % % %     allsub_RF_Neg(sub,f_scenecol,postcuestd)=results.forget.neg.bar.postcue.scene.std;
% % % % % %     allsub_RF_Neg(sub,f_facecol,postcuestd)=results.forget.neg.bar.postcue.face.std;
% % % % % %     allsub_RF_Neg(sub,f_objcol,postcuestd)=results.forget.neg.bar.postcue.object.std;
% % % % % %     allsub_RF_Neg(sub,f_restcol,postcuestd)=results.forget.neg.bar.postcue.rest.std;
% % % % % %     
% % % % % %     % % % % %
% % % % % %     % % % % %
% % % % % %     % % % % %     % Forget Post-Pre
% % % % % %     % % % % %     allsub_RF_Neg(sub,f_wordcol,postminprecol)=results.forget.neg.bar.postminuspre.word.nanmean;
% % % % % %     % % % % %     allsub_RF_Neg(sub,f_scenecol,postminprecol)=results.forget.neg.bar.postminuspre.scene.nanmean;
% % % % % %     % % % % %
% % % % % %     % % % % %     allsub_RF_Neg(sub,f_wordcol,postminprestd)=results.forget.neg.bar.postminuspre.word.std;
% % % % % %     % % % % %     allsub_RF_Neg(sub,f_scenecol,postminprestd)=results.forget.neg.bar.postminuspre.scene.std;
% % % % % %     % % % % %
% % % % % %     % % % % %
% % % % % %     %%%%%%%% Neutral %%%%%%%%%%%%%%%%%%%%%%%
% % % % % %     
% % % % % %     % Remember Precue
% % % % % %     allsub_RF_Neu(sub,r_wordcol,precuecol)=results.remember.neu.bar.precue.word.nanmean;
% % % % % %     allsub_RF_Neu(sub,r_scenecol,precuecol)=results.remember.neu.bar.precue.scene.nanmean;
% % % % % %     allsub_RF_Neu(sub,r_facecol,precuecol)=results.remember.neu.bar.precue.face.nanmean;
% % % % % %     allsub_RF_Neu(sub,r_objcol,precuecol)=results.remember.neu.bar.precue.object.nanmean;
% % % % % %     allsub_RF_Neu(sub,r_restcol,precuecol)=results.remember.neu.bar.precue.rest.nanmean;
% % % % % %     
% % % % % %     allsub_RF_Neu(sub,r_wordcol,precuestd)=results.remember.neu.bar.precue.word.std;
% % % % % %     allsub_RF_Neu(sub,r_scenecol,precuestd)=results.remember.neu.bar.precue.scene.std;
% % % % % %     allsub_RF_Neu(sub,r_facecol,precuestd)=results.remember.neu.bar.precue.face.std;
% % % % % %     allsub_RF_Neu(sub,r_objcol,precuestd)=results.remember.neu.bar.precue.object.std;
% % % % % %     allsub_RF_Neu(sub,r_restcol,precuestd)=results.remember.neu.bar.precue.rest.std;
% % % % % %     
% % % % % %     
% % % % % %     % Remember Postcue
% % % % % %     allsub_RF_Neu(sub,r_wordcol,postcuecol)=results.remember.neu.bar.postcue.word.nanmean;
% % % % % %     allsub_RF_Neu(sub,r_scenecol,postcuecol)=results.remember.neu.bar.postcue.scene.nanmean;
% % % % % %     allsub_RF_Neu(sub,r_facecol,postcuecol)=results.remember.neu.bar.postcue.face.nanmean;
% % % % % %     allsub_RF_Neu(sub,r_objcol,postcuecol)=results.remember.neu.bar.postcue.object.nanmean;
% % % % % %     allsub_RF_Neu(sub,r_restcol,postcuecol)=results.remember.neu.bar.postcue.rest.nanmean;
% % % % % %     
% % % % % %     allsub_RF_Neu(sub,r_wordcol,postcuestd)=results.remember.neu.bar.postcue.word.std;
% % % % % %     allsub_RF_Neu(sub,r_scenecol,postcuestd)=results.remember.neu.bar.postcue.scene.std;
% % % % % %     allsub_RF_Neu(sub,r_facecol,postcuestd)=results.remember.neu.bar.postcue.face.std;
% % % % % %     allsub_RF_Neu(sub,r_objcol,postcuestd)=results.remember.neu.bar.postcue.object.std;
% % % % % %     allsub_RF_Neu(sub,r_restcol,postcuestd)=results.remember.neu.bar.postcue.rest.std;
% % % % % %     
% % % % % %     
% % % % % %     % % % % %     % Remember Post-Pre
% % % % % %     % % % % %     allsub_RF_Neu(sub,r_wordcol,postminprecol)=results.remember.neu.bar.postminuspre.word.nanmean;
% % % % % %     % % % % %     allsub_RF_Neu(sub,r_scenecol,postminprecol)=results.remember.neu.bar.postminuspre.scene.nanmean;
% % % % % %     % % % % %
% % % % % %     % % % % %
% % % % % %     % % % % %     allsub_RF_Neu(sub,r_wordcol,postminprestd)=results.remember.neu.bar.postminuspre.word.std;
% % % % % %     % % % % %     allsub_RF_Neu(sub,r_scenecol,postminprestd)=results.remember.neu.bar.postminuspre.scene.std;
% % % % % %     % % % % %
% % % % % %     % % % % %
% % % % % %     
% % % % % %     % Forget Precue
% % % % % %     allsub_RF_Neu(sub,f_wordcol,precuecol)=results.forget.neu.bar.precue.word.nanmean;
% % % % % %     allsub_RF_Neu(sub,f_scenecol,precuecol)=results.forget.neu.bar.precue.scene.nanmean;
% % % % % %     allsub_RF_Neu(sub,f_facecol,precuecol)=results.forget.neu.bar.precue.face.nanmean;
% % % % % %     allsub_RF_Neu(sub,f_objcol,precuecol)=results.forget.neu.bar.precue.object.nanmean;
% % % % % %     allsub_RF_Neu(sub,f_restcol,precuecol)=results.forget.neu.bar.precue.rest.nanmean;
% % % % % %     
% % % % % %     allsub_RF_Neu(sub,f_wordcol,precuestd)=results.forget.neu.bar.precue.word.std;
% % % % % %     allsub_RF_Neu(sub,f_scenecol,precuestd)=results.forget.neu.bar.precue.scene.std;
% % % % % %     allsub_RF_Neu(sub,f_facecol,precuestd)=results.forget.neu.bar.precue.face.std;
% % % % % %     allsub_RF_Neu(sub,f_objcol,precuestd)=results.forget.neu.bar.precue.object.std;
% % % % % %     allsub_RF_Neu(sub,f_restcol,precuestd)=results.forget.neu.bar.precue.rest.std;
% % % % % %     
% % % % % %     % Forget Postcue
% % % % % %     allsub_RF_Neu(sub,f_wordcol,postcuecol)=results.forget.neu.bar.postcue.word.nanmean;
% % % % % %     allsub_RF_Neu(sub,f_scenecol,postcuecol)=results.forget.neu.bar.postcue.scene.nanmean;
% % % % % %     allsub_RF_Neu(sub,f_facecol,postcuecol)=results.forget.neu.bar.postcue.face.nanmean;
% % % % % %     allsub_RF_Neu(sub,f_objcol,postcuecol)=results.forget.neu.bar.postcue.object.nanmean;
% % % % % %     allsub_RF_Neu(sub,f_restcol,postcuecol)=results.forget.neu.bar.postcue.rest.nanmean;
% % % % % %     
% % % % % %     
% % % % % %     allsub_RF_Neu(sub,f_wordcol,postcuestd)=results.forget.neu.bar.postcue.word.std;
% % % % % %     allsub_RF_Neu(sub,f_scenecol,postcuestd)=results.forget.neu.bar.postcue.scene.std;
% % % % % %     allsub_RF_Neu(sub,f_facecol,postcuestd)=results.forget.neu.bar.postcue.face.std;
% % % % % %     allsub_RF_Neu(sub,f_objcol,postcuestd)=results.forget.neu.bar.postcue.object.std;
% % % % % %     allsub_RF_Neu(sub,f_restcol,postcuestd)=results.forget.neu.bar.postcue.rest.std;
% % % % % %     
% % % % % %     
% % % % % %     % % % % %
% % % % % %     % % % % %     % Forget Post-Pre
% % % % % %     % % % % %     allsub_RF_Neu(sub,f_wordcol,postminprecol)=results.forget.neu.bar.postminuspre.word.nanmean;
% % % % % %     % % % % %     allsub_RF_Neu(sub,f_scenecol,postminprecol)=results.forget.neu.bar.postminuspre.scene.nanmean;
% % % % % %     % % % % %
% % % % % %     % % % % %     allsub_RF_Neu(sub,f_wordcol,postminprestd)=results.forget.neu.bar.postminuspre.word.std;
% % % % % %     % % % % %     allsub_RF_Neu(sub,f_scenecol,postminprestd)=results.forget.neu.bar.postminuspre.scene.std;
% % % % % %     % % % % %
% % % % % %     % % % % %
% % % % % %     
% % % % % %     %%%%%%%%%  By TR timecourse graphs ( 1 by 7 timepoints) %%%%%%%%%
% % % % % %     
% % % % % %     %%%%%% Neg
% % % % % %     % Remember : Neg
% % % % % %     allsub_RF_Neg(sub,r_wordcol,TR1col:TR7col)=results.remember.neg.mean.word;
% % % % % %     allsub_RF_Neg(sub,r_scenecol,TR1col:TR7col)=results.remember.neg.mean.scene;
% % % % % %     allsub_RF_Neg(sub,r_facecol,TR1col:TR7col)=results.remember.neg.mean.face;
% % % % % %     allsub_RF_Neg(sub,r_objcol,TR1col:TR7col)=results.remember.neg.mean.object;
% % % % % %     allsub_RF_Neg(sub,r_restcol,TR1col:TR7col)=results.remember.neg.mean.rest;
% % % % % %     
% % % % % %     allsub_RF_Neg(sub,r_wordcol,TR1std:TR7std)=results.remember.neg.std.word;
% % % % % %     allsub_RF_Neg(sub,r_scenecol,TR1std:TR7std)=results.remember.neg.std.scene;
% % % % % %     allsub_RF_Neg(sub,r_facecol,TR1std:TR7std)=results.remember.neg.std.face;
% % % % % %     allsub_RF_Neg(sub,r_objcol,TR1std:TR7std)=results.remember.neg.std.object;
% % % % % %     allsub_RF_Neg(sub,r_restcol,TR1std:TR7std)=results.remember.neg.std.rest;
% % % % % %     
% % % % % %     
% % % % % %     
% % % % % %     
% % % % % %     % Forget Neg
% % % % % %     
% % % % % %     allsub_RF_Neg(sub,f_wordcol,TR1col:TR7col)=results.forget.neg.mean.word;
% % % % % %     allsub_RF_Neg(sub,f_scenecol,TR1col:TR7col)=results.forget.neg.mean.scene;
% % % % % %     allsub_RF_Neg(sub,f_facecol,TR1col:TR7col)=results.forget.neg.mean.face;
% % % % % %     allsub_RF_Neg(sub,f_objcol,TR1col:TR7col)=results.forget.neg.mean.object;
% % % % % %     allsub_RF_Neg(sub,f_restcol,TR1col:TR7col)=results.forget.neg.mean.rest;
% % % % % %     
% % % % % %     allsub_RF_Neg(sub,f_wordcol,TR1std:TR7std)=results.forget.neg.std.word;
% % % % % %     allsub_RF_Neg(sub,f_scenecol,TR1std:TR7std)=results.forget.neg.std.scene;
% % % % % %     allsub_RF_Neg(sub,f_facecol,TR1std:TR7std)=results.forget.neg.std.face;
% % % % % %     allsub_RF_Neg(sub,f_objcol,TR1std:TR7std)=results.forget.neg.std.object;
% % % % % %     allsub_RF_Neg(sub,f_restcol,TR1std:TR7std)=results.forget.neg.std.rest;
% % % % % %     
% % % % % %     
% % % % % %     
% % % % % %     %%%%%% Neutral
% % % % % %     % Remember Neutral
% % % % % %     allsub_RF_Neu(sub,r_wordcol,TR1col:TR7col)=results.remember.neu.mean.word;
% % % % % %     allsub_RF_Neu(sub,r_scenecol,TR1col:TR7col)=results.remember.neu.mean.scene;
% % % % % %     allsub_RF_Neu(sub,r_facecol,TR1col:TR7col)=results.remember.neu.mean.face;
% % % % % %     allsub_RF_Neu(sub,r_objcol,TR1col:TR7col)=results.remember.neu.mean.object;
% % % % % %     allsub_RF_Neu(sub,r_restcol,TR1col:TR7col)=results.remember.neu.mean.rest;
% % % % % %     
% % % % % %     allsub_RF_Neu(sub,r_wordcol,TR1std:TR7std)=results.remember.neu.std.word;
% % % % % %     allsub_RF_Neu(sub,r_scenecol,TR1std:TR7std)=results.remember.neu.std.scene;
% % % % % %     allsub_RF_Neu(sub,r_facecol,TR1std:TR7std)=results.remember.neu.std.face;
% % % % % %     allsub_RF_Neu(sub,r_objcol,TR1std:TR7std)=results.remember.neu.std.object;
% % % % % %     allsub_RF_Neu(sub,r_restcol,TR1std:TR7std)=results.remember.neu.std.rest;
% % % % % %     
% % % % % %     
% % % % % %     
% % % % % %     
% % % % % %     % forget Neutral
% % % % % %     
% % % % % %     allsub_RF_Neu(sub,f_wordcol,TR1col:TR7col)=results.forget.neu.mean.word;
% % % % % %     allsub_RF_Neu(sub,f_scenecol,TR1col:TR7col)=results.forget.neu.mean.scene;
% % % % % %     allsub_RF_Neu(sub,f_facecol,TR1col:TR7col)=results.forget.neu.mean.face;
% % % % % %     allsub_RF_Neu(sub,f_objcol,TR1col:TR7col)=results.forget.neu.mean.object;
% % % % % %     allsub_RF_Neu(sub,f_restcol,TR1col:TR7col)=results.forget.neu.mean.rest;
% % % % % %     
% % % % % %     allsub_RF_Neu(sub,f_wordcol,TR1std:TR7std)=results.forget.neu.std.word;
% % % % % %     allsub_RF_Neu(sub,f_scenecol,TR1std:TR7std)=results.forget.neu.std.scene;
% % % % % %     allsub_RF_Neu(sub,f_facecol,TR1std:TR7std)=results.forget.neu.std.face;
% % % % % %     allsub_RF_Neu(sub,f_objcol,TR1std:TR7std)=results.forget.neu.std.object;
% % % % % %     allsub_RF_Neu(sub,f_restcol,TR1std:TR7std)=results.forget.neu.std.rest;
% % % % % %     
% % % % % %     
% % % % % % end
% % % % % % 
% % % % % % 
% % % % % % % Flatten matrix to 2D to print to excel\
% % % % % % allsub_NegRF=squeeze(allsub_RF_Neg(1,:,:)); %squeeze removes singleton 1st dimension of sub.
% % % % % % allsub_NeuRF=squeeze(allsub_RF_Neu(1,:,:));
% % % % % % 
% % % % % % for sub=2:size(sublist,2)
% % % % % %     allsub_NegRF=vertcat(allsub_NegRF,squeeze(allsub_RF_Neg(sub,:,:)));
% % % % % %     allsub_NeuRF=vertcat(allsub_NeuRF,squeeze(allsub_RF_Neu(sub,:,:)));
% % % % % % end
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % % output to excel
% % % % % % cd(outdir);
% % % % % % xlsheader={'subject','Category_fsowr','Rem/For','Precue','Postcue','PostminPre','TR1','TR2','TR3','TR4','TR5','TR6','TR7','PrecueSTD','PostcueSTD','PostminPreSTD', 'TR1STD','TR2STD','TR3STD','TR4STD','TR5STD','TR6STd','TR7STD'};
% % % % % % xlswrite(outname,xlsheader,'Neg','A1');
% % % % % % xlswrite(outname,allsub_NegRF,'Neg','A2');
% % % % % % xlswrite(outname,xlsheader,'Neut','A1');
% % % % % % xlswrite(outname,allsub_NeuRF,'Neut','A2');
% % % % % % 
% % % % % % % save data in .MAT file
% % % % % % save(matname,'allsub_NegRF','allsub_NeuRF','-append');

