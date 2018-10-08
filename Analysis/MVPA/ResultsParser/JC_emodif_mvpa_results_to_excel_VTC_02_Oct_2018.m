% Wrapper file for grabbing classification results from Tracy's Parsed
% files  "phase_bar_compare.mat"
%
% First used for PHC Masks
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
% Sept. 2018
% Oct. 2018 -- Modified to include fields for Emotion



%% Before running whole script
clear all; close all; clc;

computeruse=1;  %  1 for mac, 2 for PC

%%  Input Arguments %%%%%%%



% code directory
codedirMAC='/Users/jc/Box Sync/CMF_BoxWork/Exp1_JudyMVPA/EmoDF_Scripts/EmoDF_LocGit/MVPA_Results';
codedirPC=' ' ;

% slash
slashMAC='/';
slashPC='\';

% Data directory
datadirMAC='/Users/jc/Box Sync/SharedFiles/EmoDF_LewPeaLab/emodif_data/';
datadirPC='';

% Output directory
outdirMAC='/Users/jc/EmoDFmac/Analysis/MVPA_Results';
outdirPC='';

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

outname=strcat('emodif_',phasedirstr,'_',maskdirstr,'_classEvidence_',date,'.xlsx');

% output .MAT file name
matname=strcat('EmoDiF_mvpa_',maskdirstr,'_results_to_excel_',date,'.mat');


%% Parameters for script

% row and column indices
% stim type - -allsub_data(:,%,:)
r_wordcol=1;
r_scenecol=2;
r_facecol=3;
r_objcol=4;
r_restcol=5;
f_wordcol=6;
f_scenecol=7;
f_facecol=8;
f_objcol=9;
f_restcol=10;


%Data indices allsub_data(:,:,%)
subcol=1;
precuecol=4;
postcuecol=5;
postminprecol=6;
TR1col=7;
TR2col=8;
TR3col=9;
TR4col=10;
TR5col=11;
TR6col=12;
TR7col=13;

precuestd=14;
postcuestd=15;
postminprestd=16;
TR1std=17;
TR7std=23;

% Grab data from Results folder for DELL
clear allsub_RF;
allsub_RF=nan(size(sublist,2),10,23);

for sub=1:size(sublist,2)
    clear results;
    
    subdir=strcat('emodif_', num2str(sublist(sub)),slashsym);
    datedir=datecell{datedir_ind(sub)};
    
    % load results struct from the specified phase
    % specify date by subject
    loaddir=strcat(datadir,subdir,'results',slashsym,phasedir,maskdir,datedir)
    
    load(strcat(loaddir,'results_DFencode_bar_compare.mat'));
    
    % Participant and condition info
    allsub_RF(sub,:,1)=sublist(sub);
    allsub_RF(sub,1:5,2)=[4,2,1,3,5]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
    allsub_RF(sub,6:10,2)=[4,2,1,3,5];
    
    allsub_RF(sub,1:5,3)=1;  %remember
    allsub_RF(sub,6:10,3)=2; % forget
    
    
    
    %% By main effect of Remember vs. Forget %%%%%
    
    %%%%%%%%%%  Averaged TR bar graphs ( 1 number per data point) %%%%%%%%%%%%
    
    % Remember Precue
    allsub_RF(sub,r_wordcol,precuecol)=results.remember.all.bar.precue.word.nanmean;
    allsub_RF(sub,r_scenecol,precuecol)=results.remember.all.bar.precue.scene.nanmean;
    allsub_RF(sub,r_facecol,precuecol)=results.remember.all.bar.precue.face.nanmean;
    allsub_RF(sub,r_objcol,precuecol)= results.remember.all.bar.precue.object.nanmean;
    allsub_RF(sub,r_restcol,precuecol)=results.remember.all.bar.precue.rest.nanmean;
    
    allsub_RF(sub,r_wordcol,precuestd)=results.remember.all.bar.precue.word.std;
    allsub_RF(sub,r_scenecol,precuestd)=results.remember.all.bar.precue.scene.std;
    allsub_RF(sub,r_facecol,precuestd)=results.remember.all.bar.precue.face.std;
    allsub_RF(sub,r_objcol,precuestd)=results.remember.all.bar.precue.object.std;
    allsub_RF(sub,r_restcol,precuestd)=results.remember.all.bar.precue.rest.std;
    
    % Remember Postcue
    
    allsub_RF(sub,r_wordcol,postcuecol)=results.remember.all.bar.postcue.word.nanmean;
    allsub_RF(sub,r_scenecol,postcuecol)=results.remember.all.bar.postcue.scene.nanmean;
    allsub_RF(sub,r_facecol,postcuecol)=results.remember.all.bar.postcue.face.nanmean;
    allsub_RF(sub,r_objcol,postcuecol)=results.remember.all.bar.postcue.object.nanmean;
    allsub_RF(sub,r_restcol,postcuecol)=results.remember.all.bar.postcue.rest.nanmean;
    
    allsub_RF(sub,r_wordcol,postcuestd)=results.remember.all.bar.postcue.word.std;
    allsub_RF(sub,r_scenecol,postcuestd)=results.remember.all.bar.postcue.scene.std;
    allsub_RF(sub,r_facecol,postcuestd)=results.remember.all.bar.postcue.face.std;
    allsub_RF(sub,r_objcol,postcuestd)=results.remember.all.bar.postcue.object.std;
    allsub_RF(sub,r_restcol,postcuestd)=results.remember.all.bar.postcue.rest.std;
    
    % Remember Post-Pre
    allsub_RF(sub,r_wordcol,postminprecol)=results.remember.all.bar.postminuspre.nofacenorm.word.nanmean;
    allsub_RF(sub,r_scenecol,postminprecol)=results.remember.all.bar.postminuspre.nofacenorm.scene.nanmean;
    
    
    allsub_RF(sub,r_wordcol,postminprestd)=results.remember.all.bar.postminuspre.nofacenorm.word.std;
    allsub_RF(sub,r_scenecol,postminprestd)=results.remember.all.bar.postminuspre.nofacenorm.scene.std;
    
    % Forget Precue
    
    allsub_RF(sub,f_wordcol,precuecol)=results.forget.all.bar.precue.word.nanmean;
    allsub_RF(sub,f_scenecol,precuecol)=results.forget.all.bar.precue.scene.nanmean;
    allsub_RF(sub,f_facecol,precuecol)=results.forget.all.bar.precue.face.nanmean;
    allsub_RF(sub,f_objcol,precuecol)=results.forget.all.bar.precue.object.nanmean;
    allsub_RF(sub,f_restcol,precuecol)=results.forget.all.bar.precue.rest.nanmean;
    
    allsub_RF(sub,f_wordcol,precuestd)=results.forget.all.bar.precue.word.std;
    allsub_RF(sub,f_scenecol,precuestd)=results.forget.all.bar.precue.scene.std;
    allsub_RF(sub,f_facecol,precuestd)=results.forget.all.bar.precue.face.std;
    allsub_RF(sub,f_objcol,precuestd)=results.forget.all.bar.precue.object.std;
    allsub_RF(sub,f_restcol,precuestd)=results.forget.all.bar.precue.rest.std;
    
    
    % Forget Postcue
    allsub_RF(sub,f_wordcol,postcuecol)=results.forget.all.bar.postcue.word.nanmean;
    allsub_RF(sub,f_scenecol,postcuecol)=results.forget.all.bar.postcue.scene.nanmean;
    allsub_RF(sub,f_facecol,postcuecol)=results.forget.all.bar.postcue.face.nanmean;
    allsub_RF(sub,f_objcol,postcuecol)=results.forget.all.bar.postcue.object.nanmean;
    allsub_RF(sub,f_restcol,postcuecol)=results.forget.all.bar.postcue.rest.nanmean;
    
    allsub_RF(sub,f_wordcol,postcuestd)=results.forget.all.bar.postcue.word.std;
    allsub_RF(sub,f_scenecol,postcuestd)=results.forget.all.bar.postcue.scene.std;
    allsub_RF(sub,f_facecol,postcuestd)=results.forget.all.bar.postcue.face.std;
    allsub_RF(sub,f_objcol,postcuestd)=results.forget.all.bar.postcue.object.std;
    allsub_RF(sub,f_restcol,postcuestd)=results.forget.all.bar.postcue.rest.std;
    
    
    % Forget Post-Pre
    allsub_RF(sub,f_wordcol,postminprecol)=results.forget.all.bar.postminuspre.nofacenorm.word.nanmean;
    allsub_RF(sub,f_scenecol,postminprecol)=results.forget.all.bar.postminuspre.nofacenorm.scene.nanmean;
    
    allsub_RF(sub,f_wordcol,postminprestd)=results.forget.all.bar.postminuspre.nofacenorm.word.std;
    allsub_RF(sub,f_scenecol,postminprestd)=results.forget.all.bar.postminuspre.nofacenorm.scene.std;
    
    
    %%%%%%%%%  By TR timecourse graphs ( 1 by 7 timepoints) %%%%%%%%%
    % Remember
    allsub_RF(sub,r_wordcol,TR1col:TR7col)=results.remember.all.mean.word;
    allsub_RF(sub,r_scenecol,TR1col:TR7col)=results.remember.all.mean.scene;
    allsub_RF(sub,r_facecol,TR1col:TR7col)=results.remember.all.mean.face;
    allsub_RF(sub,r_objcol,TR1col:TR7col)=results.remember.all.mean.object;
    allsub_RF(sub,r_restcol,TR1col:TR7col)=results.remember.all.mean.rest;
    
    allsub_RF(sub,r_wordcol,TR1std:TR7std)=results.remember.all.std.word;
    allsub_RF(sub,r_scenecol,TR1std:TR7std)=results.remember.all.std.scene;
    allsub_RF(sub,r_facecol,TR1std:TR7std)=results.remember.all.std.face;
    allsub_RF(sub,r_objcol,TR1std:TR7std)=results.remember.all.std.object;
    allsub_RF(sub,r_restcol,TR1std:TR7std)=results.remember.all.std.rest;
    
    % Forget
    allsub_RF(sub,f_wordcol,TR1col:TR7col)=results.forget.all.mean.word;
    allsub_RF(sub,f_scenecol,TR1col:TR7col)=results.forget.all.mean.scene;
    allsub_RF(sub,f_facecol,TR1col:TR7col)=results.forget.all.mean.face;
    allsub_RF(sub,f_objcol,TR1col:TR7col)=results.forget.all.mean.object;
    allsub_RF(sub,f_restcol,TR1col:TR7col)=results.forget.all.mean.rest;
    
    allsub_RF(sub,f_wordcol,TR1std:TR7std)=results.forget.all.std.word;
    allsub_RF(sub,f_scenecol,TR1std:TR7std)=results.forget.all.std.scene;
    allsub_RF(sub,f_facecol,TR1std:TR7std)=results.forget.all.std.face;
    allsub_RF(sub,f_objcol,TR1std:TR7std)=results.forget.all.std.object;
    allsub_RF(sub,f_restcol,TR1std:TR7std)=results.forget.all.std.rest;
end


% Flatten matrix to 2D to print to excel\
allsub_xls=squeeze(allsub_RF(1,:,:)); %squeeze removes singleton 1st dimension of sub.
for sub=2:size(sublist,2)
    allsub_xls=vertcat(allsub_xls,squeeze(allsub_RF(sub,:,:)));
end

% output to excel
cd(outdir);

% % % outname=strcat('emodif_',phasedirstr,'_',maskdirstr,'_classEvidence_',date,'.xlsx');
% % % xlsheader={'subject','Category_fsowr','Rem/For','Precue','Postcue','PostminPre','TR1','TR2','TR3','TR4','TR5','TR6','TR7','PrecueSTD','PostcueSTD','PostminPreSTD', 'TR1STD','TR2STD','TR3STD','TR4STD','TR5STD','TR6STd','TR7STD'};
% % % xlswrite(outname,xlsheader,'mainRF','A1');
% % % xlswrite(outname,allsub_xls,'mainRF','A2');

% save data in .MAT file
save(matname,'allsub_xls','datadir','codedir','outdir','sublist','datecell','datedir_ind');
cd(codedir);

%%  By Interation of DF * Accuracy  %%%%%%%%%%%%%%%%%%%

clear allsub_RF_corr;
clear allsub_RF_incorr;
allsub_RF_corr=nan(size(sublist,2),10,23);
allsub_RF_incorr=nan(size(sublist,2),10,23);

for sub=1:size(sublist,2)
    clear results;
    
    subdir=strcat('emodif_', num2str(sublist(sub)),slashsym);
    datedir=datecell{datedir_ind(sub)};
    
    % load results struct from the specified phase
    % specify date by subject
    loaddir=strcat(datadir,subdir,'results',slashsym,phasedir,maskdir,datedir)
    
    load(strcat(loaddir,'results_DFencode_bar_compare.mat'));
    
    % Participant and condition info
    allsub_RF_corr(sub,:,1)=sublist(sub);
    allsub_RF_corr(sub,1:5,2)=[4,2,1,3,5]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
    allsub_RF_corr(sub,6:10,2)=[4,2,1,3,5];
    
    allsub_RF_incorr(sub,:,1)=sublist(sub);
    allsub_RF_incorr(sub,1:5,2)=[4,2,1,3,5]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
    allsub_RF_incorr(sub,6:10,2)=[4,2,1,3,5];
    
    allsub_RF_corr(sub,1:5,3)=1;  %remember
    allsub_RF_corr(sub,6:10,3)=2; % forget
    
    allsub_RF_incorr(sub,1:5,3)=1;  %remember
    allsub_RF_incorr(sub,6:10,3)=2; % forget
    
    
    %%%%%%%%%%  Averaged TR bar graphs ( 1 number per data point) %%%%%%%%%%%%
    
    %%%%%%% correct %%%%%%%%%%%%%%%5
    % Remember Precue
    allsub_RF_corr(sub,r_wordcol,precuecol)=results.remember.merged.bar.OLD.precue.word.nanmean;
    allsub_RF_corr(sub,r_scenecol,precuecol)=results.remember.merged.bar.OLD.precue.scene.nanmean;
    allsub_RF_corr(sub,r_facecol,precuecol)=results.remember.merged.bar.OLD.precue.face.nanmean;
    allsub_RF_corr(sub,r_objcol,precuecol)=results.remember.merged.bar.OLD.precue.object.nanmean;
    allsub_RF_corr(sub,r_restcol,precuecol)=results.remember.merged.bar.OLD.precue.rest.nanmean;
    
    allsub_RF_corr(sub,r_wordcol,precuestd)=results.remember.merged.bar.OLD.precue.word.std;
    allsub_RF_corr(sub,r_scenecol,precuestd)=results.remember.merged.bar.OLD.precue.scene.std;
    allsub_RF_corr(sub,r_facecol,precuestd)=results.remember.merged.bar.OLD.precue.face.std;
    allsub_RF_corr(sub,r_objcol,precuestd)=results.remember.merged.bar.OLD.precue.object.std;
    allsub_RF_corr(sub,r_restcol,precuestd)=results.remember.merged.bar.OLD.precue.rest.std;
    
    
    % Remember Postcue
    allsub_RF_corr(sub,r_wordcol,postcuecol)=results.remember.merged.bar.OLD.postcue.word.nanmean;
    allsub_RF_corr(sub,r_scenecol,postcuecol)=results.remember.merged.bar.OLD.postcue.scene.nanmean;
    allsub_RF_corr(sub,r_facecol,postcuecol)=results.remember.merged.bar.OLD.postcue.face.nanmean;
    allsub_RF_corr(sub,r_objcol,postcuecol)=results.remember.merged.bar.OLD.postcue.object.nanmean;
    allsub_RF_corr(sub,r_restcol,postcuecol)=results.remember.merged.bar.OLD.postcue.rest.nanmean;
    
    allsub_RF_corr(sub,r_wordcol,postcuestd)=results.remember.merged.bar.OLD.postcue.word.std;
    allsub_RF_corr(sub,r_scenecol,postcuestd)=results.remember.merged.bar.OLD.postcue.scene.std;
    allsub_RF_corr(sub,r_facecol,postcuestd)=results.remember.merged.bar.OLD.postcue.face.std;
    allsub_RF_corr(sub,r_objcol,postcuestd)=results.remember.merged.bar.OLD.postcue.object.std;
    allsub_RF_corr(sub,r_restcol,postcuestd)=results.remember.merged.bar.OLD.postcue.rest.std;
    
    
    % % %         % Remember Post-Pre
    % % %         allsub_RF_corr(sub,r_wordcol,postminprecol)=results.remember.merged.bar.OLD.postminuspre.word.nanmean;
    % % %         allsub_RF_corr(sub,r_scenecol,postminprecol)=results.remember.merged.bar.OLD.postminuspre.scene.nanmean;
    % % %
    % % %
    % % %         allsub_RF_corr(sub,r_wordcol,postminprestd)=results.remember.merged.bar.OLD.postminuspre.word.std;
    % % %         allsub_RF_corr(sub,r_scenecol,postminprestd)=results.remember.merged.bar.OLD.postminuspre.scene.std;
    % % %
    % % %
    
    % Forget Precue
    allsub_RF_corr(sub,f_wordcol,precuecol)=results.forget.merged.bar.OLD.precue.word.nanmean;
    allsub_RF_corr(sub,f_scenecol,precuecol)=results.forget.merged.bar.OLD.precue.scene.nanmean;
    allsub_RF_corr(sub,f_facecol,precuecol)=results.forget.merged.bar.OLD.precue.face.nanmean;
    allsub_RF_corr(sub,f_objcol,precuecol)=results.forget.merged.bar.OLD.precue.object.nanmean;
    allsub_RF_corr(sub,f_restcol,precuecol)=results.forget.merged.bar.OLD.precue.rest.nanmean;
    
    allsub_RF_corr(sub,f_wordcol,precuestd)=results.forget.merged.bar.OLD.precue.word.std;
    allsub_RF_corr(sub,f_scenecol,precuestd)=results.forget.merged.bar.OLD.precue.scene.std;
    allsub_RF_corr(sub,f_facecol,precuestd)=results.forget.merged.bar.OLD.precue.face.std;
    allsub_RF_corr(sub,f_objcol,precuestd)=results.forget.merged.bar.OLD.precue.object.std;
    allsub_RF_corr(sub,f_restcol,precuestd)=results.forget.merged.bar.OLD.precue.rest.std;
    
    % Forget Postcue
    allsub_RF_corr(sub,f_wordcol,postcuecol)=results.forget.merged.bar.OLD.postcue.word.nanmean;
    allsub_RF_corr(sub,f_scenecol,postcuecol)=results.forget.merged.bar.OLD.postcue.scene.nanmean;
    allsub_RF_corr(sub,f_facecol,postcuecol)=results.forget.merged.bar.OLD.postcue.face.nanmean;
    allsub_RF_corr(sub,f_objcol,postcuecol)=results.forget.merged.bar.OLD.postcue.object.nanmean;
    allsub_RF_corr(sub,f_restcol,postcuecol)=results.forget.merged.bar.OLD.postcue.rest.nanmean;
    
    
    allsub_RF_corr(sub,f_wordcol,postcuestd)=results.forget.merged.bar.OLD.postcue.word.std;
    allsub_RF_corr(sub,f_scenecol,postcuestd)=results.forget.merged.bar.OLD.postcue.scene.std;
    allsub_RF_corr(sub,f_facecol,postcuestd)=results.forget.merged.bar.OLD.postcue.face.std;
    allsub_RF_corr(sub,f_objcol,postcuestd)=results.forget.merged.bar.OLD.postcue.object.std;
    allsub_RF_corr(sub,f_restcol,postcuestd)=results.forget.merged.bar.OLD.postcue.rest.std;
    
    
    % % %
    % % %         % Forget Post-Pre
    % % %         allsub_RF_corr(sub,f_wordcol,postminprecol)=results.forget.merged.bar.OLD.postminuspre.word.nanmean;
    % % %         allsub_RF_corr(sub,f_scenecol,postminprecol)=results.forget.merged.bar.OLD.postminuspre.scene.nanmean;
    % % %
    % % %         allsub_RF_corr(sub,f_wordcol,postminprestd)=results.forget.merged.bar.OLD.postminuspre.word.std;
    % % %         allsub_RF_corr(sub,f_scenecol,postminprestd)=results.forget.merged.bar.OLD.postminuspre.scene.std;
    % % %
    
    %%%%%%%% Incorrect %%%%%%%%%%%%%%%%%%%%%%%5
    
    % Remember Precue
    allsub_RF_incorr(sub,r_wordcol,precuecol)=results.remember.merged.bar.NEW.precue.word.nanmean;
    allsub_RF_incorr(sub,r_scenecol,precuecol)=results.remember.merged.bar.NEW.precue.scene.nanmean;
    allsub_RF_incorr(sub,r_facecol,precuecol)=results.remember.merged.bar.NEW.precue.face.nanmean;
    allsub_RF_incorr(sub,r_objcol,precuecol)=results.remember.merged.bar.NEW.precue.object.nanmean;
    allsub_RF_incorr(sub,r_restcol,precuecol)=results.remember.merged.bar.NEW.precue.rest.nanmean;
    
    allsub_RF_incorr(sub,r_wordcol,precuestd)=results.remember.merged.bar.NEW.precue.word.std;
    allsub_RF_incorr(sub,r_scenecol,precuestd)=results.remember.merged.bar.NEW.precue.scene.std;
    allsub_RF_incorr(sub,r_facecol,precuestd)=results.remember.merged.bar.NEW.precue.face.std;
    allsub_RF_incorr(sub,r_objcol,precuestd)=results.remember.merged.bar.NEW.precue.object.std;
    allsub_RF_incorr(sub,r_restcol,precuestd)=results.remember.merged.bar.NEW.precue.rest.std;
    
    
    % Remember Postcue
    allsub_RF_incorr(sub,r_wordcol,postcuecol)=results.remember.merged.bar.NEW.postcue.word.nanmean;
    allsub_RF_incorr(sub,r_scenecol,postcuecol)=results.remember.merged.bar.NEW.postcue.scene.nanmean;
    allsub_RF_incorr(sub,r_facecol,postcuecol)=results.remember.merged.bar.NEW.postcue.face.nanmean;
    allsub_RF_incorr(sub,r_objcol,postcuecol)=results.remember.merged.bar.NEW.postcue.object.nanmean;
    allsub_RF_incorr(sub,r_restcol,postcuecol)=results.remember.merged.bar.NEW.postcue.rest.nanmean;
    
    allsub_RF_incorr(sub,r_wordcol,postcuestd)=results.remember.merged.bar.NEW.postcue.word.std;
    allsub_RF_incorr(sub,r_scenecol,postcuestd)=results.remember.merged.bar.NEW.postcue.scene.std;
    allsub_RF_incorr(sub,r_facecol,postcuestd)=results.remember.merged.bar.NEW.postcue.face.std;
    allsub_RF_incorr(sub,r_objcol,postcuestd)=results.remember.merged.bar.NEW.postcue.object.std;
    allsub_RF_incorr(sub,r_restcol,postcuestd)=results.remember.merged.bar.NEW.postcue.rest.std;
    
    % % % %
    % % % %         % Remember Post-Pre
    % % % %         allsub_RF_incorr(sub,r_wordcol,postminprecol)=results.remember.merged.bar.NEW.postminuspre.word.nanmean;
    % % % %         allsub_RF_incorr(sub,r_scenecol,postminprecol)=results.remember.merged.bar.NEW.postminuspre.scene.nanmean;
    % % % %
    % % % %
    % % % %         allsub_RF_incorr(sub,r_wordcol,postminprestd)=results.remember.merged.bar.NEW.postminuspre.word.std;
    % % % %         allsub_RF_incorr(sub,r_scenecol,postminprestd)=results.remember.merged.bar.NEW.postminuspre.scene.std;
    
    
    
    % Forget Precue
    allsub_RF_incorr(sub,f_wordcol,precuecol)=results.forget.merged.bar.NEW.precue.word.nanmean;
    allsub_RF_incorr(sub,f_scenecol,precuecol)=results.forget.merged.bar.NEW.precue.scene.nanmean;
    allsub_RF_incorr(sub,f_facecol,precuecol)=results.forget.merged.bar.NEW.precue.face.nanmean;
    allsub_RF_incorr(sub,f_objcol,precuecol)=results.forget.merged.bar.NEW.precue.object.nanmean;
    allsub_RF_incorr(sub,f_restcol,precuecol)=results.forget.merged.bar.NEW.precue.rest.nanmean;
    
    allsub_RF_incorr(sub,f_wordcol,precuestd)=results.forget.merged.bar.NEW.precue.word.std;
    allsub_RF_incorr(sub,f_scenecol,precuestd)=results.forget.merged.bar.NEW.precue.scene.std;
    allsub_RF_incorr(sub,f_facecol,precuestd)=results.forget.merged.bar.NEW.precue.face.std;
    allsub_RF_incorr(sub,f_objcol,precuestd)=results.forget.merged.bar.NEW.precue.object.std;
    allsub_RF_incorr(sub,f_restcol,precuestd)=results.forget.merged.bar.NEW.precue.rest.std;
    
    % Forget Postcue
    allsub_RF_incorr(sub,f_wordcol,postcuecol)=results.forget.merged.bar.NEW.postcue.word.nanmean;
    allsub_RF_incorr(sub,f_scenecol,postcuecol)=results.forget.merged.bar.NEW.postcue.scene.nanmean;
    allsub_RF_incorr(sub,f_facecol,postcuecol)=results.forget.merged.bar.NEW.postcue.face.nanmean;
    allsub_RF_incorr(sub,f_objcol,postcuecol)=results.forget.merged.bar.NEW.postcue.object.nanmean;
    allsub_RF_incorr(sub,f_restcol,postcuecol)=results.forget.merged.bar.NEW.postcue.rest.nanmean;
    
    
    allsub_RF_incorr(sub,f_wordcol,postcuestd)=results.forget.merged.bar.NEW.postcue.word.std;
    allsub_RF_incorr(sub,f_scenecol,postcuestd)=results.forget.merged.bar.NEW.postcue.scene.std;
    allsub_RF_incorr(sub,f_facecol,postcuestd)=results.forget.merged.bar.NEW.postcue.face.std;
    allsub_RF_incorr(sub,f_objcol,postcuestd)=results.forget.merged.bar.NEW.postcue.object.std;
    allsub_RF_incorr(sub,f_restcol,postcuestd)=results.forget.merged.bar.NEW.postcue.rest.std;
    
    
    
    % % % %         % Forget Post-Pre
    % % % %         allsub_RF_incorr(sub,f_wordcol,postminprecol)=results.forget.merged.bar.NEW.postminuspre.word.nanmean;
    % % % %         allsub_RF_incorr(sub,f_scenecol,postminprecol)=results.forget.merged.bar.NEW.postminuspre.scene.nanmean;
    % % % %
    % % % %         allsub_RF_incorr(sub,f_wordcol,postminprestd)=results.forget.merged.bar.NEW.postminuspre.word.std;
    % % % %         allsub_RF_incorr(sub,f_scenecol,postminprestd)=results.forget.merged.bar.NEW.postminuspre.scene.std;
    % % % %
    % % % %
    
    %%%%%%%%%  By TR timecourse graphs ( 1 by 7 timepoints) %%%%%%%%%
    
    %%%%%% OLD (correct)
    % Remember : OLD (corr)
    allsub_RF_corr(sub,r_wordcol,TR1col:TR7col)=results.remember.merged.mean.OLD.word;
    allsub_RF_corr(sub,r_scenecol,TR1col:TR7col)=results.remember.merged.mean.OLD.scene;
    allsub_RF_corr(sub,r_facecol,TR1col:TR7col)=results.remember.merged.mean.OLD.face;
    allsub_RF_corr(sub,r_objcol,TR1col:TR7col)=results.remember.merged.mean.OLD.object;
    allsub_RF_corr(sub,r_restcol,TR1col:TR7col)=results.remember.merged.mean.OLD.rest;
    
    allsub_RF_corr(sub,r_wordcol,TR1std:TR7std)=results.remember.merged.std.OLD.word;
    allsub_RF_corr(sub,r_scenecol,TR1std:TR7std)=results.remember.merged.std.OLD.scene;
    allsub_RF_corr(sub,r_facecol,TR1std:TR7std)=results.remember.merged.std.OLD.face;
    allsub_RF_corr(sub,r_objcol,TR1std:TR7std)=results.remember.merged.std.OLD.object;
    allsub_RF_corr(sub,r_restcol,TR1std:TR7std)=results.remember.merged.std.OLD.rest;
    
    
    
    
    % Forget OLD (corr)
    
    allsub_RF_corr(sub,f_wordcol,TR1col:TR7col)=results.forget.merged.mean.OLD.word;
    allsub_RF_corr(sub,f_scenecol,TR1col:TR7col)=results.forget.merged.mean.OLD.scene;
    allsub_RF_corr(sub,f_facecol,TR1col:TR7col)=results.forget.merged.mean.OLD.face;
    allsub_RF_corr(sub,f_objcol,TR1col:TR7col)=results.forget.merged.mean.OLD.object;
    allsub_RF_corr(sub,f_restcol,TR1col:TR7col)=results.forget.merged.mean.OLD.rest;
    
    allsub_RF_corr(sub,f_wordcol,TR1std:TR7std)=results.forget.merged.std.OLD.word;
    allsub_RF_corr(sub,f_scenecol,TR1std:TR7std)=results.forget.merged.std.OLD.scene;
    allsub_RF_corr(sub,f_facecol,TR1std:TR7std)=results.forget.merged.std.OLD.face;
    allsub_RF_corr(sub,f_objcol,TR1std:TR7std)=results.forget.merged.std.OLD.object;
    allsub_RF_corr(sub,f_restcol,TR1std:TR7std)=results.forget.merged.std.OLD.rest;
    
    
    
    %%%%%% NEW (incorr)
    % Remember NEW (incorr)
    allsub_RF_incorr(sub,r_wordcol,TR1col:TR7col)=results.remember.merged.mean.NEW.word;
    allsub_RF_incorr(sub,r_scenecol,TR1col:TR7col)=results.remember.merged.mean.NEW.scene;
    allsub_RF_incorr(sub,r_facecol,TR1col:TR7col)=results.remember.merged.mean.NEW.face;
    allsub_RF_incorr(sub,r_objcol,TR1col:TR7col)=results.remember.merged.mean.NEW.object;
    allsub_RF_incorr(sub,r_restcol,TR1col:TR7col)=results.remember.merged.mean.NEW.rest;
    
    allsub_RF_incorr(sub,r_wordcol,TR1std:TR7std)=results.remember.merged.std.NEW.word;
    allsub_RF_incorr(sub,r_scenecol,TR1std:TR7std)=results.remember.merged.std.NEW.scene;
    allsub_RF_incorr(sub,r_facecol,TR1std:TR7std)=results.remember.merged.std.NEW.face;
    allsub_RF_incorr(sub,r_objcol,TR1std:TR7std)=results.remember.merged.std.NEW.object;
    allsub_RF_incorr(sub,r_restcol,TR1std:TR7std)=results.remember.merged.std.NEW.rest;
    
    
    
    
    % forget NEW (incorr)
    
    allsub_RF_incorr(sub,f_wordcol,TR1col:TR7col)=results.forget.merged.mean.NEW.word;
    allsub_RF_incorr(sub,f_scenecol,TR1col:TR7col)=results.forget.merged.mean.NEW.scene;
    allsub_RF_incorr(sub,f_facecol,TR1col:TR7col)=results.forget.merged.mean.NEW.face;
    allsub_RF_incorr(sub,f_objcol,TR1col:TR7col)=results.forget.merged.mean.NEW.object;
    allsub_RF_incorr(sub,f_restcol,TR1col:TR7col)=results.forget.merged.mean.NEW.rest;
    
    allsub_RF_incorr(sub,f_wordcol,TR1std:TR7std)=results.forget.merged.std.NEW.word;
    allsub_RF_incorr(sub,f_scenecol,TR1std:TR7std)=results.forget.merged.std.NEW.scene;
    allsub_RF_incorr(sub,f_facecol,TR1std:TR7std)=results.forget.merged.std.NEW.face;
    allsub_RF_incorr(sub,f_objcol,TR1std:TR7std)=results.forget.merged.std.NEW.object;
    allsub_RF_incorr(sub,f_restcol,TR1std:TR7std)=results.forget.merged.std.NEW.rest;
    
    
end


% Flatten matrix to 2D to print to excel\
allsub_corr=squeeze(allsub_RF_corr(1,:,:)); %squeeze removes singleton 1st dimension of sub.
allsub_incorr=squeeze(allsub_RF_incorr(1,:,:));

for sub=2:size(sublist,2)
    allsub_corr=vertcat(allsub_corr,squeeze(allsub_RF_corr(sub,:,:)));
    allsub_incorr=vertcat(allsub_incorr,squeeze(allsub_RF_incorr(sub,:,:)));
end



% output to excel
cd(outdir);
% % % xlsheader={'subject','Category_fsowr','Rem/For','Precue','Postcue','PostminPre','TR1','TR2','TR3','TR4','TR5','TR6','TR7','PrecueSTD','PostcueSTD','PostminPreSTD', 'TR1STD','TR2STD','TR3STD','TR4STD','TR5STD','TR6STd','TR7STD'};
% % % xlswrite(outname,xlsheader,'correct','A1');
% % % xlswrite(outname,allsub_corr,'correct','A2');
% % % xlswrite(outname,xlsheader,'incorrect','A1');
% % % xlswrite(outname,allsub_incorr,'incorrect','A2');

% save data in .MAT file
save(matname,'allsub_corr','allsub_incorr','-append');

%%  By Interation of Emotion by RF  %%%%%%%%%%%%%%%%%%%

clear allsub_RF_Neg;
clear allsub_RF_Neu;
allsub_RF_Neg=nan(size(sublist,2),10,23);
allsub_RF_Neu=nan(size(sublist,2),10,23);

for sub=1:size(sublist,2)
    clear results;
    
    subdir=strcat('emodif_', num2str(sublist(sub)),slashsym);
    datedir=datecell{datedir_ind(sub)};
    
    % load results struct from the specified phase
    % specify date by subject
    loaddir=strcat(datadir,subdir,'results',slashsym,phasedir,maskdir,datedir)
    
    load(strcat(loaddir,'results_DFencode_bar_compare.mat'));
    
    % Participant and condition info
    allsub_RF_Neg(sub,:,1)=sublist(sub);
    allsub_RF_Neg(sub,1:5,2)=[4,2,1,3,5]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
    allsub_RF_Neg(sub,6:10,2)=[4,2,1,3,5];
    
    allsub_RF_Neu(sub,:,1)=sublist(sub);
    allsub_RF_Neu(sub,1:5,2)=[4,2,1,3,5]; % stim category dummyk code 1=f, 2=s, 3=o, 4=w, 5=r;
    allsub_RF_Neu(sub,6:10,2)=[4,2,1,3,5];
    
    allsub_RF_Neg(sub,1:5,3)=1;  %remember
    allsub_RF_Neg(sub,6:10,3)=2; % forget
    
    allsub_RF_Neu(sub,1:5,3)=1;  %remember
    allsub_RF_Neu(sub,6:10,3)=2; % forget
    
    
    %%%%%%%%%%  Averaged TR bar graphs ( 1 number per data point) %%%%%%%%%%%%
    
    %%%%%%% Negative %%%%%%%%%%%%%%%
    % Remember Precue
    allsub_RF_Neg(sub,r_wordcol,precuecol)=results.remember.neg.bar.precue.word.nanmean;
    allsub_RF_Neg(sub,r_scenecol,precuecol)=results.remember.neg.bar.precue.scene.nanmean;
    allsub_RF_Neg(sub,r_facecol,precuecol)=results.remember.neg.bar.precue.face.nanmean;
    allsub_RF_Neg(sub,r_objcol,precuecol)=results.remember.neg.bar.precue.object.nanmean;
    allsub_RF_Neg(sub,r_restcol,precuecol)=results.remember.neg.bar.precue.rest.nanmean;
    
    allsub_RF_Neg(sub,r_wordcol,precuestd)=results.remember.neg.bar.precue.word.std;
    allsub_RF_Neg(sub,r_scenecol,precuestd)=results.remember.neg.bar.precue.scene.std;
    allsub_RF_Neg(sub,r_facecol,precuestd)=results.remember.neg.bar.precue.face.std;
    allsub_RF_Neg(sub,r_objcol,precuestd)=results.remember.neg.bar.precue.object.std;
    allsub_RF_Neg(sub,r_restcol,precuestd)=results.remember.neg.bar.precue.rest.std;
    
    
    % Remember Postcue
    allsub_RF_Neg(sub,r_wordcol,postcuecol)=results.remember.neg.bar.postcue.word.nanmean;
    allsub_RF_Neg(sub,r_scenecol,postcuecol)=results.remember.neg.bar.postcue.scene.nanmean;
    allsub_RF_Neg(sub,r_facecol,postcuecol)=results.remember.neg.bar.postcue.face.nanmean;
    allsub_RF_Neg(sub,r_objcol,postcuecol)=results.remember.neg.bar.postcue.object.nanmean;
    allsub_RF_Neg(sub,r_restcol,postcuecol)=results.remember.neg.bar.postcue.rest.nanmean;
    
    allsub_RF_Neg(sub,r_wordcol,postcuestd)=results.remember.neg.bar.postcue.word.std;
    allsub_RF_Neg(sub,r_scenecol,postcuestd)=results.remember.neg.bar.postcue.scene.std;
    allsub_RF_Neg(sub,r_facecol,postcuestd)=results.remember.neg.bar.postcue.face.std;
    allsub_RF_Neg(sub,r_objcol,postcuestd)=results.remember.neg.bar.postcue.object.std;
    allsub_RF_Neg(sub,r_restcol,postcuestd)=results.remember.neg.bar.postcue.rest.std;
    
    
    % % % % %     % Remember Post-Pre
    % % % % %     allsub_RF_Neg(sub,r_wordcol,postminprecol)=results.remember.neg.bar.postminuspre.word.nanmean;
    % % % % %     allsub_RF_Neg(sub,r_scenecol,postminprecol)=results.remember.neg.bar.postminuspre.scene.nanmean;
    % % % % %
    % % % % %
    % % % % %     allsub_RF_Neg(sub,r_wordcol,postminprestd)=results.remember.neg.bar.postminuspre.word.std;
    % % % % %     allsub_RF_Neg(sub,r_scenecol,postminprestd)=results.remember.neg.bar.postminuspre.scene.std;
    % % % % %
    
    
    % Forget Precue
    allsub_RF_Neg(sub,f_wordcol,precuecol)=results.forget.neg.bar.precue.word.nanmean;
    allsub_RF_Neg(sub,f_scenecol,precuecol)=results.forget.neg.bar.precue.scene.nanmean;
    allsub_RF_Neg(sub,f_facecol,precuecol)=results.forget.neg.bar.precue.face.nanmean;
    allsub_RF_Neg(sub,f_objcol,precuecol)=results.forget.neg.bar.precue.object.nanmean;
    allsub_RF_Neg(sub,f_restcol,precuecol)=results.forget.neg.bar.precue.rest.nanmean;
    
    allsub_RF_Neg(sub,f_wordcol,precuestd)=results.forget.neg.bar.precue.word.std;
    allsub_RF_Neg(sub,f_scenecol,precuestd)=results.forget.neg.bar.precue.scene.std;
    allsub_RF_Neg(sub,f_facecol,precuestd)=results.forget.neg.bar.precue.face.std;
    allsub_RF_Neg(sub,f_objcol,precuestd)=results.forget.neg.bar.precue.object.std;
    allsub_RF_Neg(sub,f_restcol,precuestd)=results.forget.neg.bar.precue.rest.std;
    
    % Forget Postcue
    allsub_RF_Neg(sub,f_wordcol,postcuecol)=results.forget.neg.bar.postcue.word.nanmean;
    allsub_RF_Neg(sub,f_scenecol,postcuecol)=results.forget.neg.bar.postcue.scene.nanmean;
    allsub_RF_Neg(sub,f_facecol,postcuecol)=results.forget.neg.bar.postcue.face.nanmean;
    allsub_RF_Neg(sub,f_objcol,postcuecol)=results.forget.neg.bar.postcue.object.nanmean;
    allsub_RF_Neg(sub,f_restcol,postcuecol)=results.forget.neg.bar.postcue.rest.nanmean;
    
    
    allsub_RF_Neg(sub,f_wordcol,postcuestd)=results.forget.neg.bar.postcue.word.std;
    allsub_RF_Neg(sub,f_scenecol,postcuestd)=results.forget.neg.bar.postcue.scene.std;
    allsub_RF_Neg(sub,f_facecol,postcuestd)=results.forget.neg.bar.postcue.face.std;
    allsub_RF_Neg(sub,f_objcol,postcuestd)=results.forget.neg.bar.postcue.object.std;
    allsub_RF_Neg(sub,f_restcol,postcuestd)=results.forget.neg.bar.postcue.rest.std;
    
    % % % % %
    % % % % %
    % % % % %     % Forget Post-Pre
    % % % % %     allsub_RF_Neg(sub,f_wordcol,postminprecol)=results.forget.neg.bar.postminuspre.word.nanmean;
    % % % % %     allsub_RF_Neg(sub,f_scenecol,postminprecol)=results.forget.neg.bar.postminuspre.scene.nanmean;
    % % % % %
    % % % % %     allsub_RF_Neg(sub,f_wordcol,postminprestd)=results.forget.neg.bar.postminuspre.word.std;
    % % % % %     allsub_RF_Neg(sub,f_scenecol,postminprestd)=results.forget.neg.bar.postminuspre.scene.std;
    % % % % %
    % % % % %
    %%%%%%%% Neutral %%%%%%%%%%%%%%%%%%%%%%%
    
    % Remember Precue
    allsub_RF_Neu(sub,r_wordcol,precuecol)=results.remember.neu.bar.precue.word.nanmean;
    allsub_RF_Neu(sub,r_scenecol,precuecol)=results.remember.neu.bar.precue.scene.nanmean;
    allsub_RF_Neu(sub,r_facecol,precuecol)=results.remember.neu.bar.precue.face.nanmean;
    allsub_RF_Neu(sub,r_objcol,precuecol)=results.remember.neu.bar.precue.object.nanmean;
    allsub_RF_Neu(sub,r_restcol,precuecol)=results.remember.neu.bar.precue.rest.nanmean;
    
    allsub_RF_Neu(sub,r_wordcol,precuestd)=results.remember.neu.bar.precue.word.std;
    allsub_RF_Neu(sub,r_scenecol,precuestd)=results.remember.neu.bar.precue.scene.std;
    allsub_RF_Neu(sub,r_facecol,precuestd)=results.remember.neu.bar.precue.face.std;
    allsub_RF_Neu(sub,r_objcol,precuestd)=results.remember.neu.bar.precue.object.std;
    allsub_RF_Neu(sub,r_restcol,precuestd)=results.remember.neu.bar.precue.rest.std;
    
    
    % Remember Postcue
    allsub_RF_Neu(sub,r_wordcol,postcuecol)=results.remember.neu.bar.postcue.word.nanmean;
    allsub_RF_Neu(sub,r_scenecol,postcuecol)=results.remember.neu.bar.postcue.scene.nanmean;
    allsub_RF_Neu(sub,r_facecol,postcuecol)=results.remember.neu.bar.postcue.face.nanmean;
    allsub_RF_Neu(sub,r_objcol,postcuecol)=results.remember.neu.bar.postcue.object.nanmean;
    allsub_RF_Neu(sub,r_restcol,postcuecol)=results.remember.neu.bar.postcue.rest.nanmean;
    
    allsub_RF_Neu(sub,r_wordcol,postcuestd)=results.remember.neu.bar.postcue.word.std;
    allsub_RF_Neu(sub,r_scenecol,postcuestd)=results.remember.neu.bar.postcue.scene.std;
    allsub_RF_Neu(sub,r_facecol,postcuestd)=results.remember.neu.bar.postcue.face.std;
    allsub_RF_Neu(sub,r_objcol,postcuestd)=results.remember.neu.bar.postcue.object.std;
    allsub_RF_Neu(sub,r_restcol,postcuestd)=results.remember.neu.bar.postcue.rest.std;
    
    
    % % % % %     % Remember Post-Pre
    % % % % %     allsub_RF_Neu(sub,r_wordcol,postminprecol)=results.remember.neu.bar.postminuspre.word.nanmean;
    % % % % %     allsub_RF_Neu(sub,r_scenecol,postminprecol)=results.remember.neu.bar.postminuspre.scene.nanmean;
    % % % % %
    % % % % %
    % % % % %     allsub_RF_Neu(sub,r_wordcol,postminprestd)=results.remember.neu.bar.postminuspre.word.std;
    % % % % %     allsub_RF_Neu(sub,r_scenecol,postminprestd)=results.remember.neu.bar.postminuspre.scene.std;
    % % % % %
    % % % % %
    
    % Forget Precue
    allsub_RF_Neu(sub,f_wordcol,precuecol)=results.forget.neu.bar.precue.word.nanmean;
    allsub_RF_Neu(sub,f_scenecol,precuecol)=results.forget.neu.bar.precue.scene.nanmean;
    allsub_RF_Neu(sub,f_facecol,precuecol)=results.forget.neu.bar.precue.face.nanmean;
    allsub_RF_Neu(sub,f_objcol,precuecol)=results.forget.neu.bar.precue.object.nanmean;
    allsub_RF_Neu(sub,f_restcol,precuecol)=results.forget.neu.bar.precue.rest.nanmean;
    
    allsub_RF_Neu(sub,f_wordcol,precuestd)=results.forget.neu.bar.precue.word.std;
    allsub_RF_Neu(sub,f_scenecol,precuestd)=results.forget.neu.bar.precue.scene.std;
    allsub_RF_Neu(sub,f_facecol,precuestd)=results.forget.neu.bar.precue.face.std;
    allsub_RF_Neu(sub,f_objcol,precuestd)=results.forget.neu.bar.precue.object.std;
    allsub_RF_Neu(sub,f_restcol,precuestd)=results.forget.neu.bar.precue.rest.std;
    
    % Forget Postcue
    allsub_RF_Neu(sub,f_wordcol,postcuecol)=results.forget.neu.bar.postcue.word.nanmean;
    allsub_RF_Neu(sub,f_scenecol,postcuecol)=results.forget.neu.bar.postcue.scene.nanmean;
    allsub_RF_Neu(sub,f_facecol,postcuecol)=results.forget.neu.bar.postcue.face.nanmean;
    allsub_RF_Neu(sub,f_objcol,postcuecol)=results.forget.neu.bar.postcue.object.nanmean;
    allsub_RF_Neu(sub,f_restcol,postcuecol)=results.forget.neu.bar.postcue.rest.nanmean;
    
    
    allsub_RF_Neu(sub,f_wordcol,postcuestd)=results.forget.neu.bar.postcue.word.std;
    allsub_RF_Neu(sub,f_scenecol,postcuestd)=results.forget.neu.bar.postcue.scene.std;
    allsub_RF_Neu(sub,f_facecol,postcuestd)=results.forget.neu.bar.postcue.face.std;
    allsub_RF_Neu(sub,f_objcol,postcuestd)=results.forget.neu.bar.postcue.object.std;
    allsub_RF_Neu(sub,f_restcol,postcuestd)=results.forget.neu.bar.postcue.rest.std;
    
    
    % % % % %
    % % % % %     % Forget Post-Pre
    % % % % %     allsub_RF_Neu(sub,f_wordcol,postminprecol)=results.forget.neu.bar.postminuspre.word.nanmean;
    % % % % %     allsub_RF_Neu(sub,f_scenecol,postminprecol)=results.forget.neu.bar.postminuspre.scene.nanmean;
    % % % % %
    % % % % %     allsub_RF_Neu(sub,f_wordcol,postminprestd)=results.forget.neu.bar.postminuspre.word.std;
    % % % % %     allsub_RF_Neu(sub,f_scenecol,postminprestd)=results.forget.neu.bar.postminuspre.scene.std;
    % % % % %
    % % % % %
    
    %%%%%%%%%  By TR timecourse graphs ( 1 by 7 timepoints) %%%%%%%%%
    
    %%%%%% Neg
    % Remember : Neg
    allsub_RF_Neg(sub,r_wordcol,TR1col:TR7col)=results.remember.neg.mean.word;
    allsub_RF_Neg(sub,r_scenecol,TR1col:TR7col)=results.remember.neg.mean.scene;
    allsub_RF_Neg(sub,r_facecol,TR1col:TR7col)=results.remember.neg.mean.face;
    allsub_RF_Neg(sub,r_objcol,TR1col:TR7col)=results.remember.neg.mean.object;
    allsub_RF_Neg(sub,r_restcol,TR1col:TR7col)=results.remember.neg.mean.rest;
    
    allsub_RF_Neg(sub,r_wordcol,TR1std:TR7std)=results.remember.neg.std.word;
    allsub_RF_Neg(sub,r_scenecol,TR1std:TR7std)=results.remember.neg.std.scene;
    allsub_RF_Neg(sub,r_facecol,TR1std:TR7std)=results.remember.neg.std.face;
    allsub_RF_Neg(sub,r_objcol,TR1std:TR7std)=results.remember.neg.std.object;
    allsub_RF_Neg(sub,r_restcol,TR1std:TR7std)=results.remember.neg.std.rest;
    
    
    
    
    % Forget Neg
    
    allsub_RF_Neg(sub,f_wordcol,TR1col:TR7col)=results.forget.neg.mean.word;
    allsub_RF_Neg(sub,f_scenecol,TR1col:TR7col)=results.forget.neg.mean.scene;
    allsub_RF_Neg(sub,f_facecol,TR1col:TR7col)=results.forget.neg.mean.face;
    allsub_RF_Neg(sub,f_objcol,TR1col:TR7col)=results.forget.neg.mean.object;
    allsub_RF_Neg(sub,f_restcol,TR1col:TR7col)=results.forget.neg.mean.rest;
    
    allsub_RF_Neg(sub,f_wordcol,TR1std:TR7std)=results.forget.neg.std.word;
    allsub_RF_Neg(sub,f_scenecol,TR1std:TR7std)=results.forget.neg.std.scene;
    allsub_RF_Neg(sub,f_facecol,TR1std:TR7std)=results.forget.neg.std.face;
    allsub_RF_Neg(sub,f_objcol,TR1std:TR7std)=results.forget.neg.std.object;
    allsub_RF_Neg(sub,f_restcol,TR1std:TR7std)=results.forget.neg.std.rest;
    
    
    
    %%%%%% Neutral
    % Remember Neutral
    allsub_RF_Neu(sub,r_wordcol,TR1col:TR7col)=results.remember.neu.mean.word;
    allsub_RF_Neu(sub,r_scenecol,TR1col:TR7col)=results.remember.neu.mean.scene;
    allsub_RF_Neu(sub,r_facecol,TR1col:TR7col)=results.remember.neu.mean.face;
    allsub_RF_Neu(sub,r_objcol,TR1col:TR7col)=results.remember.neu.mean.object;
    allsub_RF_Neu(sub,r_restcol,TR1col:TR7col)=results.remember.neu.mean.rest;
    
    allsub_RF_Neu(sub,r_wordcol,TR1std:TR7std)=results.remember.neu.std.word;
    allsub_RF_Neu(sub,r_scenecol,TR1std:TR7std)=results.remember.neu.std.scene;
    allsub_RF_Neu(sub,r_facecol,TR1std:TR7std)=results.remember.neu.std.face;
    allsub_RF_Neu(sub,r_objcol,TR1std:TR7std)=results.remember.neu.std.object;
    allsub_RF_Neu(sub,r_restcol,TR1std:TR7std)=results.remember.neu.std.rest;
    
    
    
    
    % forget Neutral
    
    allsub_RF_Neu(sub,f_wordcol,TR1col:TR7col)=results.forget.neu.mean.word;
    allsub_RF_Neu(sub,f_scenecol,TR1col:TR7col)=results.forget.neu.mean.scene;
    allsub_RF_Neu(sub,f_facecol,TR1col:TR7col)=results.forget.neu.mean.face;
    allsub_RF_Neu(sub,f_objcol,TR1col:TR7col)=results.forget.neu.mean.object;
    allsub_RF_Neu(sub,f_restcol,TR1col:TR7col)=results.forget.neu.mean.rest;
    
    allsub_RF_Neu(sub,f_wordcol,TR1std:TR7std)=results.forget.neu.std.word;
    allsub_RF_Neu(sub,f_scenecol,TR1std:TR7std)=results.forget.neu.std.scene;
    allsub_RF_Neu(sub,f_facecol,TR1std:TR7std)=results.forget.neu.std.face;
    allsub_RF_Neu(sub,f_objcol,TR1std:TR7std)=results.forget.neu.std.object;
    allsub_RF_Neu(sub,f_restcol,TR1std:TR7std)=results.forget.neu.std.rest;
    
    
end


% Flatten matrix to 2D to print to excel\
allsub_NegRF=squeeze(allsub_RF_Neg(1,:,:)); %squeeze removes singleton 1st dimension of sub.
allsub_NeuRF=squeeze(allsub_RF_Neu(1,:,:));

for sub=2:size(sublist,2)
    allsub_NegRF=vertcat(allsub_NegRF,squeeze(allsub_RF_Neg(sub,:,:)));
    allsub_NeuRF=vertcat(allsub_NeuRF,squeeze(allsub_RF_Neu(sub,:,:)));
end



% output to excel
cd(outdir);
xlsheader={'subject','Category_fsowr','Rem/For','Precue','Postcue','PostminPre','TR1','TR2','TR3','TR4','TR5','TR6','TR7','PrecueSTD','PostcueSTD','PostminPreSTD', 'TR1STD','TR2STD','TR3STD','TR4STD','TR5STD','TR6STd','TR7STD'};
xlswrite(outname,xlsheader,'Neg','A1');
xlswrite(outname,allsub_NegRF,'Neg','A2');
xlswrite(outname,xlsheader,'Neut','A1');
xlswrite(outname,allsub_NeuRF,'Neut','A2');

% save data in .MAT file
save(matname,'allsub_NegRF','allsub_NeuRF','-append');

