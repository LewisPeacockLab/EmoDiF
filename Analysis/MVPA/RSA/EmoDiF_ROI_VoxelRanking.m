% EmoDif Extract top X number ranking voxels from specified contrast
%
% Ref: https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/spm/johnsgems2/
%
% Judy Chiu
% 15/10/2018

clear all; close all; clc;

computeruse=2;  %1= mac, 2=pc


%% Directory information
winscriptdir= 'C:\Users\ychiu4\Box Sync\CMF_BoxWork\Exp1_JudyMVPA\EmoDF_Scripts\EmoDF_LocGit'
windatadir='C:\Users\ychiu4\Box Sync\SharedFiles\EmoDF_LewPeaLab\emodif_data\'


%output folder location
winoutfold='C:\Users\ychiu4\Box Sync\CMF_BoxWork\Exp1_JudyMVPA\Pilot_Programs\JudyfMRI_Pilot\MR_PilotData\Analyzed\ROI_voxselect\';


if computeruse==1
elseif computeruse==2
    scriptdir=winscriptdir;
    datadir=windatadir;
    
else
    error 'input used computer; 1=mac, 2=pc'; 
end

%% Use GUI to specify contrast RESULTS thresholds; then use below script to extract the ranked voxel table

s101xSPM=xSPM;


% Full image size


% Get full Contrast table results for Scene> All 
rankedroi.subject='s101';
rankedroi.contrast=s101xSPM.title;
rankedroi.contrastindex=s101xSPM.Ic;
rankedroi.stat=s101xSPM.STATstr;
rankedroi.heightthresh=s101xSPM.u;
rankedroi.voxelextent=s101xSPM.k;
rankedroi.sigvox=vertcat(s101xSPM.XYZmm,s101xSPM.XYZ,s101xSPM.Z);

% Get full Contrast table results for WORD >Scene
rankedroi_word.subject='s101';
rankedroi_word.contrast=s101wordspm.title;
rankedroi_word.contrastindex=s101wordspm.Ic;
rankedroi_word.stat=s101wordspm.STATstr;
rankedroi_word.heightthresh=s101wordspm.u;
rankedroi_word.voxelextent=s101wordspm.k;
rankedroi_word.sigvox=vertcat(s101wordspm.XYZmm,s101wordspm.XYZ,s101wordspm.Z);



% Rank each Subject's contrast results

% sort along the Z value and get index for sorting
zvals=s101xSPM.Z;
[sorted_Z,sort_ind]=sort(zvals,'descend');
rankedroi.voxrank=rankedroi.sigvox(:,sort_ind);

word_zvals=s101wordspm.Z;
[sorted_wordZ,sort_word_ind]=sort(word_zvals,'descend');
rankedroi_word.voxrank=rankedroi_word.sigvox(:,sort_word_ind);


% Extracttop ### voxels
totvox=350;  % change for diff/ numbers of voxels to select;
if size(rankedroi.voxrank,2)>350
    totvox=350
else 
    totvox=size(ranked.voxrank(:,2))
end
rankedroi.voxselect=rankedroi.voxrank(:,1:totvox);


if size(rankedroi_word.voxrank,2)>350
    totvox_word=350
else
    totvox_word=size(rankedroi_word.voxrank,2)
end
ranakedroi_word.voxselect=rankedroi_word.voxrank(:,1:totvox_word);


%% Create Mask from XYZ coordinate 

% Extract a 3 by 350 matrix of our wanted Voxels
maskXYZ=rankedroi.voxrank(1:3,:);  % grab the XYZmm information

% Extract header info of a functional image
Vin=spm_vol(spm_select(1,'image','select functional image'));
[RefY,refXYZ]=spm_read_vols(Vin);
mask=zeros(Vin.dim(1:3));

% loop over all voxels, set corresponding location in the mask to 1
for v=1:size(maskXYZ,2)
    mask(find(round(refXYZ(1,:),1)== round(maskXYZ(1,v),1)& round(refXYZ(2,:),1)== round(maskXYZ(2,v),1) & round(refXYZ(3,:),1)== round(maskXYZ(3,v),1)))=1;
end

% output as an image
Vout=Vin;
Vout.fname='C:\Users\CMF_JC\Box Sync\SharedFiles\EmoDF_LewPeaLab\emodif_data\emodif_101\SPM_GLM\Scene350.img';
spm_write_vol(Vout,mask);

%% Word Mask %%%%%%%%%%%
maskXYZ_word=rankedroi_word.voxrank(1:3,:);

% Extract header info of a functional image
Vin_Word=spm_vol(spm_select(1,'image','select functional image'));
[RefY_word,refXYZ_word]=spm_read_vols(Vin_Word);
mask_Word=zeros(Vin_Word.dim(1:3));

% loop over all voxels, set corresponding location in the mask to 1
for v=1:size(maskXYZ_word,2)
    mask_Word(find(round(refXYZ_word(1,:),1)== round(maskXYZ_word(1,v),1)& round(refXYZ_word(2,:),1)== round(maskXYZ_word(2,v),1) & round(refXYZ_word(3,:),1)== round(maskXYZ_word(3,v),1)))=1;
end

% output as an image
Vout_Word=Vin_Word;
Vout_Word.fname='C:\Users\CMF_JC\Box Sync\SharedFiles\EmoDF_LewPeaLab\emodif_data\emodif_101\SPM_GLM\Word350.img';
spm_write_vol(Vout_Word,mask_Word);

