% EmoDif Extract top X number ranking voxels from specified contrast
%
% Ref: https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/xSPM/johnsgems2/
%
% Judy Chiu
% 15/10/2018
%modified Tracy Wang 02/20/19



%% Use GUI to specify contrast RESULTS thresholds; then use below script to extract the ranked voxel table - for emodif this is .005

rankedroi.xSPM = xSPM;


% Full image size


% Get full Contrast table results for Scene> All (you have to pull this
% contrast up first
rankedroi.subj_dir=xSPM.swd;
rankedroi.contrast=xSPM.title;
rankedroi.contrastindex=xSPM.Ic;
rankedroi.stat=xSPM.STATstr;
rankedroi.heightthresh=xSPM.u;
rankedroi.voxelextent=xSPM.k;
rankedroi.sigvox=vertcat(xSPM.XYZmm,xSPM.XYZ,xSPM.Z);

% Rank each Subject's contrast results

% sort along the Z value and get index for sorting
zvals=xSPM.Z;
[sorted_Z,sort_ind]=sort(zvals,'descend');
rankedroi.voxrank=rankedroi.sigvox(:,sort_ind);


% Extract top ### voxels - for 100, 200 and 300

vox100=100;  % change for diff/ numbers of voxels to select;
vox200=200;
vox300=300;

rankedroi.voxselect300=rankedroi.voxrank(:,1:vox300);
rankedroi.voxselect200=rankedroi.voxrank(:,1:vox200);
rankedroi.voxselect100=rankedroi.voxrank(:,1:vox100);

voxel_sizes = [100 200 300];

%% Create Mask from XYZ coordinate 

% Extract a 3 by N matrix of our wanted Voxels
maskXYZ_100=rankedroi.voxselect100(1:3,:);  % grab the XYZmm information
maskXYZ_200=rankedroi.voxselect200(1:3,:);  % grab the XYZmm information
maskXYZ_300=rankedroi.voxselect300(1:3,:);  % grab the XYZmm information

% Extract header info of a functional image
Vin=spm_vol(spm_select(1,'image','select functional image'));
[RefY,refXYZ]=spm_read_vols(Vin);
mask=zeros(Vin.dim(1:3));

% loop over all voxels, set corresponding location in the mask to 1
mask100 = mask;
mask200 = mask;
mask300 = mask;
for v=1:size(maskXYZ_100,2)
    mask100(find(round(refXYZ(1,:),1)== round(maskXYZ_100(1,v),1)& round(refXYZ(2,:),1)== round(maskXYZ_100(2,v),1) & round(refXYZ(3,:),1)== round(maskXYZ_100(3,v),1)))=1;
end

for v=1:size(maskXYZ_200,2)
   
    mask200(find(round(refXYZ(1,:),1)== round(maskXYZ_200(1,v),1)& round(refXYZ(2,:),1)== round(maskXYZ_200(2,v),1) & round(refXYZ(3,:),1)== round(maskXYZ_200(3,v),1)))=1;
end

for v=1:size(maskXYZ_300,2)
  
    mask300(find(round(refXYZ(1,:),1)== round(maskXYZ_300(1,v),1)& round(refXYZ(2,:),1)== round(maskXYZ_300(2,v),1) & round(refXYZ(3,:),1)== round(maskXYZ_300(3,v),1)))=1;
end

cd(rankedroi.subj_dir);

% output as an image
Vout=Vin;
Vout.fname='Scene100.nii';
spm_write_vol(Vout,mask100);

Vout=Vin;
Vout.fname='Scene200.nii';
spm_write_vol(Vout,mask200);

Vout=Vin;
Vout.fname='Scene300.nii';
spm_write_vol(Vout,mask300);

foutname = 'rankedroiselection.mat';
save(foutname, 'rankedroi');

clear all





