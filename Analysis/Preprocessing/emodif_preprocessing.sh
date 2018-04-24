#!/bin/bash
# MODIFIED BY T.WANG twang.work@gmail.com 1.12.16 to include all preprocessing and creation of hippocampal, amygdalar AND ventral temporal lobe masks
# Usage: >./emodif 'imdif_study_epi_17'
#

# subject number
# SUBNO=101
#MIDRUN = name of the middle run (include localizer, preview  and study together)
#read in middle run (reference run) from the second argument in the script
#SUBNO=$1
SUBNO='101'

# base directory where experiments are located on SCRATCH
BASEDIR='/Users/tw24955/emodif_data'
SOFTWAREDIR='/Users/tw24955/EmoDiF'
#makes functional average directory for coregistrations.


# name of the experiment
STUDYNAME='emodif'
STUDY_DATA='EmoDiF'

# code used for subjects when scanning
SUBCODE=${STUDYNAME}_${SUBNO}

# CHANGING STUDYNAME BECAUSE OF PERMISSION ERROR
STUDYNAME=${STUDYNAME}

#SCRIPTDIR=${BASEDIR}/${STUDY_DATA}/batch/scripts
#mkdir -p ${SCRIPTDIR}

SUBDIR=${BASEDIR}/${SUBCODE}
# /scratch/03032/twang04/imdif/imdif_1506191/
#cd ${SUBDIR}

ST_MASK_DIR=${BASEDIR}/std_masks

MASK_DIR=${SUBDIR}/mask
mkdir -p ${MASK_DIR}
FUNC_DIR=${SUBDIR}/BOLD
ANAT_DIR=${SUBDIR}/anatomical

# NOT YET FIELDMAP_DIR= ${SUBDIR}/fieldmaps

FUNC_AVG_DIR=${FUNC_DIR}/avg_func_ref
mkdir ${FUNC_AVG_DIR}

##### VARIABLES THAT CAN CHANGE!!!!####

ST_MASK_HIP=${ST_MASK_DIR}/hippocampus_thr50_MNI152.nii 
ST_MASK=${ST_MASK_DIR}/tempoccfusi_pHg_combined_bin_MNI152.nii
ST_TEMPLATE=${ST_MASK_DIR}/MNI152_T1_1mm_brain.nii
ST_TEMPLATE_HEAD=${ST_MASK_DIR}/MNI152_T1_1mm.nii
ST_TEMPLATE_MASK=${ST_MASK_DIR}/MNI152_T1_1mm_brain_mask_dil.nii
MASK=${MASK_DIR}/tempoccfusi_pHg_combined_epi_space.nii
MASK_HIP=${MASK_DIR}/hippocampus_thr50_epi_space.nii
###########################################
#Create Middle Run Mean = this is your functional reference image fRI
midrun_bold=${FUNC_DIR}/DF_encoding_1_ref_mcf
pre_fRI=${FUNC_DIR}/DF_encoding_1_avg_mcf
fRI_head=${FUNC_AVG_DIR}/DF_encoding_1_avg_mcf
fRI=${FUNC_AVG_DIR}/DF_encoding_1_avg_mcf_brain
aRI=${ANAT_DIR}/T1_hires_brain
aRI_head=${ANAT_DIR}/T1_hires

slicetimecorr=${BASEDIR}/slicetiming

# take the bold.nii from the middle run, mcflirt and mean
echo movement correction on middle run
mcflirt -in ${FUNC_DIR}/DF_encoding_1 -out ${midrun_bold} -mats -plots;
echo fslmaths ${midrun_bold} -Tmean ${pre_fRI} 
fslmaths ${midrun_bold} -Tmean ${pre_fRI} 
echo copy middle run average to functional acerag directory and brain extraction
cp ${pre_fRI}.nii.gz ${FUNC_AVG_DIR}
cd ${FUNC_AVG_DIR}
bet DF_encoding_1_avg_mcf DF_encoding_1_avg_mcf_brain -R -m
gunzip ${FUNC_AVG_DIR}/*.gz

# motion correct (coreg) all functionals the to average of the middle run; brain extract and bptf high pass all runs

cd ${FUNC_DIR}

#scan_list=`ls ${FUNC_DIR} | grep imdif`
#for x in $scan_list; do cd ${FUNC_DIR}/$x; echo $x;
echo 'slice time correction'; slicetimer -i Preview1.nii -o Preview1_corr.nii -r 2; 
echo 'movement correction'; mcflirt -in Preview1_corr.nii -out Preview1_corr_mcf -reffile ${fRI_head} -mats -plots;
echo 'epi brain extraction'; bet Preview1_corr_mcf.nii Preview1_corr_mcf_brain.nii -F;
echo 'creating temporary mean'; fslmaths Preview1_corr_mcf_brain.nii -Tmean tempmean
echo 'high pass temporal filter';fslmaths Preview1_corr_mcf_brain.nii -bptf 32 -1 -add tempmean Preview1_corr_dt_mcf_brain.nii; rm tempmean.nii.gz; #128s high pass filter

echo 'slice time correction'; slicetimer -i Preview2.nii -o Preview2_corr.nii -r 2;
echo 'movement correction'; mcflirt -in Preview2_corr.nii -out Preview2_corr_mcf -reffile ${fRI_head} -mats -plots;
echo 'epi brain extraction'; bet Preview2_corr_mcf.nii Preview2_corr_mcf_brain.nii -F;
echo 'creating temporary mean'; fslmaths Preview2_corr_mcf_brain.nii -Tmean tempmean
echo 'high pass temporal filter';fslmaths Preview2_corr_mcf_brain.nii -bptf 32 -1 -add tempmean Preview2_corr_dt_mcf_brain.nii; rm tempmean.nii.gz; #128s high pass filter

echo 'slice time correction'; slicetimer -i DF_encoding_1.nii -o DF_encoding_1_corr.nii -r 2;
echo 'movement correction'; mcflirt -in DF_encoding_1_corr.nii -out DF_encoding_1_corr_mcf -reffile ${fRI_head} -mats -plots;
echo 'epi brain extraction'; bet DF_encoding_1_corr_mcf.nii DF_encoding_1_corr_mcf_brain.nii -F;
echo 'creating temporary mean'; fslmaths DF_encoding_1_corr_mcf_brain.nii -Tmean tempmean
echo 'high pass temporal filter';fslmaths DF_encoding_1_corr_mcf_brain.nii -bptf 32 -1 -add tempmean DF_encoding_1_corr_dt_mcf_brain.nii; rm tempmean.nii.gz; #128s high pass filter

echo 'slice time correction'; slicetimer -i DF_encoding_2.nii -o DF_encoding_2_corr.nii -r 2;
echo 'movement correction'; mcflirt -in DF_encoding_2_corr.nii -out DF_encoding_2_corr_mcf -reffile ${fRI_head} -mats -plots;
echo 'epi brain extraction'; bet DF_encoding_2_corr_mcf.nii DF_encoding_2_corr_mcf_brain.nii -F;
echo 'creating temporary mean'; fslmaths DF_encoding_2_corr_mcf_brain.nii -Tmean tempmean
echo 'high pass temporal filter';fslmaths DF_encoding_2_corr_mcf_brain.nii -bptf 32 -1 -add tempmean DF_encoding_2_corr_dt_mcf_brain.nii; rm tempmean.nii.gz; #128s high pass filter

echo 'slice time correction'; slicetimer -i MVPA_training_1.nii -o MVPA_training_1_corr.nii -r 2;
echo 'movement correction'; mcflirt -in MVPA_training_1_corr.nii -out MVPA_training_1_corr_mcf -reffile ${fRI_head} -mats -plots;
echo 'epi brain extraction'; bet MVPA_training_1_corr_mcf.nii MVPA_training_1_corr_mcf_brain.nii -F;
echo 'creating temporary mean'; fslmaths MVPA_training_1_corr_mcf_brain.nii -Tmean tempmean
echo 'high pass temporal filter';fslmaths MVPA_training_1_corr_mcf_brain.nii -bptf 32 -1 -add tempmean MVPA_training_1_corr_dt_mcf_brain.nii; rm tempmean.nii.gz; #128s high pass filter

echo 'slice time correction'; slicetimer -i MVPA_training_2.nii -o MVPA_training_2_corr.nii -r 2;
echo 'movement correction'; mcflirt -in MVPA_training_2_corr.nii -out MVPA_training_2_corr_mcf -reffile ${fRI_head} -mats -plots;
echo 'epi brain extraction'; bet MVPA_training_2_corr_mcf.nii MVPA_training_2_corr_mcf_brain.nii -F;
echo 'creating temporary mean'; fslmaths MVPA_training_2_corr_mcf_brain.nii -Tmean tempmean
echo 'high pass temporal filter';fslmaths MVPA_training_2_corr_mcf_brain.nii -bptf 32 -1 -add tempmean MVPA_training_2_corr_dt_mcf_brain.nii; rm tempmean.nii.gz; #128s high pass filter

echo 'unzipping all files'; gunzip *; 


#BET your structural with -R (recursive) and -B (neck removal)

#temporarily pre- skull stripped from anonymizing
echo 'bet ${aRI_head} ${aRI} -R -B'
bet ${aRI_head} ${aRI} -B -f 0.3

#Register your fRI to the anatomical scan (MPRAGE)

echo epi_reg --epi=${fRI} --t1=${aRI_head} --t1brain=${aRI} --out=${FUNC_AVG_DIR}/bold_co_avg_mcf_brain.nii.gz 
epi_reg --epi=${fRI} --t1=${aRI_head} --t1brain=${aRI} --out=${FUNC_AVG_DIR}/bold_co_avg_mcf_brain.nii.gz 


##### CREATE VENTRAL TEMPORAL MASK #####
#you already have the fRI to MPRAGE from epi_reg prior. 
#Coregister the MPRAGE to the standard (T1 MNI152)

cd ${FUNC_AVG_DIR}
echo "creating AFFINE transform for non-linear registration"
flirt -ref ${ST_TEMPLATE} -in ${aRI} -omat affineT12MNI_4fnirt.mat
echo "FNIRT in progress, non-linear registration of MPRAGE to MNI152_1mm.nii"
fnirt --ref=${ST_TEMPLATE_HEAD} --in=${aRI_head} --refmask=${ST_TEMPLATE_MASK} --aff=affineT12MNI_4fnirt.mat --cout=T12MNI --iout=highres_T12MNI_image  

#compute inverse transform (standard to MPRAGE)
MNI2T1=${FUNC_AVG_DIR}/MNI2T1
#set the fRI to MPRAGE mat file 
fRI2T1=${FUNC_AVG_DIR}/bold_co_avg_mcf_brain.mat
fRI2STD=${FUNC_AVG_DIR}/fRI2STD
T12fRI=${FUNC_AVG_DIR}/T12fRI
T12MNI=${FUNC_AVG_DIR}/T12MNI
 
echo "compute inverse transform MNI to T1"
invwarp --ref=${aRI_head} --warp=${T12MNI} --out=MNI2T1
#compute inverse transform (T1 to fRI)
echo "compute inverse transform T1 to fRI"
convert_xfm -omat T12fRI -inverse ${fRI2T1}
#concatenate both mat and warp files to achieve fRI to standard
echo "apply concateonated warps"
# to ventral temporal mask
applywarp --ref=${fRI} --in=${ST_MASK} --warp=${MNI2T1} --postmat=${T12fRI} --out=${MASK}
# to hippocampus mask
#applywarp --ref=${fRI} --in=${ST_MASK_HIP} --warp=${MNI2T1} --postmat=${T12fRI} --out=${MASK_HIP}
#unzip mask
echo "unzipping subject-specific mask for MVPA decoding"
gunzip ${MASK}
#gunzip ${MASK_HIP}
gunzip ${FUNC_DIR}/*/*.gz

#return to launch
cd ${SCRIPTDIR}
