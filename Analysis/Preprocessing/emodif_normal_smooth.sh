#!/bin/bash
# MODIFIED BY T.WANG twang.work@gmail.com 10.16.18 to normalise and smoothin in FSL # Usage: >./emodif 'imdif_study_epi_17'
#

# subject number
# SUBNO=101
#MIDRUN = name of the middle run (include localizer, preview  and study together)
#read in middle run (reference run) from the second argument in the script
#SUBNO=$1
SUBNO=$1
S_kernal='2.5' #translates to about 6mm FWHM

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

FUNC_DIR=${SUBDIR}/BOLD
ANAT_DIR=${SUBDIR}/anatomical

# NOT YET FIELDMAP_DIR= ${SUBDIR}/fieldmaps

FUNC_AVG_DIR=${FUNC_DIR}/avg_func_ref
mkdir ${FUNC_AVG_DIR}

##### VARIABLES THAT CAN CHANGE!!!!####
ST_TEMPLATE=${ST_MASK_DIR}/MNI152_T1_1mm_brain.nii
ST_TEMPLATE_HEAD=${ST_MASK_DIR}/MNI152_T1_1mm.nii
ST_TEMPLATE_MASK=${ST_MASK_DIR}/MNI152_T1_1mm_brain_mask_dil.nii
###########################################
#Create Middle Run Mean = this is your functional reference image fRI
midrun_bold=${FUNC_DIR}/DF_encoding_1_ref_mcf
pre_fRI=${FUNC_DIR}/DF_encoding_1_avg_mcf
fRI_head=${FUNC_AVG_DIR}/DF_encoding_1_avg_mcf
fRI=${FUNC_AVG_DIR}/DF_encoding_1_avg_mcf_brain
aRI=${ANAT_DIR}/T1_hires_brain
aRI_head=${ANAT_DIR}/T1_hires
func_STD=${ST_MASK_DIR}/MNI152_T1_2mm_brain.nii


#compute inverse transform (standard to MPRAGE)
MNI2T1=${FUNC_AVG_DIR}/MNI2T1 #set the fRI to MPRAGE mat file 
fRI2T1=${FUNC_AVG_DIR}/bold_co_avg_mcf_brain.mat
fRI2STD=${FUNC_AVG_DIR}/fRI2STD
T12fRI=${FUNC_AVG_DIR}/T12fRI
T12MNI=${FUNC_AVG_DIR}/T12MNI


cd ${FUNC_DIR}

# normalises all functions to MNI then smooths
echo 'normalizing Preview'

flirt -in Preview1_corr_mcf_brain.nii -ref ${func_STD} -out Preview1_corr_mcf_brain_mni.nii -omat Preview1_corr_mcf_brain_mni.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

applyxfm4d Preview1_corr_mcf_brain.nii Preview1_corr_mcf_brain_mni.nii Preview1_corr_mcf_brain_mni_4D Preview1_corr_mcf_brain_mni.mat -singlematrix

flirt -in Preview2_corr_mcf_brain.nii -ref ${func_STD} -out Preview2_corr_mcf_brain_mni.nii -omat Preview2_corr_mcf_brain_mni.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

applyxfm4d Preview2_corr_mcf_brain.nii Preview2_corr_mcf_brain_mni.nii Preview2_corr_mcf_brain_mni_4D Preview2_corr_mcf_brain_mni.mat -singlematrix

echo 'normalizing DFencode'
flirt -in DF_encoding_1_corr_mcf_brain.nii -ref ${func_STD} -out DF_encoding_1_corr_mcf_brain_mni.nii -omat DF_encoding_1_corr_mcf_brain_mni.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

applyxfm4d DF_encoding_1_corr_mcf_brain.nii DF_encoding_1_corr_mcf_brain_mni.nii DF_encoding_1_corr_mcf_brain_mni_4D DF_encoding_1_corr_mcf_brain_mni.mat -singlematrix

flirt -in DF_encoding_2_corr_mcf_brain.nii -ref ${func_STD} -out DF_encoding_2_corr_mcf_brain_mni.nii -omat DF_encoding_2_corr_mcf_brain_mni.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

applyxfm4d DF_encoding_2_corr_mcf_brain.nii DF_encoding_2_corr_mcf_brain_mni.nii DF_encoding_2_corr_mcf_brain_mni_4D DF_encoding_2_corr_mcf_brain_mni.mat -singlematrix

echo 'normalizing localizer'

flirt -in MVPA_training_1_corr_mcf_brain.nii -ref ${func_STD} -out MVPA_training_1_corr_mcf_brain_mni.nii -omat MVPA_training_1_corr_mcf_brain_mni.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

applyxfm4d MVPA_training_1_corr_mcf_brain.nii MVPA_training_1_corr_mcf_brain_mni.nii MVPA_training_1_corr_mcf_brain_mni_4D MVPA_training_1_corr_mcf_brain_mni.mat -singlematrix

flirt -in MVPA_training_2_corr_mcf_brain.nii -ref ${func_STD} -out MVPA_training_2_corr_mcf_brain_mni.nii -omat MVPA_training_2_corr_mcf_brain_mni.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

applyxfm4d MVPA_training_2_corr_mcf_brain.nii MVPA_training_2_corr_mcf_brain_mni.nii MVPA_training_2_corr_mcf_brain_mni_4D MVPA_training_2_corr_mcf_brain_mni.mat -singlematrix

echo 'smoothing Preview'
fslmaths Preview1_corr_mcf_brain_mni_4D -s ${S_kernal} Preview1_corr_mcf_brain_mni_4D_s6.nii
gunzip Preview1_corr_mcf_brain_mni_4D_s6.nii

fslmaths Preview2_corr_mcf_brain_mni_4D -s ${S_kernal} Preview2_corr_mcf_brain_mni_4D_s6.nii
gunzip Preview2_corr_mcf_brain_mni_4D_s6.nii

echo 'smoothing DF_encoding'
fslmaths DF_encoding_1_corr_mcf_brain_mni_4D -s ${S_kernal} DF_encoding_1_corr_mcf_brain_mni_4D_s6.nii
gunzip DF_encoding_1_corr_mcf_brain_mni_4D_s6.nii

fslmaths DF_encoding_2_corr_mcf_brain_mni_4D -s ${S_kernal} DF_encoding_2_corr_mcf_brain_mni_4D_s6.nii
gunzip DF_encoding_2_corr_mcf_brain_mni_4D_s6.nii

echo 'smoothing localizer'
fslmaths MVPA_training_1_corr_mcf_brain_4D_mni -s ${S_kernal} MVPA_training_1_corr_mcf_brain_mni_4D_s6.nii
gunzip MVPA_training_1_corr_mcf_brain_mni_4D_s6.nii

fslmaths MVPA_training_2_corr_mcf_brain_4D_mni -s ${S_kernal} MVPA_training_2_corr_mcf_brain_mni_4D_s6.nii
gunzip MVPA_training_2_corr_mcf_brain_mni_4D_s6.nii

#return to launch
cd ${SCRIPTDIR}


