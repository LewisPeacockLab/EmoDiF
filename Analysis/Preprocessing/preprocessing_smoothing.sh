#!/bin/bash
#created twang.work@gmail.com to smooth emodif images by different smoothing kernals. FSL - fslmaths smooths in sigma - so conversion has to happen within the script. 

#subject number
#SUBNO = 101
#S_kernal = is in sigma  example 4mm = 1.70
SUBNO=$1
S_kernal='1.74'

# base directory where experiments are located on SCRATCH
BASEDIR='/Users/tw24955/emodif_data'
SOFTWAREDIR='/Users/tw24955/EmoDiF'
#makes functional average directory for coregistrations.


# name of the experiment
STUDYNAME='emodif'
STUDY_DATA='EmoDiF'

# code used for subjects when scanning
SUBCODE=${STUDYNAME}_${SUBNO}

#where the data is
SUBDIR=${BASEDIR}/${SUBCODE}

#functional data
FUNC_DIR=${SUBDIR}/BOLD

cd ${FUNC_DIR}

echo 'smoothing'

fslmaths Preview1_corr_mcf_brain.nii -s ${S_kernal} Preview1_corr_mcf_brain_sm4.nii
fslmaths Preview2_corr_mcf_brain.nii -s ${S_kernal} Preview2_corr_mcf_brain_sm4.nii
fslmaths DF_encoding_1_corr_mcf_brain.nii -s ${S_kernal} DF_encoding_1_corr_mcf_brain_sm4.nii
fslmaths DF_encoding_2_corr_mcf_brain.nii -s ${S_kernal} DF_encoding_2_corr_mcf_brain_sm4.nii
fslmaths MVPA_training_1_corr_mcf_brain.nii -s ${S_kernal} MVPA_training_1_corr_mcf_brain_sm4.nii
fslmaths MVPA_training_2_corr_mcf_brain.nii -s ${S_kernal} MVPA_training_2_corr_mcf_brain_sm4.nii
gunzip *

cd ${SOFWAREDIR}
