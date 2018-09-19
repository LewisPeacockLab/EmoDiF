#!/bin/bash
# MODIFIED BY T.WANG twang.work@gmail.com 1.12.16 to include all preprocessing and creation of hippocampal, amygdalar AND ventral temporal lobe masks
# Usage: >./emodif 'imdif_study_epi_17'
#

# subject number
# SUBNO=101
#MIDRUN = name of the middle run (include localizer, preview  and study together)
#read in middle run (reference run) from the second argument in the script
#SUBNO=$1
SUBNO='106'
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
ST_MASK_2=${ST_MASK_DIR}/tempoccfusi_pHg_LOC_combined_bin_MNI152.nii
ST_MASK_HIP=${ST_MASK_DIR}/JC_AMYHP_MNI.nii
ST_MASK=${ST_MASK_DIR}/tempoccfusi_pHg_combined_bin_MNI152.nii
ST_TEMPLATE=${ST_MASK_DIR}/MNI152_T1_1mm_brain.nii
ST_TEMPLATE_HEAD=${ST_MASK_DIR}/MNI152_T1_1mm.nii
ST_TEMPLATE_MASK=${ST_MASK_DIR}/MNI152_T1_1mm_brain_mask_dil.nii
MASK=${MASK_DIR}/tempoccfusi_pHg_combined_epi_space.nii
MASK_HIP=${MASK_DIR}/JC_AMYHP_epi_space.nii
MASK_2=${MASK_DIR}/tempoccfusi_pHg_LOC_combined_epi_space.nii
###########################################
#Create Middle Run Mean = this is your functional reference image fRI
midrun_bold=${FUNC_DIR}/DF_encoding_1_ref_mcf
pre_fRI=${FUNC_DIR}/DF_encoding_1_avg_mcf
fRI_head=${FUNC_AVG_DIR}/DF_encoding_1_avg_mcf
fRI=${FUNC_AVG_DIR}/DF_encoding_1_avg_mcf_brain
aRI=${ANAT_DIR}/T1_hires_brain
aRI_head=${ANAT_DIR}/T1_hires

slicetimecorr=${BASEDIR}/slicetiming

#compute inverse transform (standard to MPRAGE)
MNI2T1=${FUNC_AVG_DIR}/MNI2T1 #set the fRI to MPRAGE mat file 
fRI2T1=${FUNC_AVG_DIR}/bold_co_avg_mcf_brain.mat
fRI2STD=${FUNC_AVG_DIR}/fRI2STD
T12fRI=${FUNC_AVG_DIR}/T12fRI
T12MNI=${FUNC_AVG_DIR}/T12MNI


##### This is where the magic happens ####

# to ventral temporal mask
applywarp --ref=${fRI} --in=${ST_MASK_2} --warp=${MNI2T1} --postmat=${T12fRI} --out=${MASK_2}
# to hippocampus mask
#applywarp --ref=${fRI} --in=${ST_MASK_HIP} --warp=${MNI2T1} --postmat=${T12fRI} --out=${MASK_HIP}
#unzip mask
echo "unzipping subject-specific mask for MVPA decoding"
gunzip ${MASK}
#gunzip ${MASK_HIP}
gunzip ${FUNC_DIR}/*/*.gz

#return to launch
cd ${SCRIPTDIR}

