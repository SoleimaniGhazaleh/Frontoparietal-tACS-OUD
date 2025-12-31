#!/bin/bash

subjects=("AU880" "AO136" "AO931" "AP385" "AQ090" "AQ822" "AS590" "AY910" "BH088" "BI343" "BI346" "BI360" "BI392" "BI449" "BI460" "BI461" "BI491" "BI513" "BI550" "BI565" "BI632" "BI634" "BI699" "BI875" "BJ130" "BJ131" "BJ162" "BJ386" "BJ469" "BJ493" "BJ510" "BJ542" "BJ601" "BJ612" "BJ643" "BJ660" "BL307" "BL501" "BL749" "BL779" "BL884" "BM282" "BM419" "BN573" "BN643" "BP433" "BP446" "BP880" "BP917" "BQ069" "BQ850" "BQ976" "BR275" "BR276" "BR527" "BR665" "BR952" "BR953")

subjectsTotal=("${subjects[@]}")

data_dir="/Volumes/ExtremeSSDD/LIBR_tACS/Pre"
stim_folder="/Volumes/ExtremeSSDD/LIBR_tACS/Preprocessing_pre/stim_folder_tACS"
session_id="MCR1"

# 1D files
opioid_1D="${stim_folder}/opioid-R1.1D"
neutral_1D="${stim_folder}/neutral-R1.1D"
yellow_1D="${stim_folder}/yellow-R1.1D"

for SUBJECT in ${subjectsTotal[@]}
do
    task_id=$SUBJECT-$session_id

    afni_proc.py -subj_id $task_id \
        -dsets $data_dir/$SUBJECT/func/sub-${SUBJECT}_ses-v1_task-opioidcuereactivity1_bold.nii.gz \
        -blocks despike tshift align tlrc volreg blur mask scale regress \
        -anat_uniform_method unifize \
        -despike_opts_3dDes -ignore 3 \
        -volreg_align_e2a \
        -volreg_tlrc_warp \
        -volreg_warp_dxyz 2.0 \
        -regress_censor_motion 0.3 \
        -regress_censor_outliers 0.1 \
        -tlrc_base /Volumes/ExtremeSSDD/LIBR_tACS/Template/MNI152_T1_2009c+tlrc \
        -copy_anat $data_dir/$SUBJECT/anat/sub-${SUBJECT}_ses-v1_run-1_T1w.nii.gz \
        -regress_stim_times  $opioid_1D  $neutral_1D  $yellow_1D \
        -regress_stim_labels O N Y \
        -regress_basis_multi 'BLOCK(31,1)' 'BLOCK(31,1)' 'BLOCK(5,1)' \
        -regress_opts_3dD -bout -nfirst 3 \
        -gltsym 'SYM: O -N' -glt_label 1 "OvN" \
        -execute \
        -remove_preproc_files
done
