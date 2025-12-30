#!/bin/bash

# Define the path to the preprocessing directory
preproc_dir="/Volumes/ExtremeSSDD/LIBR_tACS/Preprocessing_pre"

# Define the array of subjects
subjects=("AU880" "AO136" "AO931" "AP385" "AQ090" "AQ822" "AS590" "AY910" "BH088" "BI343" "BI346" "BI360" "BI392" "BI449" "BI460" "BI461" "BI491" "BI513" "BI550" "BI565" "BI632" "BI634" "BI699" "BI875" "BJ130" "BJ131" "BJ162" "BJ386" "BJ469" "BJ493" "BJ510" "BJ542" "BJ601" "BJ643" "BJ660" "BL307" "BL501" "BL749" "BL779" "BL884" "BM282" "BM419" "BN573" "BN643" "BP433" "BP446" "BP880" "BP917" "BQ069" "BQ850" "BQ976" "BR275" "BR276" "BR527" "BR665" "BR952" "BR953")



group_a_subjects=("AQ090" "BI343" "BI346" "BI460" "BI491" "BI513" "BI565" "BI632" "BI699" "BI875" "BJ131" "BJ386" "BJ469" "BJ510" "BJ643" "BJ660" "BL307" "BL501" "BL749" "BL884" "BM419" "BP433" "BP446" "BP880" "BP917" "BQ069" "BR275" "BR952" "BR953")

group_b_subjects=("AU880" "AO136" "AO931" "AP385" "AQ822" "AS590" "AY910" "BH088" "BI360" "BI392" "BI449" "BI461" "BI550" "BI634" "BJ130" "BJ162" "BJ493" "BJ542" "BJ601" "BJ612" "BL779" "BM282" "BN573" "BN643" "BQ850" "BQ976" "BR276" "BR527" "BR665")


# Initialize a variable to hold paths for the group analysis
group_paths=""

for subj in ${subjects[@]}; do
    session_id="MCR1"
    stats_file="${preproc_dir}/${subj}-${session_id}.results/stats.${subj}-${session_id}+tlrc.BRIK"
    
    # Check if the BRIK file exists without the sub-brick index
    if [ -e "${stats_file}" ]; then
        # Append the path with the sub-brick index for 3dttest++
        group_paths+="'${stats_file}[33]' "  # Correctly append the sub-brick index for AFNI
    else
        echo "BRIK file not found for subject ${subj}: ${stats_file}"
    fi
done

# Ensure there's no leading or trailing whitespace in group_paths
group_paths=$(echo $group_paths | xargs)

# Check if group_paths is empty
if [ -z "$group_paths" ]; then
    echo "No valid datasets were found. Check the file paths and permissions."
    exit 1
fi

# Proceed with 3dttest++
3dttest++ -prefix group_comparison_OvN_Total_Updated -setA $group_paths
