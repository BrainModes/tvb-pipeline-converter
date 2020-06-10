#!/bin/bash
# some path to where results are stored
input_dir=${1}
output_dir=${2}
mrtrix_output=${3}
fmriprep_output=${4}
fmriprep_workdir=${5}
tvb_output=${6}
tvb_workdir=${7}

# recon all dir 
recon_all_dir=${8}
recon_all_name=${9}

# participant_label, which subject to process from BIDS dataset, i.e. sub-<participant_label>
participant_label=${10}

# name of task fmri, as in the BIDS_dataset
task_name=${11}

# parcellation for SC and FC
parcellation=${12}
weights_path=${13}
tracts_path=${14}

# number of cpus to use during parallel processing
n_cpus=${15}

# clean the fmri data with "nonaggressive" (i.e. only removing variance unique to the noise components) denoising using AROMA components
echo "Clean fMRI data"
fsl_regfilt -i ${fmriprep_output}/fmriprep/sub-${participant_label}/func/sub-${participant_label}_task-${task_name}_space-T1w_desc-preproc_bold.nii.gz \
    -f $(cat ${fmriprep_output}/fmriprep/sub-${participant_label}/func/sub-${participant_label}_task-${task_name}_AROMAnoiseICs.csv) \
    -d ${fmriprep_output}/fmriprep/sub-${participant_label}/func/sub-${participant_label}_task-${task_name}_desc-MELODIC_mixing.tsv \
    -o ${fmriprep_output}/fmriprep/sub-${participant_label}/func/sub-${participant_label}_task-${task_name}_space-T1w_AromaNonAggressiveDenoised.nii.gz

# resample parcellated image from mrtrix docker pipeline to fMRI resolution
flirt -in ${mrtrix_output}/sub-${participant_label}/anat/sub-${participant_label}_parc-${parcellation}_indices.nii.gz \
      -out ${mrtrix_output}/sub-${participant_label}/anat/sub-${participant_label}_parc-${parcellation}_indices_resample2bold.nii.gz \
      -ref ${fmriprep_output}/fmriprep/sub-${participant_label}/func/sub-${participant_label}_task-${task_name}_space-T1w_boldref.nii.gz \
      -nosearch -interp nearestneighbour

# extract fmri ROI timeseries using the parcellation
echo "Extract fmri ROI timeseries" 
fMRI_ROI_ts=${tvb_output}/sub-${participant_label}_task-${task_name}_parc-${parcellation}_ROI_timeseries.txt
fslmeants -o ${fMRI_ROI_ts} \
          -i ${fmriprep_output}/fmriprep/sub-${participant_label}/func/sub-${participant_label}_task-${task_name}_space-T1w_AromaNonAggressiveDenoised.nii.gz \
          --label=${mrtrix_output}/sub-${participant_label}/anat/sub-${participant_label}_parc-${parcellation}_indices_resample2bold.nii.gz


# deal with region mapping for different parcellations
# for parcellations defined only on the volume create the region map
# using wb_command 
declare -a parcellations=("aal" "aal2" "craddock200" "craddock400" "perry512")

if [[ " ${parcellations[*]} " == *"${parcellation}"* ]];
then
    echo "Create region map for parcellation "${parcellation}
    # from the parcellated image create a volume label image for wb_command,
    # with specifiying empty '' label-list-file
    wb_command -volume-label-import \
                ${mrtrix_output}/sub-${participant_label}/anat/sub-${participant_label}_parc-${parcellation}_indices.nii.gz \
                '' ${tvb_workdir}/sub-${participant_label}_parc-${parcellation}_indices_label_volume.nii.gz

    for Hemisphere in rh lh; do
        # get the surfaces in scanner space
        for surface in white pial; do
            mris_convert --to-scanner ${recon_all_dir}"/"${recon_all_name}"/surf/"${Hemisphere}"."${surface} \
                        ${tvb_workdir}"/"${Hemisphere}"."${surface}".surf.gii"
        done
         
        # NOTE: sometimes left and right cortical surfaces may overlap, especially in frontal areas
        # this will assign some region to both hemispheres
        # also nonlinear registration can result in a bad match, having some regions might reach over to the wrong hemisphere
        # can be crucial for surface simulations, TODO: Think of a fix for that !

        # do the mapping
        wb_command -volume-label-to-surface-mapping \
            ${tvb_workdir}/sub-${participant_label}_parc-${parcellation}_indices_label_volume.nii.gz \
            ${tvb_workdir}"/"${Hemisphere}".pial.surf.gii" \
            ${tvb_workdir}"/"${Hemisphere}"."${parcellation}".pial.label.gii" \
            -ribbon-constrained ${tvb_workdir}"/"${Hemisphere}".white.surf.gii" ${tvb_workdir}"/"${Hemisphere}".pial.surf.gii"
        
        # dilate if some vertices weren't assigned a label
        # this will assign "subcortical vertices" to cortical regions
        # but they will be filtered out afterwards during the python script
        wb_command -label-dilate ${tvb_workdir}"/"${Hemisphere}"."${parcellation}".pial.label.gii" \
                    ${tvb_workdir}"/"${Hemisphere}".pial.surf.gii" \
                    40 ${tvb_workdir}"/"${Hemisphere}"."${parcellation}".pial.label.gii"

        
    done
fi 


echo "Run the TVB converter py script."
source activate mne 
python /convert2TVB_format.py ${recon_all_name} ${recon_all_dir} ${mrtrix_output} ${participant_label} ${parcellation} \
                             ${mrtrix_output}/sub-${participant_label}/anat/sub-${participant_label}_parc-${parcellation}_indices.nii.gz \
                             ${tvb_output} ${tvb_workdir} ${n_cpus} ${weights_path} ${tracts_path} ${input_dir} ${fMRI_ROI_ts} ${task_name}


