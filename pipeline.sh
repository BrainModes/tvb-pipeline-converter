#!/bin/bash
# input_dir, folder with BIDS data set
input_dir="/Users/andreas/Desktop/Brain_Modes/Projects/DATA_structure/BIDS_project/BIDS_test"

# output_dir, where to store results
output_dir="/Users/andreas/Desktop/Brain_Modes/MR_processing/Docker/output/test01"

# participant_label, which subject to process from BIDS dataset
participant_label="IS20120809"

# parcellation for SC and FC
parcellation="desikan"

# n_cpus, number of cpus for parallel processing
n_cpus=4
export OMP_NUM_THREADS=${n_cpus}
# mem_mb, memory for fmriprep in MB
mem_mb=15000

# freesurf_license, path to a freesurfer license for fmriprep
freesurf_license="/Users/andreas/Desktop/Brain_Modes/MR_processing/Docker/license.txt"

# name of task fmri, as in the BIDS_dataset
task_name="rest"

# create subfolder for pipeline to store results
mrtrix_output=${output_dir}"/mrtrix_output"
fmriprep_output=${output_dir}"/fmriprep_output"
fmriprep_workdir=${fmriprep_output}"/tmp"
tvb_output=${output_dir}"/TVB_output"
tvb_workdir=${tvb_output}"/tmp"
mkdir -p ${mrtrix_output} ${fmriprep_output} ${fmriprep_workdir} ${tvb_output} ${tvb_workdir}


###########################################################################################
###########################################################################################
###########################################################################################
# run mrtrix3 connectome to get a structural connectome
# the following command run the bids/mrtrix3_connectome
# a short workaround is implemented to overwrite the cleanup setting of mrtrix3_connectome.py script
# this way the freesurfer recon-all results are not deleted after the pipeline has finsihed
docker run -i -v ${input_dir}:${input_dir} \
        -v ${mrtrix_output}:${mrtrix_output} \
        --entrypoint python bids/mrtrix3_connectome -c "from mrtrix3 import app; app.cleanup=False; \
        import sys; sys.argv='/mrtrix3_connectome.py ${input_dir} ${mrtrix_output} participant \
        --participant_label ${participant_label} --parcellation ${parcellation} \
        --output_verbosity 3 --template_reg ants --n_cpus ${n_cpus} --debug \
        '.split(); execfile('/mrtrix3_connectome.py')"


###########################################################################################
###########################################################################################
###########################################################################################
# run fmriprep to preprocess fmri data
# use the freesurfer recon-all results from mrtrix3_connectome container if available
# only these parcellations run recon-all in the mrtrix3_connectome pipeline
declare -a parcellations=("desikan" "destrieux" "hcpmmp1") 
if [[ " ${parcellations[*]} " == *"${parcellation}"* ]];
then
        # find recon-all results in temporary folder and copy them to fmriprep_output
        mrtrix_recon_all_dir=`find ${mrtrix_output} -name "mrtrix3_connectome.py*" -type d`
        mrtrix_recon_all_name="freesurfer"
        mkdir ${fmriprep_output}"/freesufer"
        cp -r ${mrtrix_recon_all_dir}"/"${mrtrix_recon_all_name} ${fmriprep_output}"/freesurfer/sub-"${participant_label}
fi 


docker run -i --rm -v ${input_dir}:${input_dir} \
        -v ${fmriprep_output}:${fmriprep_output} \
        -v ${fmriprep_workdir}:${fmriprep_workdir} \
        -v ${freesurf_license}:/opt/freesurfer/license.txt \
        poldracklab/fmriprep:latest ${input_dir} \
        ${fmriprep_output} participant --bold2t1w-dof 6 --nthreads ${n_cpus} --omp-nthreads ${n_cpus} --use-aroma \
        --mem_mb ${mem_mb} --output-spaces T1w MNI152NLin6Asym:res-2 fsaverage5 --participant_label ${participant_label}  \
        -w ${fmriprep_workdir} --low-mem

###########################################################################################
###########################################################################################
###########################################################################################
# run TVBconverter to generate TVB format data, additionally compute the EEG leadfield matrix

# use the freesurfer folder of mrtrix3_connectome or fmriprep pipeline
if [ -z ${mrtrix_recon_all_dir+x} ]; 
then  # mrtrix_recon_all_dir is unset and recon-all has been performed within fmriprep
        recon_all_dir=${fmriprep_output}"/freesurfer/"
        recon_all_subject_name="sub-"${participant_label}
else # mrtrix_recon_all_dir is set
        recon_all_dir=${mrtrix_recon_all_dir}
        recon_all_subject_name=${mrtrix_recon_all_name}
        
fi

weighs_path=${mrtrix_output}"/sub-"${participant_label}"/connectome/sub-"${participant_label}"_parc-"${parcellation}"_level-participant_connectome.csv"
tracts_path=${mrtrix_output}"/sub-"${participant_label}"/connectome/sub-"${participant_label}"_parc-"${parcellation}"_meanlength.csv"

docker run -i --rm -v ${output_dir}:${output_dir} -v ${input_dir}:${input_dir} -v ${freesurf_license}:/opt/freesurfer/license.txt \
        thevirtualbrain/tvb_converter:latest ${input_dir} ${output_dir} ${mrtrix_output} ${fmriprep_output} ${fmriprep_workdir} \
        ${tvb_output} ${tvb_workdir} ${recon_all_dir} ${recon_all_subject_name} ${participant_label} ${task_name} \
        ${parcellation} ${weighs_path} ${tracts_path} ${n_cpus} 
