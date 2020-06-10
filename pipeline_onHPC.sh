#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH -C mc
#SBATCH --output=out_pipeline_desikan.txt
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=normal
#SBATCH --hint=nomultithread

export OMP_NUM_THREADS=36
module load /apps/daint/UES/easybuild/modulefiles/daint-mc
module load /apps/daint/system/modulefiles/shifter-ng/18.06.0

# input_dir, folder with BIDS data set
input_dir="/scratch/snx3000/bp000228/BIDS_dataset"

# output_dir, where to store results
output_dir="/scratch/snx3000/bp000228/test_desikan"

# participant_label, which subject to process from BIDS dataset
participant_label="QL20120814"

# parcellation for SC and FC
n="desikan"

# n_cpus, number of cpus for parallel processing
n_cpus=36

# freesurf_license, path to a freesurfer license for fmriprep
# put it into your $HOME 
freesurf_license=/users/bp000228/freesurfer_license/license.txt 

# name of task fmri, as in the BIDS_dataset
task_name="rest"

# --aroma-melodic-dimensionality
aroma_melodic_dimensionality=-120

# create subfolder for pipeline to store results
mrtrix_output=${output_dir}"/mrtrix_output"
fmriprep_output=${output_dir}"/fmriprep_output"
fmriprep_workdir=${fmriprep_output}"/tmp"
tvb_output=${output_dir}"/TVB_output"
tvb_workdir=${tvb_output}"/tmp"
mkdir -p ${output_dir} ${mrtrix_output} ${fmriprep_output} ${fmriprep_workdir} ${tvb_output} ${tvb_workdir}



# run mrtrix3 connectome to get a structural connectome
# the following command run the bids/mrtrix3_connectome
# a short workaround is implemented to overwrite the cleanup setting of mrtrix3_connectome.py script
# this way the freesurfer recon-all results are not deleted after the pipeline has finsihed
 srun shifter run --mount=type=bind,source=$HOME,destination=$HOME \
                --mount=type=bind,source=${input_dir},destination=/BIDS_dataset \
                --mount=type=bind,source=${mrtrix_output},destination=/mrtrix3_out \
                --writable-volatile=/home \
                bids/mrtrix3_connectome:latest python -c "from mrtrix3 import app; app.cleanup=False; \
                import sys; sys.argv='/mrtrix3_connectome.py /BIDS_dataset /mrtrix3_out participant \
                --participant_label ${participant_label} --parcellation ${parcellation} \
                --output_verbosity 2 --template_reg ants --n_cpus ${n_cpus} --debug \
                '.split(); execfile('/mrtrix3_connectome.py')"



# use the freesurfer recon-all results from mrtrix3_connectome container if available
# fmriprep will use them, when they are stored in the fmriprep_output dir with the subject name
# only these parcellations run recon-all in the mrtrix3_connectome pipeline
declare -a parcellations=("desikan" "destrieux" "hcpmmp1") 
if [[ " ${parcellations[*]} " == *"${parcellation}"* ]];
then
        # find recon-all results in temporary folder and copy them to fmriprep_output
        mrtrix_recon_all_dir=`find ${mrtrix_output} -name "mrtrix3_connectome.py*" -type d`
        mrtrix_recon_all_name="freesurfer"
        mkdir ${fmriprep_output}"/freesurfer"
        cp -r ${mrtrix_recon_all_dir}"/"${mrtrix_recon_all_name} ${fmriprep_output}"/freesurfer/sub-"${participant_label}
fi 

# run fmriprep to preprocess fmri data
srun shifter run --mount=type=bind,source=$HOME,destination=$HOME \
                --mount=type=bind,source=${input_dir},destination=/BIDS_dataset \
                --mount=type=bind,source=${fmriprep_output},destination=/fmriprep_out/ \
		 --mount=type=bind,source=${fmriprep_workdir},destination=/fmriprep_workdir/ \
                --writable-volatile=/home \
                poldracklab/fmriprep:latest /usr/local/miniconda/bin/fmriprep /BIDS_dataset /fmriprep_out/ participant \
                --use-aroma --bold2t1w-dof 6 --nthreads $SLURM_CPUS_PER_TASK --omp-nthreads $SLURM_CPUS_PER_TASK \
                --output-spaces T1w MNI152NLin6Asym:res-2 fsaverage5 --participant_label ${participant_label} \
                --fs-license-file ${freesurf_license} --aroma-melodic-dimensionality ${aroma_melodic_dimensionality} -w /fmriprep_workdir

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

weights_path=${mrtrix_output}"/sub-"${participant_label}"/connectome/sub-"${participant_label}"_parc-"${parcellation}"_level-participant_connectome.csv"
tracts_path=${mrtrix_output}"/sub-"${participant_label}"/connectome/sub-"${participant_label}"_parc-"${parcellation}"_meanlength.csv"
# workarounds, because of using shifter instead of docker
# mne.bem.make_watershed_bem seems to behave differently in shifter than in docker, therefore execute the py script inside the bem directory
# shifter won't let you mount single file into a directory, therefore mount the dir containing the freesurfer license file and than copy it from there into /opt/freesurfer
srun shifter run   --mount=type=bind,source=`dirname ${freesurf_license}`,destination=/freesurfer_license_dir/ \
		   --mount=type=bind,source=${output_dir},destination=/output_dir --mount=type=bind,source=${mrtrix_output},destination=/mrtrix3_out \
		   --mount=type=bind,source=${fmriprep_output},destination=/fmriprep_out --mount=type=bind,source=${fmriprep_workdir},destination=/fmriprep_workdir \
		   --mount=type=bind,source=${tvb_output},destination=/tvb_out --mount=type=bind,source=${tvb_workdir},destination=/tvb_workdir \
		   --mount=type=bind,source=${recon_all_dir},destination=/recon_all_dir \
	    	   --writable-volatile=/opt \
                   triebkjp/tvb_converter:latest /bin/bash -c "cp /freesurfer_license_dir/license.txt /opt/freesurfer/ && \
								mkdir -p /recon_all_dir/${recon_all_subject_name}/bem && \
								cd /recon_all_dir/${recon_all_subject_name}/bem && \
								/tvb_converter_pipeline.sh /output_dir /mrtrix3_out /fmriprep_out \
 						                /fmriprep_workdir /tvb_out /tvb_workdir /recon_all_dir ${recon_all_subject_name} \
						                ${participant_label} ${task_name} ${parcellation} ${weights_path} ${tracts_path} ${n_cpus}"

