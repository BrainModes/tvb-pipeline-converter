#!/bin/bash

######## MAIN TVB PIPELINE ORCHESTRATOR SCRIPT
# Generic frontend to coordinate containers
# Authors: Michael Schirner 
#          Charité University Medicine, Berlin, Germany
#          Brain Simulation Section (PI: Petra Ritter)
version="1.0.0" # Sets version variable
# HISTORY:
#
# * Early 2020 - v1.0.0  - Initial Version: clean fMRI, convert to TVB format (developed by Paul Triebkorn)
# * Late 2020  - v2.0.0  - tvb_converter is now orchestrator of the pipeline and supports data  privacy
#
# NOTES:
#        The script uses bwrap to create a sandboxed bash and output directory.
#        Personal data is only decrypted from the sandboxed bash into the sanboxed folder.
#        Control is never given back to an unsandboxed process until computations that
#        inolve unencrypted personal data finish.
#
# For information on Shell script best practices, please see:
# https://github.com/ralish/bash-script-template/blob/stable/template.sh
# and
# https://natelandau.com/boilerplate-shell-script-template/
# ##################################################

# Enable xtrace if the DEBUG environment variable is set
if [[ ${DEBUG-} =~ ^1|yes|true$ ]]; then
    set -o xtrace       # Trace the execution of the script (debug)
fi

# A better class of script...
set -o errexit          # Exit on most errors (see the manual)
set -o errtrace         # Make sure any error trap is inherited
set -o nounset          # Disallow expansion of unset variables
set -o pipefail         # Use last non-zero exit code in a pipeline

declare -A secrets=()

log(){
    local msg="$1"
    timeAndDate=`date`
    echo "[$timeAndDate] $msg" >> $logFile
}

# DESC: Handler for unexpected errors
# ARGS: $1 (optional): Exit code (defaults to 1)
# OUTS: None
script_trap_err() {
    log "script_trap_err(): Trapped error. Pipeline will now exit."
    local exit_code=1

    # Disable the error trap handler to prevent potential recursion
    trap - ERR

    # Consider any further errors non-fatal to ensure we run to completion
    set +o errexit
    set +o pipefail

    # Remove tmp directory
    if [ ! -z ${tmpDir+x} ]; then
        if [ -d "${tmpDir}" ]; then
            rm -r "${tmpDir}"
            log "script_trap_err(): temporal directory removed: ${tmpDir}"
        fi
    fi

    # Validate any provided exit code
    if [[ ${1-} =~ ^[0-9]+$ ]]; then
        exit_code="$1"
    fi

    # Exit with failure status
    exit "$exit_code"
}

# DESC: Handler for exiting the script
# ARGS: None
# OUTS: None
script_trap_exit() {
    log "script_trap_exit(): Trapped exit. Pipeline will now exit."
    # Remove tmp directory
    if [ ! -z ${tmpDir+x} ]; then
        if [ -d "${tmpDir}" ]; then
            rm -r "${tmpDir}"
            log "script_trap_exit(): temporal directory removed: ${tmpDir}"
        fi
    fi

    cd "$orig_cwd"
}

# DESC: Exit script with the given message
# ARGS: $1 (required): Message to print on exit
#       $2 (optional): Exit code (defaults to 0)
# OUTS: None
# NOTE: The convention used in this script for exit codes is:
#       0: Normal exit
#       1: Abnormal exit due to external error
#       2: Abnormal exit due to script error
script_exit() {
    log "script_exit(): Pipeline will now exit."
    # Remove tmp directory
    if [ ! -z ${tmpDir+x} ]; then
        if [ -d "${tmpDir}" ]; then
            rm -r "${tmpDir}"
            log "script_exit(): temporal directory removed: ${tmpDir}"
        fi
    fi

    if [[ $# -eq 1 ]]; then
        printf '%s\n' "$1"
        exit 0
    fi

    if [[ ${2-} =~ ^[0-9]+$ ]]; then
        printf '%b\n' "$1"
        # If we've been provided a non-zero exit code run the error trap
        if [[ $2 -ne 0 ]]; then
            script_trap_err "$2"
        else
            exit 0
        fi
    fi

    script_exit 'Missing required argument to script_exit()!' 2
}



# DESC: Generic script initialisation
# ARGS: $@ (optional): Arguments provided to the script
# OUTS: $orig_cwd: The current working directory when the script was run
#       $script_path: The full path to the script
#       $script_dir: The directory path of the script
#       $script_name: The file name of the script
#       $script_params: The original parameters provided to the script
# NOTE: $script_path only contains the path that was used to call the script
#       and will not resolve any symlinks which may be present in the path.
#       You can use a tool like realpath to obtain the "true" path. The same
#       caveat applies to both the $script_dir and $script_name variables.
# shellcheck disable=SC2034
script_init() {
    # Useful paths
    readonly orig_cwd="$PWD"
    readonly script_path="${BASH_SOURCE[0]}"
    readonly script_dir="$(dirname "$script_path")"
    readonly script_name="$(basename "$script_path")"
    readonly script_params="$*"
    echo "script_init(): TVB Pipeline initialized"
}


# The Pipeline has 2 modes
# Mode 1: Pull containers
# Mode 2: Start daemon, generate keys, submit job, inject pass into running job
script_usage() {
    cat << EOF
The TVB Processing Pipeline has two usage modes:
Mode 1: Pull containers.
Mode 2: Start login-node daemon, generate input data keys.
Mode 3: Create and submit main job.
Mode 4: Main Job mode.

Usage Mode 1: Pull containers
    ./$script_name -m 1

    -m 1| Usage mode 1
    -h|   Displays this help

Usage Mode 2: Start login-node daemon, generate input data keys.
    ./$script_name -m 2 -i <path/to/input> -o <path/to/output>

    -m 2| Usage mode 2
    -i|   Path to input folder (BIDS format)
    -o|   Path to output folder
    -h|   Displays this help

Usage Mode 3: Create and submit main job.
    ./$script_name -m 3 -i <path/to/input> -o <path/to/output>

    -m 3| Usage mode 3
    -i|   Path to input folder (BIDS format)
    -o|   Path to output folder
    -h|   Displays this help

Usage Mode 4: Run main job.
    ./$script_name -m 4 -i <path/to/input> -o <path/to/output>

    -m 4| Usage mode 4
    -i|   Path to input folder (BIDS format)
    -o|   Path to output folder
    -h|   Displays this help

EOF
}



# DESC: Parameter parser
# ARGS: $@ (optional): Arguments provided to the script
# OUTS: Variables indicating command-line parameters and options
# The Pipeline has 2 modes
# Mode 1: Pull containers
# Mode 2: Start login-node daemon, generate keys for input, inject password into running job
# Mode 3: Create and submit batch for main job
# Mode 4: main job (runs on compute node)
# Mode 5: sandboxed main job (runs on compute node)
parse_params() {
    mode=-1
    input_path_flag=false
    output_path_flag=false

    # specify options and whether there are arguments expected (:)
    options='m:i:o:h'

    while getopts $options option
    do
        case "$option" in
            m  ) mode=$OPTARG;;
            i  ) input_path_flag=true; input_dir=$OPTARG;;
            o  ) output_path_flag=true; output_dir=$OPTARG;;
            h  ) script_usage; exit 0;;
            \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
            :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
            *  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
        esac
    done

    # No arguments were provided. Exit.
    if [ $# -lt 1 ]; then
        echo "No arguments provided. Aborting." >&2
        script_usage
        exit 1
    fi

    # check if mode number is valid (must be from 1-4)
    if ((mode > 5 || mode < 1)); then
        echo "Invalid usage mode "${mode}". Mode must be a number between 1 and 4." >&2
        script_usage
        exit 1
    fi

    # Make sure that input and output directories were specified and exist.
    if ((mode > 1)); then
        if ! $input_path_flag; then
            echo "Error: Usage mode "${mode}" was selected, but no input data path specified. Aborting." >&2
            script_usage
            exit 1
        fi
        if ! $output_path_flag; then
            echo "Error: Usage mode "${mode}" was selected, but no output data path specified. Aborting." >&2
            script_usage
            exit 1
        fi
        if [ ! -d "${input_dir}" ]; then
            echo "Input directory ${input_dir} does not exist. Aborting." >&2
            exit 1
        fi
        if [ ! -d "${output_dir}" ]; then
            echo "Output directory ${output_dir} does not exist. Aborting." >&2
            exit 1
        fi
    fi

    # Input argument parsing successful, start logging.
    log "******THE VIRTUAL BRAIN PROCESSING PIPELINE******"
    log "parse_params(): All arguments successfully parsed."
    log "parse_params(): Starting mode ${mode}..."
    if ((mode > 1)); then
        log "parse_params(): Input directory: ${input_dir}"
        log "parse_params(): Output directory: ${output_dir}"
    fi
}

# Generate RSA keys and stores them in two files. The private key file is protected by a password and encrypted with AES128-CBC and 'scrypt' to thwart dictionary attacks.
generate_keys() {
    log "generate_keys(): generating keys..."
    secrets[$1]="${RANDOM}_${RANDOM}_${RANDOM}_${RANDOM}" # This password must remain secret!

    # Load modules and run tvb_converter container to generate keys
    module load daint-mc
    module load sarus
#    srun -C mc --account=ich012 sarus run --mount=type=bind,source=${input_dir}/keys,destination=/keys michamischa/tvb-pipeline-converter:1.0 /tvb_converter.sh -k $passphrase -o /keys
#    sarus run --mount=type=bind,source=${input_dir}/keys,destination=/keys michamischa/tvb-pipeline-converter:1.0 /tvb_converter.sh -k $passphrase -o /keys

    sarus run --entrypoint "" --mount=type=bind,source=${input_dir}/keys,destination=/keys michamischa/tvb-pipeline-converter:1.2 python generateKeys.py ${secrets[$1]} /keys

#    ls -ltr ${input_dir}/keys
    mv ${input_dir}/keys/private_key.bin ${input_dir}/keys/private_key_${1}.bin
    mv ${input_dir}/keys/public_key.pem ${input_dir}/keys/public_key_${1}.pem

    log "generate_keys(): $1 keys generated."
}

slurm_header() {
    printf "#!/bin/bash -l\n#SBATCH --account=ich012\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --output=job-${1}.out\n#SBATCH --error=job-${1}.err\n#SBATCH --time=${2}\n#SBATCH --constraint=mc\nmodule load daint-mc\nmodule load sarus\n"
}

# Pull the three pipeline containers
pull_containers() {
    echo "pull_containers(): Generating batch files to pull containers." >&2
    log "pull_containers(): Generating batch files to pull containers."

    # specify header of slurm batch job files
    head1=$(slurm_header "SC" "05:00:00")
    head2=$(slurm_header "FC" "05:00:00")
    head3=$(slurm_header "Co" "05:00:00")

    # create slurm batch job files
    cat <<EOF > pull_c1.sh
${head1}
#srun sarus pull thevirtualbrain/tvb-pipeline-sc:1.0
srun sarus pull bids/mrtrix3_connectome
EOF

    cat <<EOF > pull_c2.sh
${head2}
srun sarus pull nipreps/fmriprep
#srun sarus pull thevirtualbrain/tvb-pipeline-fmriprep:1.0
EOF

    cat <<EOF > pull_c3.sh
${head3}
srun sarus pull michamischa/tvb-pipeline-converter:1.2
EOF


    echo "pull_containers(): Slurm batch files generated. Submitting jobs..." >&2
    log "pull_containers(): Slurm batch files generated. Submitting jobs..."

    sbatch pull_c1.sh
    sbatch pull_c2.sh
    sbatch pull_c3.sh

    echo "pull_containers(): Jobs submitted. Removing batch files." >&2
    log "pull_containers(): Jobs submitted. Removing batch files."

    rm pull_c1.sh
    rm pull_c2.sh
    rm pull_c3.sh

    exit 0
}


# Create main job batch scripts
submit_main_job() {
    echo "submit_main_job(): Generating batch files for main job." >&2
    log "submit_main_job(): Generating batch files for main job."

    # specify header of slurm batch job files
    head1=$(slurm_header "main" "00:05:00")

    # create slurm batch job files
    cat <<EOF > main_job.sh
${head1}
srun ${script_path} -m 4 -i ${input_dir} -o ${output_dir}
EOF

    sbatch main_job.sh

    echo "submit_main_job(): Jobs submitted. Removing batch files." >&2
    log "submit_main_job(): Jobs submitted. Removing batch files."

    rm main_job.sh
}


# Create Sandbox (currently unused)
# -----------------------------------
# Here, we spawn a Bash process that behaves exactly as outside the sandbox but
# additionally mounts a sandboxed temp directory.
# This directory will contain the unencrypted personal data.
# The temp directory contains three random numbers and the process ID
# in the name. This directory is removed automatically at exit.
# -----------------------------------
start_sandbox() {
    tmpDir="${output_dir}/tmp.$RANDOM.$RANDOM.$RANDOM.$$"
    (umask 077 && mkdir "${tmpDir}") || {
        echo "start_sandbox(): Could not generate temporal directory. Aborting." >&2
        log "start_sandbox(): Could not generate temporal directory. Aborting."
        exit 1
    }
    log "start_sandbox(): temporal directory created: ${tmpDir}"
    log "start_sandbox(): Spawning sandbox and re-starting pipeline."
    echo "start_sandbox(): Spawning sandbox: ${tmpDir}"
    bwrap --die-with-parent --dev-bind / / --tmpfs ${tmpDir} ${script_path} -m 5 -i ${input_dir} -o ${output_dir} -t ${tmpDir}
    echo "start_sandbox(): Sandbox returned. Removing tmp folder. Stopping."
    rm -rf ${tmpDir} # remove tmp dir (also in trap functions above)
    exit 0
}

start_sandbox_job() {
    tmpDir="${output_dir}/tmp.$RANDOM.$RANDOM.$RANDOM.$$"
    (umask 077 && mkdir "${tmpDir}") || {
        echo "start_sandbox_job(): Could not generate temporal directory. Aborting." >&2
        log "start_sandbox_job(): Could not generate temporal directory. Aborting."
        exit 1
    }
    log "start_sandbox_job(): Creating main pipeline job. (${tmpDir})"
    echo "start_sandbox_job(): Creating main pipeline job. (${tmpDir})"

    head=$(slurm_header "main" "00:10:00")
    cat <<EOF > main_job.sh
${head}
bwrap --die-with-parent --dev-bind / / --tmpfs ${tmpDir} ${script_path} -m ${mode} -i ${input_dir} -o ${output_dir} -s ${tmpDir}
EOF

    sbatch main_job.sh
    echo "start_sandbox_job(): Main pipeline job submitted. Stopping."
    log "start_sandbox_job(): Main pipeline job submitted. Stopping."
    exit 0
}



# This function lets the login-node daemon wait until
# the key from the compute node appears in the keys folder
# After 12 hours without sync the script exits.
sync_with_compute_node() {
    log "sync_with_compute_node(): Waiting for compute node..."
    echo "sync_with_compute_node(): Waiting for compute node..."
    SECONDS=0

    rm -f ${input_dir}/keys/public_key_compute.pem
    rm -f ${input_dir}/keys/encrypted_pwd_computenode.bin
    while ! test -f "${input_dir}/keys/public_key_compute.pem"; do
        sleep 5
        if (($SECONDS > 43200)); then
            log "sync_with_compute_node(): Synchronization failed after ${SECONDS} seconds. Stopping now."
            exit 1
        fi
    done
    log "sync_with_compute_node(): ...received public key."
    echo "sync_with_compute_node(): ...received public key."
}


# Here the login node daemon encrypts the password for the private key
# for decrypting the input data on the compute node using the public
# key that was just produced on the compute node.
encrypt_password() {
    sarus run --entrypoint "" --mount=type=bind,source=${input_dir}/keys,destination=/keys michamischa/tvb-pipeline-converter:1.2 python encrypt_secret.py ${secrets[$1]} /keys /keys/public_key_compute.pem

    log "encrypt_password(): Password for compute node encrypted."
    echo "encrypt_password(): Password for compute node encrypted."
}

# This function decrypts the freshly produced ${input_dir}/keys/private_key_input_pwd.bin
# which contains the password for decrypting the private key (private_key_input.bin)
# for decrypting the input data (encrypted_input_data.aes).
decrypt_data() {
    log "decrypt_data(): Waiting for login node to encrypt password..."
    echo "decrypt_data(): Waiting for login node to encrypt password..."
    SECONDS=0

    while ! test -f "${input_dir}/keys/private_key_input_pwd.bin"; do
        sleep 5
        if (($SECONDS > 43200)); then
            log "decrypt_data(): Synchronization failed after ${SECONDS} seconds. Stopping now."
            exit 1
        fi
    done
    log "decrypt_data(): ...received encrypted password."
    echo "decrypt_data(): ...received encrypted password."

    # The encrypted password for the input data private key arrived, now let's first decrypt this password,
    # then the private key file, and lastly the data.
    tmpDir=$(mktemp -d -p /scratch/snx3000/bp000225/)  || exit 1
    module load daint-mc
    module load sarus
    errormessage=$( sarus run --entrypoint "" --mount=type=bind,source=${input_dir},destination=/input --mount=type=bind,source=${tmpDir},destination=/data michamischa/tvb-pipeline-converter:1.2 python decrypt_data.py /input/keys/private_key_compute.bin ${secrets["compute"]} /input/keys/private_key_input_pwd.bin /input/keys/private_key_input.bin /input/encrypted_password.bin /input/input_data.zip.aes 2>&1 )


    echo $errormessage
    log $errormessage
    log "decrypt_data(): Data decrypted."
    echo "decrypt_data(): Data decrypted."

    cd ${tmpDir}
    unzip data.zip

    log "decrypt_data(): Data unzipped."
    echo "decrypt_data(): Data unzipped."
}


# This function starts the pipeline workflow with pipeline containers like MRtrix3_connectome, fmriprep, TVB pipeline converter
pipeline_workflow() {
    log "pipeline_workflow(): Started..."
    echo "pipeline_workflow(): Started..."

    module load daint-mc
    module load sarus

    mkdir ${tmpDir}/output
    # sarus run --mount=type=bind,source=${tmpDir},target=/BIDS_dataset --mount=type=bind,source=${tmpDir}/output,target=/output bids/mrtrix3_connectome /BIDS_dataset /output preproc

    # Add pipeline containers as needed
    # e.g. sarus run fmriprep etc.


    log "pipeline_workflow(): MRI processing finished."
    echo "pipeline_workflow(): MRI processing finished."
}


# This function encrypts the results with the public key from the data controllers computer after the
# processing of MRI data was finished.
encrypt_results() {
    log "encrypt_results(): Started..."
    echo "encrypt_results(): Started..."

    # compress results into archive
    tar -zcvf ${tmpDir}/results.tar.gz ${tmpDir}/output

    errormessage=$( sarus run --entrypoint "" --mount=type=bind,source=${input_dir}/keys,destination=/key --mount=type=bind,source=${tmpDir},destination=/input --mount=type=bind,source=${output_dir},destination=/output michamischa/tvb-pipeline-converter:1.2 python encrypt_results.py /keys/public_key_results.pem /input/results.tar.gz /output/results.tar.gz.aes /output/results_password.bin 2>&1 )

    echo $errormessage
    log $errormessage

    log "encrypt_results(): Results encrypted."
    echo "encrypt_results(): Results encrypted."
}

# This function deletes all unencrypted outputs
delete_unencrypted_output() {

    # remove the entire temporary folder
    rm -r ${tmpDir}

    log "delete_unencrypted_output(): Unencrypted results deleted."
    echo "delete_unencrypted_output(): Unencrypted results deleted."
}

# DESC: Main control flow
# ARGS: $@ (optional): Arguments provided to the script
# OUTS: None
main() {
    trap script_trap_err ERR
    trap script_trap_exit EXIT

    script_init "$@"
    readonly logFile="${script_dir}/${script_name}.log"
    parse_params "$@"

    if ((mode == 1)); then
        pull_containers
        exit 0
    fi

    if ((mode == 2)); then
        rm -rf "${input_dir}/keys"
        mkdir "${input_dir}/keys"
        generate_keys "input"
        sync_with_compute_node
        encrypt_password "input"
        exit 0
    fi

    if ((mode == 3)); then
        submit_main_job
        exit 0
    fi

    if ((mode == 4)); then
        #start_sandbox
        generate_keys "compute"
        decrypt_data
        pipeline_workflow
        encrypt_results
        delete_unencrypted_output
        exit 0
    fi
}

# Start script
main "$@"
exit 0
