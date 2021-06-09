FROM ubuntu:18.04

# some install steps are taken from bids/mrtrix3_connectome/dockerfile

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git gnupg2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# add the NeuroDebian repository to the native package management system, 
# than we can install neuroscience software the same way as any other package
RUN wget -O- http://neuro.debian.net/lists/bionic.de-m.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9 && apt-get update

# FSL installs tzdata as a dependency; however its installer is interactive
# Therefore we need to do some shenanigans here to force it though
RUN DEBIAN_FRONTEND=noninteractive \
    apt-get install -y tzdata

# install connectome workbench and fsl
RUN apt-get update && apt-get install -y connectome-workbench fsl

ENV FSLDIR=/usr/share/fsl/5.0
ENV PATH=$PATH:$FSLDIR/bin
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0:/usr/share/fsl/5.0/bin

#simulate . ${FSLDIR}/etc/fslconf/fsl.sh
ENV FSLBROWSER=/etc/alternatives/x-www-browser
ENV FSLCLUSTER_MAILOPTS=n
ENV FSLLOCKDIR=
ENV FSLMACHINELIST=
ENV FSLMULTIFILEQUIT=TRUE
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV FSLREMOTECALL=
ENV FSLTCLSH=/usr/bin/tclsh
ENV FSLWISH=/usr/bin/wish
ENV POSSUMDIR=/usr/share/fsl/5.0

# install conda3 from https://hub.docker.com/r/continuumio/miniconda3/dockerfile
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

RUN apt-get update && apt-get install vim libgl1-mesa-glx -y

# install mne
COPY mne_envorinment_wo_visuals.yml /mne_envorinment_wo_visuals.yml
RUN conda env update --file mne_envorinment_wo_visuals.yml
# Q: Why can't I install mayavi with the *.yml file above.
# A: I don't know why it gives an error. But this way works. 
RUN /bin/bash -c "source activate mne && pip install mayavi"

# install freesurfer
RUN wget -O- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz | \
    tar zx -C /opt \
    --exclude='freesurfer/trctrain' \
    --exclude='freesurfer/subjects/fsaverage_sym' \
    --exclude='freesurfer/subjects/fsaverage3' \
    --exclude='freesurfer/subjects/fsaverage4' \
    --exclude='freesurfer/subjects/fsaverage5' \
    --exclude='freesurfer/subjects/fsaverage6' \
    --exclude='freesurfer/subjects/cvs_avg35' \
    --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
    --exclude='freesurfer/subjects/bert' \
    --exclude='freesurfer/subjects/V1_average' \
    --exclude='freesurfer/average/mult-comp-cor' \
    --exclude='freesurfer/lib/cuda' \
    --exclude='freesurfer/lib/qt'

# Make FreeSurfer happy
ENV PATH=/opt/freesurfer/bin:/opt/freesurfer/mni/bin:$PATH
ENV OS Linux
ENV SUBJECTS_DIR /opt/freesurfer/subjects
ENV FSF_OUTPUT_FORMAT nii.gz
ENV MNI_DIR /opt/freesurfer/mni
ENV LOCAL_DIR /opt/freesurfer/local
ENV FREESURFER_HOME /opt/freesurfer
ENV FSFAST_HOME /opt/freesurfer/fsfast
ENV MINC_BIN_DIR /opt/freesurfer/mni/bin
ENV MINC_LIB_DIR /opt/freesurfer/mni/lib
ENV MNI_DATAPATH /opt/freesurfer/mni/data
ENV FMRI_ANALYSIS_DIR /opt/freesurfer/fsfast
ENV PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5
ENV MNI_PERL5LIB /opt/freesurfer/mni/lib/perl5/5.8.5


# add toolboxes for data privacy
RUN /bin/bash -c "pip install --upgrade pip"
RUN /bin/bash -c "pip install --upgrade pycryptodome"
RUN /bin/bash -c "pip install --upgrade pyAesCrypt"


# add scripts
COPY convert2TVB_format.py /convert2TVB_format.py
COPY encrypt_results.py /encrypt_results.py
COPY generateKeys.py /generateKeys.py
COPY encrypt_secret.py /encrypt_secret.py
COPY decrypt_data.py /decrypt_data.py
COPY tvb_converter.sh /tvb_converter.sh
COPY tvb_image_processing_pipeline.sh /tvb_image_processing_pipeline.sh
RUN chmod 775 /tvb_converter.sh
RUN chmod 775 /decrypt_data.py
RUN chmod 775 /encrypt_secret.py
RUN chmod 775 /generateKeys.py

# apt cleanup to recover as much space as possible
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


ENTRYPOINT [ "/tvb_converter.sh" ]
