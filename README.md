# tvb-pipeline-converter
Final leg of TVB processing pipeline. Converts outputs from other workflows into TVB and BIDS format.
Please see the [TVB PIPELINE EBRAINS Collab](https://wiki.ebrains.eu/bin/view/Collabs/tvb-pipeline/) for further reference on how to run the containerized converter on a Fenix site and the starter script [pipeline.sh](https://github.com/BrainModes/tvb-pipeline-converter/blob/master/pipeline.sh) for the call signature of the container.

To run the converter locally without containerization, the [convert2TVB_format.py](https://github.com/BrainModes/tvb-pipeline-converter/blob/master/convert2TVB_format.py) Python script can be used directly.

# References
Schirner, M., Domide, L., Perdikis, D., Triebkorn, P., Stefanovski, L., Pai, R., ... & Ritter, P. (2022). Brain simulation as a cloud service: The Virtual Brain on EBRAINS. NeuroImage, 118973.

# Acknowledgments

This  research  has  received  funding  from  the  European  Union’s  Horizon  2020  Framework  Programme  for  Research  and  Innovation  under  the  Specific  Grant  Agreement  Nos.  785907  (Human  Brain Project SGA2),  945539  (Human  Brain Project SGA3), ICEI 800858, VirtualBrainCloud 826421 and ERC 683049.
