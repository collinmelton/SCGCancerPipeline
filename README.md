# SCGCancerPipeline


## example run
`source ~/cancerPatients/pythonenvs/forDrmaa/bin/activate
python /srv/gsfs0/clinical/cancerPatientAnno/SCGCancerPipeline/JobManagementSoftware/PipelineWrapper_v2.py --P EP-88B --O /srv/gsfs0/scratch/cmelton/EP-88B/ --C /srv/gsfs0/clinical/cancerPatientAnno/Esplin/EP-88B/pr.csv --CFQ1 /srv/gsfs0/clinical/cancerPatientAnno/Esplin/FastqFiles/EP-88B/EP-88B_R1.fastq.gz --CFQ2 /srv/gsfs0/clinical/cancerPatientAnno/Esplin/FastqFiles/EP-88B/EP-88B_R2.fastq.gz --NFQ1 /srv/gsfs0/clinical/cancerPatientAnno/Esplin/FastqFiles/EP-107-NL/EP-107-NL_R1.fastq --NFQ2 /srv/gsfs0/clinical/cancerPatientAnno/Esplin/FastqFiles/EP-107-NL/EP-107-NL_R2.fastq --NL /srv/gsfs0/scratch/cmelton/EP-Adeno/EP-AdenoCa.merged.normal.RG.sorted.dedup.bam --RN 'F'`
