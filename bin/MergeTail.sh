sampleName=$1

cmd=""
for i in *.tails.sam;
do
    cmd="$cmd I=$i"
done

/hgsc_software/java/latest/bin/java -Xmx64g -jar /hgsc_software/picard/picard-tools-1.84/MergeSamFiles.jar TMP_DIR=/space1/tmp/ VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true O=$sampleName.bam $cmd
