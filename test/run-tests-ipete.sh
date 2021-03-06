parent=$(dirname $PWD)
export NXF_CLASSPATH=${parent}/nextflow/lib/

runpipeline() {
    CASEDIR=$BATS_TEST_DIRNAME/$1
    rm -rf "$CASEDIR/testoutput"
    mkdir "$CASEDIR/testoutput"
    cd "$CASEDIR/testoutput"
    curdir=$(pwd)
    cp ../SampleSheet.csv .

    cat ../config.txt | envsubst > ./config.txt

    nextflow run $BATS_TEST_DIRNAME/../nextflow/segment2-ipete.nf \
        -ansi-log false \
	    -profile ipete_docker \
	    -with-dag ipete.png \
	    -c config.txt \
	    > nextflow.log 2>&1 
}
