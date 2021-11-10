module load nextflow_latest/19.07.0

parent=$(dirname $PWD)
export NXF_CLASSPATH=${parent}/nextflow/lib/

runpipeline() {
    CASEDIR=$BATS_TEST_DIRNAME/$1
    rm -rf "$CASEDIR/testoutput"
    mkdir "$CASEDIR/testoutput"
    cd "$CASEDIR/testoutput"
    curdir=$(pwd)

    cat ../config.txt | envsubst > ./config.txt

    nextflow run $BATS_TEST_DIRNAME/../nextflow/segment3-ipete.nf \
        -ansi-log false \
	    -profile ipete_docker \
	    -with-dag ipete.png \
	    -c config.txt \
	    > nextflow.log 2>&1 
}
