#!/usr/bin/env bats

source $BATS_TEST_DIRNAME/run-demux-ipete.sh

@test "Testing sample workflow for complete-without-errors (demux-ipete)..." {
    run runpipeline demux-ipete
    [[ $status -eq 0 ]]
}
