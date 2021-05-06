#!/usr/bin/env bats

source $BATS_TEST_DIRNAME/run-seg3-ipete.sh

@test "Testing sample workflow for complete-without-errors (segment3-ipete)..." {
    run runpipeline segment3-ipete
    [[ $status -eq 0 ]]
}
