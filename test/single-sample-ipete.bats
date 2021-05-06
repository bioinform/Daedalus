#!/usr/bin/env bats

source $BATS_TEST_DIRNAME/run-tests-ipete.sh

@test "Testing sample workflow for complete-without-errors (single-sample-ipete)..." {
    run runpipeline single-sample-ipete
    [[ $status -eq 0 ]]
}
