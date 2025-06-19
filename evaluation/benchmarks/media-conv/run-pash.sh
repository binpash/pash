#!/bin/bash
cd $(dirname $0)

{ time ENTRIES=150 IN=media-conv/inputs/wav_full_heavy/wav/ OUT=media-conv/outputs/to-mp3-heavy-pash-core4-limit4/ SERVERLESS_PASH=1 $PASH_TOP/pa.sh -w1 --parallel_pipelines --parallel_pipelines_limit 4 scripts/to_mp3_heavy.sh ; } 2>outputs/to-mp3-heavy-pash-core4-limit4-time.log