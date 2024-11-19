#!/bin/bash

cd "$(realpath $(dirname "$0"))"
IN="inputs"
mkdir -p inputs
cd inputs


if [ ! -d "wav" ]; then
    WAV_DATA_FILES=120
    mkdir wav_full
    wget https://atlas-group.cs.brown.edu/data/wav.zip -O wav.zip
    unzip wav.zip -d wav_full
    cd wav_full/wav
    for f in *.wav; do
        for (( i = 0; i <= $WAV_DATA_FILES; i++ )); do
            echo copying to $f$i.wav
            cp "$f" "$f$i.wav"
        done
    done
    cd -
    echo "WAV Generated"
fi

# Check if the wav_small directory does not exist
if [ ! -d "wav_small" ]; then
    WAV_DATA_FILES=20
    wget https://atlas-group.cs.brown.edu/data/wav.zip -O wav_small.zip
    mkdir wav_small
    unzip wav_small.zip -d wav_small
    cd wav_small/wav
    for f in *.wav; do
        for (( i = 0; i <= $WAV_DATA_FILES; i++ )); do
            echo copying to $f$i.wav
            cp "$f" "$f$i.wav"
        done
    done
    cd -
    echo "WAV_small Generated"
fi



# Check if the directories don't exist and handle full/jpg.zip
if [ ! -d "jpg" ]; then
    JPG_DATA_LINK=https://atlas-group.cs.brown.edu/data/full/jpg.zip
    wget $JPG_DATA_LINK -O jpg_full.zip
    unzip jpg_full.zip -d jpg_full
    echo "JPG Generated"
    rm -rf jpg_full.zip
fi

# Check if the small/jpg directory doesn't exist and handle small/jpg.zip
if [ ! -d "jpg_small" ]; then
    JPG_DATA_LINK=https://atlas-group.cs.brown.edu/data/small/jpg.zip
    wget $JPG_DATA_LINK -O jpg_small.zip
    unzip jpg_small.zip -d jpg_small
    echo "JPG_small Generated"
    rm -rf jpg_small.zip
fi
