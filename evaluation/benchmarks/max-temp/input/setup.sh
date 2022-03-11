#!/bin/bash
setup_dataset() {
  echo 'This experiment is expected to fetch data from a remote server'
  echo 'To fetch the original dataset, use an FTP client'
  echo 'e.g., "lftp ftp://ftp.ncdc.noaa.gov/pub/data/noaa"'
} 

source_var() {
  export IN=
}
