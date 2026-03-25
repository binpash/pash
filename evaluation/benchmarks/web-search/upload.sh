#!/bin/bash

cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

[ -z "$AWS_BUCKET" ] && {
  echo "AWS_BUCKET not set"
  exit
}

BENCHMARK_DIR="$PASH_TOP/evaluation/benchmarks/web-search"
INPUTS_DIR="$BENCHMARK_DIR/inputs"

S3_BUCKET_PREFIX="s3://$AWS_BUCKET"
S3_BENCHMARK_DIR="web-search"
S3_INPUTS_DIR="$S3_BENCHMARK_DIR/inputs"

# Upload stopwords
STOPWORDS="$INPUTS_DIR/stopwords.txt"
if [ -f "$STOPWORDS" ]; then
  S3_URI="$S3_BUCKET_PREFIX/$S3_INPUTS_DIR/stopwords.txt"
  echo "Uploading $STOPWORDS to $S3_URI"
  aws s3 cp "$STOPWORDS" "$S3_URI"
fi

# Upload index files
INDEX_FILES=(
  # index_min.txt
  index_small.txt
  # index.txt
)

for INDEX in "${INDEX_FILES[@]}"; do
  INDEX_PATH="$INPUTS_DIR/$INDEX"
  S3_URI="$S3_BUCKET_PREFIX/$S3_INPUTS_DIR/$INDEX"
  echo "Uploading $INDEX_PATH to $S3_URI"
  aws s3 cp "$INDEX_PATH" "$S3_URI"
done

# Upload article directories
ARTICLE_DIRS=(
  # articles_min
  articles_small
  # articles
)

for ARTICLES in "${ARTICLE_DIRS[@]}"; do
  ARTICLES_PATH="$INPUTS_DIR/$ARTICLES"
  S3_URI="$S3_BUCKET_PREFIX/$S3_INPUTS_DIR/$ARTICLES"
  echo "Uploading $ARTICLES_PATH to $S3_URI"
  aws s3 cp --recursive "$ARTICLES_PATH" "$S3_URI"
done
