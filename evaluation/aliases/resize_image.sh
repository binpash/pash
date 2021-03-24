# Resize image 
resize_image()
(
    cd $1
    find . -name "*.jpg" -o -name "*.png" -o -name "*.jpeg" -o -name "*.JPG" -o -name "*.PNG" \
        -o -name "*.JPEG" | xargs -P 16 -I {} sh -c "convert -resize 70% {} {}.70; mv {}.70 $2" sh {}
)

resize_image input/jpg $PWD/output/
