python=$1
# Where to store intermediate files
intermediate_directory=$2
train_loader=./dataloader/train.pt

optimizer=${intermediate_directory}optimizer.pt
model=${intermediate_directory}model.pt
criterion=${intermediate_directory}criterion.pt

batches_dir=${intermediate_directory}batches/
labels_dir=${intermediate_directory}labels/

train_batch() {
    $python zero_grad.py $optimizer
    $python feed_batch_to_model.py $model $1 # model, batch
    $python calc_loss.py output.pt $2 $criterion #output, labels, criterion
    $python step_optimizer.py $optimizer
}

$python network_generator.py $intermediate_directory
echo Finished generating network elements

n_batches=`$python batchify.py $train_loader | tr -dc [:digit:]`
echo Finished exporting batches

for epoch in {1..15}
do   
    for i in $(seq 1 $n_batches)
    do
        train_batch ${batches_dir}batch_$i.pt ${labels_dir}labels_$i.pt
    done
done