n_epoch=$1
intermediate_directory=$2
python=`which python`

optimizer=${intermediate_directory}optimizer.pt
model=${intermediate_directory}model.pt
criterion=${intermediate_directory}criterion.pt
output=${intermediate_directory}output.pt
loss=${intermediate_directory}loss.pt

batches_dir=${intermediate_directory}batches/
labels_dir=${intermediate_directory}labels/

scripts_dir=./python_pipeline/

n_batches=`ls -v $batches_dir | tail -1 | tr -dc [:digit:]`

train_batch() {
    $python ${scripts_dir}zero_grad.py $optimizer
    $python ${scripts_dir}feed_batch_to_model.py $model $1 $output
    $python ${scripts_dir}calc_loss.py $output $2 $criterion $loss
    $python ${scripts_dir}step_optimizer.py $optimizer
}

for epoch in $(seq 1 $n_epoch)
do   
    echo Epoch $epoch
    for i in $(seq 1 $n_batches)
    do
        echo $i
        train_batch ${batches_dir}batch_$i.pt ${labels_dir}labels_$i.pt
    done
done