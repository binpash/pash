intermediate_directory=$1
train_loader=./dataloader/train.pt
python=`which python`

$python network_generator.py $intermediate_directory
echo Finished generating network elements

$python batchify.py $train_loader