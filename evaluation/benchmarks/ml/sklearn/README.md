Running fit.sh will generate temporary files in a ./tmp folder

To parallelize, we want one-vs-rest classification, where we generate multiple models.
Additionally, the forest cover dataset has much more samples than it has features.
This makes the Newton-Cholesky solver ideal for this task.