import sys
import pickle

with open(sys.argv[1], 'rb') as file:
    model = pickle.load(file)

n_classes = len(model.classes_)
classes = model.classes_
if n_classes < 2:
    raise ValueError(
        "This solver needs samples of at least 2 classes"
        " in the data, but the data contains only one"
        " class: %r"
        % classes[0]
    )
    # Maybe have exit here?
if len(classes) == 2:
    n_classes = 1
    classes = classes[1:]    

with open(sys.argv[2], 'w+b') as file:
    pickle.dump(classes, file)
print(n_classes)
    