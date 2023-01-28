import sys
import pickle

with open(sys.argv[1], 'r+b') as file:
    classes = pickle.load(file)
    n_classes = len(classes)
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
        pickle.dump(classes, file) # only save if classes is changed

    print(n_classes)
    