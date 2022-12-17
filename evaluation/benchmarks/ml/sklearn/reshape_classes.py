def set_classes(self_classes):
    n_classes = len(self_classes)
    classes_ = self_classes
    if n_classes < 2:
        raise ValueError(
            "This solver needs samples of at least 2 classes"
            " in the data, but the data contains only one"
            " class: %r"
            % classes_[0]
        )

    if len(self_classes) == 2:
        n_classes = 1
        classes_ = classes_[1:]
    
    return n_classes, classes_