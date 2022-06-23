from util import log

PATH_ANNOTATION_REPO=None

def get_path_annotation_repo():
    if PATH_ANNOTATION_REPO is None:
        log("No path for annotation repository given! Specify it in compiler/definitions/definition_path_for_annotation_repo.py")
        raise Exception("No path for annotation repository given! Specify it in compiler/definitions/definition_path_for_annotation_repo.py")
    return PATH_ANNOTATION_REPO