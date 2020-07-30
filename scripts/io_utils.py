import pathlib
import warnings


def check_if_write(file_path, force, throw_warning=False):
    if file_path.exists():
        if force:
            if throw_warning:
                warnings.warn(f"{file_path} exists, overwriting")
            return True
        else:
            if throw_warning:
                warnings.warn(f"{file_path} exists, NOT overwriting")
            return False
    else:
        return True
