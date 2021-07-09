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


def read_csvs_with_chunksize(filename,chunksize=10000,**kwargs):
    """
    Read a CSV with an optionally passed chunksize to make reading large files easier.
    Re-raises any exceptions (likely mostly going to be FileNotFound errors) so they 
    can continue to be handled how they currently are in the various locations.
    """
    try:
        dflist = []
        for chunk in pd.read_csv(filename,chunksize=chunksize,**kwargs):
            dflist.append(chunk)
        df = pd.concat(dflist)
        return df
    except:
        raise