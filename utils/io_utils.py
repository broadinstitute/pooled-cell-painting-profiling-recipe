import glob
import os
import yaml
import json
import pandas as pd


def load_configs(path_to_defaults_config, path_to_experiment_config):
    with open(path_to_defaults_config, "r") as stream:
        defaults = yaml.safe_load(stream)

    with open(path_to_experiment_config, "r") as stream:
        config = json.load(stream)
    return defaults, config


def parse_data_set(data_set):
    """
    All data must be in a single same batch
    """
    if len(data_set) != 1:
        raise ValueError(
            f"Each data set must contain exactly one batch, but found {len(data_set)}: {list(data_set.keys())}"
        )
    batch = list(data_set.keys())[0]
    plate_list = list(data_set[batch].keys())
    list_plate_well_tuples = []
    for plate in plate_list:
        well_list = data_set[batch][plate]
        for well in well_list:
            list_plate_well_tuples.append((plate, well))
    return batch, plate_list, list_plate_well_tuples


def site_folders_for_plate_well(folderlist, plate, well):
    """Site folders are named {plate}-{well}-{site}; match plate and well exactly."""
    return [x for x in folderlist if x.rsplit("-", 1)[0] == f"{plate}-{well}"]


def site_folders_for_plate(folderlist, plate):
    """Site folders are named {plate}-{well}-{site}; match plate exactly."""
    return [x for x in folderlist if x.rsplit("-", 2)[0] == plate]


def find_file_with_suffix(folder, suffix, exclude_suffix=None):
    """
    Find the file in folder whose name ends in suffix, allowing any prefix.
    Files ending in exclude_suffix are ignored (e.g. *BarcodeFoci.csv when looking for *Foci.csv).
    If nothing matches, returns folder/suffix so downstream load errors name a sensible file.
    """
    matches = glob.glob(os.path.join(folder, f"*{suffix}"))
    if exclude_suffix:
        matches = [x for x in matches if not x.endswith(exclude_suffix)]
    if len(matches) > 1:
        raise ValueError(f"Multiple files match *{suffix} in {folder}: {sorted(matches)}")
    if matches:
        return matches[0]
    return os.path.join(folder, suffix)


def read_csvs_with_chunksize(filename, chunksize=10000, **kwargs):
    """
    Read a CSV with an optionally passed chunksize to make reading large files easier.
    If Pandas ParserError, tries a second time.
    Re-raises any exceptions (likely mostly going to be FileNotFound errors) so they
    can continue to be handled how they currently are in the various locations.
    """
    try:
        with pd.read_csv(filename, chunksize=chunksize, **kwargs) as reader:
            dflist = []
            for chunk in reader:
                dflist.append(chunk)
            df = pd.concat(dflist)
        return df
    except pd.errors.ParserError or OSError:
        print(f"Error reading {filename}")
        try:
            with pd.read_csv(filename, chunksize=chunksize, **kwargs) as reader:
                dflist = []
                for chunk in reader:
                    dflist.append(chunk)
                df = pd.concat(dflist)
            print(f"Read {filename} on second try.")
            return df
        except:
            print(f"Error reading {filename} on second try.")
            raise
    except:
        raise
