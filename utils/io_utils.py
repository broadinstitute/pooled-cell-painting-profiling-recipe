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
    batch = list(data_set.keys())[0]
    plate_list = list(data_set[batch].keys())
    list_plate_well_tuples = []
    for plate in plate_list:
        well_list = data_set[batch][plate]
        for well in well_list:
            list_plate_well_tuples.append((plate, well))
    return batch, plate_list, list_plate_well_tuples


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
