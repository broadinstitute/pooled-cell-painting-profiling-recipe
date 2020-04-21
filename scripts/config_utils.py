import pathlib
import yaml


def load_config(config, subset=None):
    data = {}
    with open(config, "r") as stream:
        for document in yaml.load_all(stream, Loader=yaml.FullLoader):
            data.update(document)
    return data


def get_step(config):
    # Determines which configuration to process
    data = load_config(config)
    return data["core"]["step"]


def make_batch_path(config, load_data=True):
    if load_data:
        data = load_config(config)
    else:
        data = config
    core = data["core"]
    batch_dir = pathlib.PurePath(f'{core["project_dir"]}/{core["data_dir"]}')
    return batch_dir


def preprocess_sites_config(config):
    data = load_config(config)
    batch_dir = make_batch_path(data, load_data=False)
    data["core"]["batch_dir"] = batch_dir
    return data


def generate_profiles_config(config):
    data = load_config(config)
    batch_dir = make_batch_path(data, load_data=False)
    data["core"]["batch_dir"] = batch_dir
    return data


def get_config(config):
    step = get_step(config)
    # Processes specified configuration file
    if step == "preprocess-sites":
        return preprocess_sites_config(config)
    elif step == "generate-profiles":
        return generate_profiles_config(config)
