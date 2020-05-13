import pathlib
import yaml


def load_config(config):
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
        config = load_config(config)
    core = config["core"]
    batch_dir = pathlib.Path(
        core["project_dir"], core["project"], "workspace", "analysis", core["batch"],
    )
    return batch_dir


def preprocess_sites_config(config):
    config = load_config(config)
    batch_dir = make_batch_path(config, load_data=False)
    config["core"]["batch_dir"] = batch_dir

    # Build paths in the prefilter yaml document
    config["prefilter"]["prefilter_file"] = pathlib.Path(
        config["prefilter"]["output_basedir"],
        config["core"]["batch"],
        "feature_prefilter.tsv",
    )
    config["prefilter"]["example_site_dir"] = pathlib.Path(
        batch_dir, config["prefilter"]["example_site"]
    )

    # Build paths in the process-spots yaml document
    config["process-spots"]["output_spotdir"] = pathlib.Path(
        config["process-spots"]["output_basedir"], config["core"]["batch"], "spots"
    )

    return config


def generate_profiles_config(config):
    config = load_config(config)
    batch_dir = make_batch_path(config, load_data=False)
    config["core"]["batch_dir"] = batch_dir
    return config


def process_config_file(config):
    step = get_step(config)
    # Processes specified configuration file
    if step == "preprocess-sites":
        return preprocess_sites_config(config)
    elif step == "generate-profiles":
        return generate_profiles_config(config)
