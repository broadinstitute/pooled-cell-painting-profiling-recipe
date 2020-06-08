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
    return data["master_config"]["step"]


def make_batch_path(config, load_data=True):
    if load_data:
        config = load_config(config)
    master = config["master_config"]
    core = config["core"]
    batch_dir = pathlib.Path(
        core["project_dir"], master["project_tag"], "workspace", "analysis", core["batch"],
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

    # Build visualization information
    if config["core"]["categorize_cell_quality"] == "simple":
        config["summarize-cells"]["cell_category_order"] = [
            "Perfect", "Great", "Imperfect", "Bad", "Empty"
        ]
        config["summarize-cells"]["cell_category_colors"] = [
            "#DB5F57", "#91DB57", "#57D3DB", "#A157DB", "#776244",
        ]
    elif config["core"]["categorize_cell_quality"] == "simple_plus":
        config["summarize-cells"]["cell_category_order"] = [
            "Perfect", "Great", "Imperfect-High", "Imperfect-Low", "Bad", "Empty"
        ]
        config["summarize-cells"]["cell_category_colors"] = [
            "#DB5F57", "#91DB57", "#57D3DB", "#556FD4", "#A157DB", "#776244",
        ]

    return config


def generate_profiles_config(config):
    config = load_config(config)
    batch_dir = make_batch_path(config, load_data=False)
    config["core"]["batch_dir"] = batch_dir
    batch = config["core"]["batch"]
    # Build paths to single cell yaml document of the profiling pipeline
    config["single_cell"]["prefilter_file"] = pathlib.Path(
        config["core"]["site_dir"], batch, "feature_prefilter.tsv"
    )

    config["single_cell"]["spot_metadata_dir"] = pathlib.Path(
        config["core"]["site_dir"], batch, "spots"
    )

    config["single_cell"]["paint_metadata_dir"] = pathlib.Path(
        config["core"]["site_dir"], batch, "paint"
    )

    config["single_cell"]["single_cell_output_dir"] = pathlib.Path(
        config["single_cell"]["output_basedir"], batch
    )

    # This file is only used if single_file_only flag is used in 0.merge-single-cells.py
    config["single_cell"]["single_file_only_output_file"] = pathlib.Path(
        config["single_cell"]["single_cell_output_dir"],
        f"{batch}_single_cell_profiles.csv.gz",
    )

    # Build paths to aggregate yaml document
    config["aggregate"]["aggregate_output_dir"] = pathlib.Path(
        config["aggregate"]["output_basedir"], batch
    )

    return config


def process_config_file(config):
    step = get_step(config)
    # Processes specified configuration file
    if step == "preprocess-sites":
        return preprocess_sites_config(config)
    elif step == "generate-profiles":
        return generate_profiles_config(config)
