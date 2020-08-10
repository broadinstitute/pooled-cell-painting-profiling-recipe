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
    return data["main_config"]["step"]


def get_output_path(config, load_data=True, mkdir=False, add_batch=False):
    if load_data:
        config = load_config(config)

    core = config["core"]
    output_basedir = pathlib.Path(core["output_basedir"])
    if add_batch:
        output_basedir = pathlib.Path(output_basedir / core["batch"])
    if mkdir:
        output_basedir.mkdir(exist_ok=True)
    return output_basedir


def make_input_batch_path(config, load_data=True):
    if load_data:
        config = load_config(config)
    main = config["main_config"]
    core = config["core"]
    input_batch_dir = pathlib.Path(
        core["project_dir"],
        main["project_tag"],
        "workspace",
        "analysis",
        core["batch"],
    )
    return input_batch_dir


def preprocess_sites_config(config):
    config = load_config(config)
    input_batch_dir = make_input_batch_path(config, load_data=False)
    output_basedir = get_output_path(config, load_data=False, add_batch=True)

    config["core"]["batch_dir"] = input_batch_dir

    # Build paths in the prefilter yaml document
    config["prefilter"]["prefilter_file"] = pathlib.Path(
        output_basedir, config["prefilter"]["output_dir"], "feature_prefilter.tsv",
    )

    config["prefilter"]["example_site_dir"] = pathlib.Path(
        input_batch_dir, config["prefilter"]["example_site"]
    )

    # Build paths in the process-spots yaml document
    config["process-spots"]["output_spotdir"] = pathlib.Path(
        output_basedir, config["process-spots"]["output_dir"],
    )

    config["process-spots"]["image_file"] = pathlib.Path(
        output_basedir,
        config["process-spots"]["image_output_dir"],
        "image_metadata.tsv",
    )

    config["process-cells"]["output_paintdir"] = pathlib.Path(
        output_basedir, config["process-cells"]["output_dir"],
    )

    config["summarize-cells"]["output_summary_resultsdir"] = pathlib.Path(
        output_basedir, config["summarize-cells"]["output_resultsdir"],
    )

    config["summarize-cells"]["output_summary_figuresdir"] = pathlib.Path(
        output_basedir, config["summarize-cells"]["output_figuresdir"],
    )

    config["summarize-cells"]["cell_count_file"] = pathlib.Path(
        config["summarize-cells"]["output_summary_resultsdir"], "cell_count.tsv"
    )
    config["summarize-cells"]["total_cell_count_file"] = pathlib.Path(
        config["summarize-cells"]["output_summary_resultsdir"], "total_cell_count.tsv"
    )

    # Build visualization information
    if config["core"]["categorize_cell_quality"] == "simple":
        config["summarize-cells"]["cell_category_order"] = [
            "Perfect",
            "Great",
            "Imperfect",
            "Bad",
            "Empty",
        ]
        config["summarize-cells"]["cell_category_colors"] = [
            "#DB5F57",
            "#91DB57",
            "#57D3DB",
            "#A157DB",
            "#776244",
        ]
    elif config["core"]["categorize_cell_quality"] == "simple_plus":
        config["summarize-cells"]["cell_category_order"] = [
            "Perfect",
            "Great",
            "Imperfect-High",
            "Imperfect-Low",
            "Bad",
            "Empty",
        ]
        config["summarize-cells"]["cell_category_colors"] = [
            "#DB5F57",
            "#91DB57",
            "#57D3DB",
            "#556FD4",
            "#A157DB",
            "#776244",
        ]

    return config


def generate_profiles_config(config, load_data=True):
    if load_data:
        config = load_config(config)

    input_batch_dir = make_input_batch_path(config, load_data=False)
    output_basedir = get_output_path(config, load_data=False, add_batch=True)

    config["core"]["batch_dir"] = input_batch_dir
    core = config["core"]
    batch = core["batch"]
    site_dir = core["site_dir"]
    ignore_files = core["ignore_files"]

    # Build paths to single cell yaml document of the profiling pipeline
    config["single_cell"]["prefilter_file"] = pathlib.Path(
        core["site_dir"], batch, core["prefilter_dir"], "feature_prefilter.tsv",
    )

    config["single_cell"]["single_cell_output_dir"] = pathlib.Path(
        output_basedir, config["single_cell"]["output_dir"]
    )

    config["single_cell"]["spot_metadata_dir"] = pathlib.Path(
        core["site_dir"], batch, core["spot_dir"]
    )
    config["single_cell"]["paint_metadata_dir"] = pathlib.Path(
        core["site_dir"], batch, core["paint_dir"]
    )

    # This file is only used if single_file_only flag is used in 0.merge-single-cells.py
    config["single_cell"]["single_file_only_output_file"] = pathlib.Path(
        config["single_cell"]["single_cell_output_dir"],
        f"{batch}_single_cell_profiles.csv.gz",
    )

    # Build single cell site files
    sites = [
        x.name
        for x in config["single_cell"]["spot_metadata_dir"].iterdir()
        if x.name not in ignore_files
    ]

    config["single_cell"]["site_files"] = {}
    for site in sites:
        # Define single cell output directory and files
        site_output_dir = pathlib.Path(
            config["single_cell"]["single_cell_output_dir"], site
        )
        config["single_cell"]["site_files"][site] = pathlib.Path(
            site_output_dir, f"{site}_single_cell.csv.gz"
        )

    # Build paths to aggregate yaml document
    config["aggregate"]["aggregate_output_dir"] = pathlib.Path(
        output_basedir, config["aggregate"]["output_dir"]
    )

    # Build aggregated output files
    config["aggregate"]["aggregate_output_files"] = {}
    for aggregate_level, aggregate_columns in config["aggregate"]["levels"].items():
        config["aggregate"]["aggregate_output_files"][aggregate_level] = pathlib.Path(
            config["aggregate"]["aggregate_output_dir"],
            f"{batch}_{aggregate_level}.csv.gz",
        )

    config["aggregate"]["aggregate_output_files"]["single_cell"] = config[
        "single_cell"
    ]["single_file_only_output_file"]

    # Build paths to normalize yaml document
    config["normalize"]["normalize_output_dir"] = pathlib.Path(
        output_basedir, config["normalize"]["output_dir"]
    )

    # Build normalized output files
    config["normalize"]["normalize_output_files"] = {}
    for normalize_level in config["normalize"]["levels"]:
        config["normalize"]["normalize_output_files"][normalize_level] = pathlib.Path(
            config["normalize"]["normalize_output_dir"],
            f"{batch}_{normalize_level}_normalized.csv.gz",
        )

    # Build paths to normalize yaml document
    config["feature_select"]["feature_select_output_dir"] = pathlib.Path(
        output_basedir, config["feature_select"]["output_dir"]
    )

    # Build feature select output files
    config["feature_select"]["feature_select_output_files"] = {}
    for feature_select_level in config["feature_select"]["levels"]:
        config["feature_select"]["feature_select_output_files"][
            feature_select_level
        ] = pathlib.Path(
            config["feature_select"]["feature_select_output_dir"],
            f"{batch}_{feature_select_level}_normalized_feature_select.csv.gz",
        )

    return config


def process_config_file(config):
    step = get_step(config)
    # Processes specified configuration file
    if step == "preprocess-sites":
        return preprocess_sites_config(config)
    elif step == "generate-profiles":
        return generate_profiles_config(config)
