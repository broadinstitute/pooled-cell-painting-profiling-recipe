import os
import sys
import argparse
import yaml
import json
import subprocess

from utils.io_utils import load_configs


def check_experiment_config(experiment_config):
    required_keys = [
        'file_location',
        'control_barcodes',
        'drop_barcodes',
        'data_sets',
        'compartments',
        'cell_quality_method',
        'keep_cell_qualities',
        'overwrite_files',
        'perform_process_qc',
        'perform_process_SBS',
        'perform_process_single_cells',
        'perform_summarize_SBS',
        'perform_aggregate',
        'perform_normalize',
        'perform_feature_select',
        'match_to_library',
        'library_structure',
        'library_location',
    ]
    for key in required_keys:
        if key not in experiment_config:
            raise ValueError(f"Missing required key '{key}' in experiment config")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pooled Phenotypic Screening Profiling")
    parser.add_argument(
        "--defaults-config",
        dest="defaults_config_path",
        default="config/defaults.yaml",
        help="Path to defaults config file",
    )
    parser.add_argument(
        "--experiment-config",
        dest="experiment_config_path",
        default="config/experiment.json",
        help="Path to experiment config file",
    )

    args = parser.parse_args()
    location = os.path.dirname(os.path.abspath(__file__))
    script_location = os.path.join(location, "scripts")

    defaults_config, experiment_config = load_configs(args.defaults_config_path, args.experiment_config_path)

    if experiment_config['perform_process_qc']:
        p = subprocess.Popen(
            ["python", os.path.join(script_location, "0.image-qc.py"), args.defaults_config_path, args.experiment_config_path],
        )
        p.communicate()
        if p.returncode != 0:
            print(f"Step 0.image-qc failed unexpectedly")
            sys.exit(1)

    if experiment_config['perform_process_SBS']:
        p = subprocess.Popen(
            ["python", os.path.join(script_location, "1.process-SBS.py"),args.defaults_config_path, args.experiment_config_path],
        )
        p.communicate()
        if p.returncode != 0:
            print(f"Step 1.process-SBS failed unexpectedly")
            sys.exit(1)

    if experiment_config['perform_process_single_cells']:
        p = subprocess.Popen(
            ["python", os.path.join(script_location, "2.merge-single-cells.py"),args.defaults_config_path, args.experiment_config_path],
        )
        p.communicate()
        if p.returncode != 0:
            print(f"Step 2.merge-single-cells failed unexpectedly")
            sys.exit(1)

    if experiment_config['perform_summarize_SBS']:
        p = subprocess.Popen(
            ["python", os.path.join(script_location, "3.summarize-SBS.py"),args.defaults_config_path, args.experiment_config_path],
        )
        p.communicate()
        if p.returncode != 0:
            print(f"Step 3.summarize-SBS failed unexpectedly")
            sys.exit(1)

    if experiment_config['perform_aggregate']:
        p = subprocess.Popen(
            ["python", os.path.join(script_location, "4.aggregate.py"),args.defaults_config_path, args.experiment_config_path],
        )
        p.communicate()
        if p.returncode != 0:
            print(f"Step 4.aggregate failed unexpectedly")
            sys.exit(1)

    if experiment_config['perform_normalize']:
        p = subprocess.Popen(
            ["python", os.path.join(script_location, "5.normalize.py"),args.defaults_config_path, args.experiment_config_path],
        )
        p.communicate()
        if p.returncode != 0:
            print(f"Step 5.normalize failed unexpectedly")
            sys.exit(1)

    if experiment_config['perform_feature_select']:
        p = subprocess.Popen(
            ["python", os.path.join(script_location, "6.feature-select.py"),args.defaults_config_path, args.experiment_config_path],
        )
        p.communicate()
        if p.returncode != 0:
            print(f"Step 6.feature_select failed unexpectedly")
            sys.exit(1)

    if experiment_config['perform_explore']:
        p = subprocess.Popen(
            ["python", os.path.join(script_location, "7.explore.py"),args.defaults_config_path, args.experiment_config_path],
        )
        p.communicate()
        if p.returncode != 0:
            print(f"Step 7.explore failed unexpectedly")
            sys.exit(1)