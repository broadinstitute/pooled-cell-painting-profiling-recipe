import argparse


def parse_command_args(config_file):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config_file",
        help="configuration yaml file for preprocessing pipeline",
        default=config_file,
    )
    parser.add_argument(
        "--force", help="force overwriting of feature data", action="store_true"
    )
    args = parser.parse_args()
    return args
