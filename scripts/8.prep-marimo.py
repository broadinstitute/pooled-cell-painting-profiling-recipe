import os
import shutil
import logging
import traceback
import glob
import sys
import json

recipe_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(recipe_path, "utils"))
from io_utils import load_configs

# Configure logging
logfolder = os.path.join(recipe_path, "logs")
if not os.path.isdir(logfolder):
    os.mkdir(logfolder)
logging.basicConfig(
    filename=os.path.join(logfolder, "6.feature-select.log"),
    level=logging.INFO,
)


def handle_excepthook(exc_type, exc_value, exc_traceback):
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    traceback_details = "\n".join(traceback.extract_tb(exc_traceback).format())
    print(f"Uncaught Exception: {traceback_details}")


sys.excepthook = handle_excepthook


def printandlog(msg, type="info"):
    print(msg)
    if type == "warning":
        logging.warning(msg)
    else:
        logging.info(msg)

def rename_to_public_folder(root,filename,marimo_folder_loc):
    name_without_root = filename[len(root):]
    #print(name_without_root)
    name_without_root = name_without_root.replace('/','_')
    shutil.copy(filename,os.path.join(marimo_folder_loc,name_without_root))

def prep_marimo(path_to_defaults_config, path_to_experiment_config):
    printandlog("Starting 8.prep-marimo")
    defaults_config, experiment_config = load_configs(
        path_to_defaults_config, path_to_experiment_config
    )

    # define experiment variables
    data_sets = experiment_config["data_sets"]

    # define variables
    out_root = defaults_config["directory_structure"]["root"]
    qc_location = experiment_config["qc_report_loc"]

    marimo_folder_loc = os.path.join(recipe_path, 'scripts','public')

    if not os.path.exists(marimo_folder_loc):
        os.makedirs(marimo_folder_loc,exist_ok=True)

    image_recipe_glob = glob.glob(f"{out_root}/**/*.png",recursive=True)
    image_qcfile_glob = glob.glob(f"{qc_location}/**/*.png",recursive=True)
    profiles_glob = glob.glob(f"{out_root}/**/*_normalized_feature_selected_byall.csv.gz",recursive=True)
    mAP_glob = glob.glob(f"{out_root}/**/*_feature_selected_mAP_activity.csv.gz",recursive=True)
    stats_glob = glob.glob(f"{out_root}/**/Data_Stats*.json",recursive=True)

    all_recipe_globs = image_recipe_glob + profiles_glob + mAP_glob + stats_glob
    for each_glob in all_recipe_globs:
        rename_to_public_folder(out_root, each_glob, marimo_folder_loc)
    for each_glob in image_qcfile_glob:
        rename_to_public_folder(qc_location, each_glob, marimo_folder_loc)


    with open(os.path.join(marimo_folder_loc,'data_sets.json'),'w') as f:
        json.dump(data_sets,f)


#################################
# MAIN USER INTERACTION
#################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Pass paths to defaults config and experiment config.")
        sys.exit()
    path_to_defaults_config = sys.argv[1]
    path_to_experiment_config = sys.argv[2]
    prep_marimo(path_to_defaults_config, path_to_experiment_config)