import pathlib
import pandas as pd


def get_compartment_file(compartment, example_dir):
    compartment = compartment.capitalize()
    compart_file = pathlib.PurePath(f"{example_dir}/{compartment}.csv")
    return compart_file


def load_compartments(core, example_dir):
    compartments = core["compartments"]

    data = {}
    for compartment in compartments:

        compart_file = get_compartment_file(compartment, example_dir)
        df = pd.read_csv(compart_file)
        df = recode_cols(df, core, compartment)

        data[compartment] = df

    return data


def recode_cols(df, core, compartment):
    df.columns = [f"{compartment}_{x}" for x in df.columns]

    rename_dict = {}
    recode_cols = [f"{compartment}_{x}" for x in core["id_cols"]]
    if compartment.lower() in core["parent_cols"]:
        recode_cols += [
            f"{compartment}_{x}" for x in core["parent_cols"][compartment.lower()]
        ]

    for recode_col in recode_cols:
        rename_dict[recode_col] = f"Metadata_{recode_col}"

    df = df.rename(rename_dict, axis="columns")

    features = df.columns.tolist()
    meta_cols = [features.index(x) for x in features if x.startswith("Metadata")]
    [features.insert(0, features.pop(x)) for x in meta_cols]

    df = df.reindex(features, axis="columns")

    return df


def load_features(core, example_dir):
    data = load_compartments(core, example_dir)

    features = {}
    for compartment in data:
        df = data[compartment]
        feature_df = pd.DataFrame(df.columns.tolist(), columns=["feature_name"]).assign(
            compartment=compartment
        )
        features[compartment] = feature_df

    feature_df = (
        pd.concat(features).sort_values(by="feature_name").reset_index(drop=True)
    )

    # Reorder rows for Metadata to exist first
    meta_df = feature_df.loc[feature_df.feature_name.str.contains("Metadata"), :]
    feature_df = feature_df.query("feature_name not in @meta_df.feature_name")
    feature_df = pd.concat([meta_df, feature_df]).reset_index(drop=True)

    return feature_df


def flag_features(df, flags):
    flag_cols = []
    for flag in flags:
        flag_cols += [x for x in df.feature_name if flag in x]

    flag_cols = list(set(flag_cols))
    df = df.assign(prefilter_column=df.feature_name.isin(flag_cols))
    return df


def prefilter_features(core, example_dir, flag_columns):
    feature_df = load_features(core, example_dir)
    feature_df = flag_features(feature_df, flag_columns)
    return feature_df
