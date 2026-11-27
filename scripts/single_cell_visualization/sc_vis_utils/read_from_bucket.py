from io import BytesIO
from matplotlib import image as mpimg
import pandas as pd

_s3_client = None


def get_s3_client():
    """Create the S3 client on first use, so importing this module needs no AWS setup."""
    global _s3_client
    if _s3_client is None:
        import boto3

        # uses boto3's default credential chain (env vars, AWS_PROFILE, ~/.aws/credentials, etc.)
        _s3_client = boto3.client('s3')
    return _s3_client


def split_path(path):
    """Return (bucket, key) for a full s3://bucket/key URI."""
    if not path.startswith("s3://"):
        raise ValueError(f"Expected an s3://bucket/key URI, got {path!r}")
    bucket, key = path[len("s3://"):].split("/", 1)
    return bucket, key


def read_csv_gzip(path):

    bucket, key = split_path(path)
    data = get_s3_client().get_object(Bucket=bucket, Key=key)
    df_p_s = pd.read_csv(data["Body"], compression="gzip")

    return df_p_s


def read_image(path):

    bucket, key = split_path(path)
    key = key.replace('//','/') #sometimes the script making this accidentally has a //, but s3 paths should not
    response = get_s3_client().get_object(Bucket=bucket, Key=key)
    format_ext = path.split(".")[-1]
    image = mpimg.imread(BytesIO(response["Body"].read()), format=format_ext)

    return image
