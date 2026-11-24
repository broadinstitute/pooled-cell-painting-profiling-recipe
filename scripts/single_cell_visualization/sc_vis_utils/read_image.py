"""Read an image from an S3 bucket or from local disk, depending on the path."""

import skimage.io


def is_s3_uri(path):
    return str(path).startswith("s3://")


def read_image(path):
    """
    path: full path to the image. s3://bucket/key URIs are read from S3
    (see read_from_bucket.py); anything else is read from local disk.
    """
    if is_s3_uri(path):
        from sc_vis_utils import read_from_bucket

        return read_from_bucket.read_image(path)
    return skimage.io.imread(path)
