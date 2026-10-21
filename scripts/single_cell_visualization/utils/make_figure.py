import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import skimage
import time

from utils import read_from_bucket

def CP_to_RGB_single(im_cp, channels, percentile=100):
    """
    This function takes a cell painting image (as channels last array) and converts it to RGB

    Inputs:
    ++ im_cp   (np array) size --> (width)x(height)x(channels):
    input dataframe contains single cells profiles as rows (make sure it has "Nuclei_Location_Center_X"or"Y" columns)

    ++ channels (dtype: list): list of channels to be displayed as columns of output image
           example: channels=['Mito','AGP','Brightfield','ER','DNA','Outline']

    Returns:
    colorImage0 (np array): dims-> width , height , 3 (RGB channels)

    """

    depth = 65535
    channels_colormap = {
        "DNA": "Blue",
        "ER": "Green",
        "RNA": "Yellow",
        "AGP": "Red",
        "Mito": "Magenta",
        "DAPI": "Blue",
        "WGA": "Red",
        "Phalloidin": "Green",
        "ZO1": "Magenta",
        "outline": "White"
    }

    channel_colors = [
        np.array(mcolors.to_rgb(channels_colormap[c])) * depth for c in channels
    ]

    comb_pars = [1 / im_cp.shape[2]] * im_cp.shape[2]
    colorImagesList = []
    for i in range(im_cp.shape[2]):
        image_gray = im_cp[:, :, i]
        image_gray = skimage.exposure.rescale_intensity(
                image_gray, in_range=(image_gray.min(), np.percentile(image_gray, percentile))
            )
        image_color = (
            (skimage.color.gray2rgb(image_gray).astype(float) / depth)
            * channel_colors[i]
            * comb_pars[i]
        )
        colorImagesList.append(image_color)

    colorImage0 = sum(colorImagesList)

    colorImage0 = skimage.exposure.rescale_intensity(
        colorImage0, out_range=(0, 255)
    ).astype(np.uint8)

    colorImagesList2 = [
        skimage.exposure.rescale_intensity(colim, out_range=(0, 255)).astype(np.uint8)
        for colim in colorImagesList
    ]

    return colorImage0, colorImagesList2

def crop_single_cell_image(image, xCenter, yCenter, halfBoxSize):

    im_h, im_w = image.shape
    before_y_pad = 0
    after_y_pad = 0
    before_x_pad = 0
    after_x_pad = 0

    if xCenter - halfBoxSize < 0:
        before_x_pad = abs(xCenter - halfBoxSize)

    if yCenter - halfBoxSize < 0:
        before_y_pad = abs(yCenter - halfBoxSize)

    if xCenter + halfBoxSize > im_w:
        after_x_pad = abs(im_w - xCenter - halfBoxSize)

    if yCenter + halfBoxSize > im_h:
        after_y_pad = abs(im_h - yCenter - halfBoxSize)

    image_cropped = image[
        np.maximum(yCenter - halfBoxSize, 0) : np.minimum(yCenter + halfBoxSize, im_h),
        np.maximum(xCenter - halfBoxSize, 0) : np.minimum(xCenter + halfBoxSize, im_w),
    ]
    if np.max([before_y_pad, after_y_pad, before_x_pad, after_x_pad]) > 0:
        image_cropped = np.pad(
            image_cropped,
            ((before_y_pad, after_y_pad), (before_x_pad, after_x_pad)),
            "minimum",
        )

    return image_cropped


def make_figure(gene, df, im_dict, cell_selection_method, box_size,
                objects, color, bc_nuclei, bc_nuclei_name, percentile
                ):
    im_channels = list(im_dict.keys())
    col_channels = list(im_dict.values())
    title=gene + "_" + cell_selection_method
    halfbox_size = int(box_size / 2)

    columns_count = len(im_channels) + 2*int(objects) + int(color) + int(bc_nuclei) + int(objects and color)
    rows_count = df.shape[0]

    f, axarr = plt.subplots(
        rows_count, columns_count, figsize=(columns_count * 2, rows_count * 2)
    )
    f.suptitle(title)
    f.subplots_adjust(hspace=0, wspace=0)

    for index in range(rows_count):

        xCenter = xCenterC = int(df.loc[index, "Nuclei_AreaShape_Center_X"])
        yCenter = yCenterC = int(df.loc[index, "Nuclei_AreaShape_Center_Y"])

        sc_collage_row = np.zeros((box_size, box_size, columns_count))

        num_channel_cols = len(im_channels)
        clim_max = []
        for ci in range(num_channel_cols):
            if index == 0:
                axarr[index, ci].set_title(im_channels[ci])

            ch_fName = df.loc[index, "FileName_Corr" + im_channels[ci]]
            ch_pName = df.loc[index, "PathName_Corr" + im_channels[ci]]
            
            xCenterC = xCenter + int(
                df.loc[
                    index, "Align_Xshift_" + col_channels[ci]]
            )
            yCenterC = yCenter + int(
                df.loc[
                    index, "Align_Yshift_" + col_channels[ci]]
            )

            imPath = os.path.join(ch_pName,ch_fName)
            image = np.squeeze(read_from_bucket.read_image(imPath))

            # make images brighter for visibility
            image = skimage.exposure.rescale_intensity(
                image, in_range=(image.min(), np.percentile(image, percentile))
            )

            clim_max.append(image.max())
            image_cropped = crop_single_cell_image(
                image, xCenterC, yCenterC, halfbox_size
            )

            sc_collage_row[:, :, ci] = image_cropped

        if bc_nuclei:
            if index == 0:
                axarr[index, num_channel_cols].set_title("SBS Nuclei")
            # different naming, don't need alignment shift
            ch_fName = df.loc[index, "FileName_" + bc_nuclei_name]
            ch_pName = df.loc[index, "PathName_" + bc_nuclei_name]
            imPath = ch_pName + "/" + ch_fName
            image = np.squeeze(read_from_bucket.read_image(imPath))
            # make images brighter for visibility
            image = skimage.exposure.rescale_intensity(
                image, in_range=(image.min(), np.percentile(image, percentile))
            )

            clim_max.append(image.max())
            image_cropped = crop_single_cell_image(
                image, xCenter, yCenter, halfbox_size
            )

            sc_collage_row[:, :, num_channel_cols] = image_cropped
            num_channel_cols += 1

        for c in range(num_channel_cols):
            axarr[index, c].imshow(
                sc_collage_row[:, :, c],
                interpolation=None,
                cmap="gray",
                clim=(0, clim_max[c]),
            )

        if objects:
            c += 1
            if index == 0:
                axarr[index, num_channel_cols].set_title("Cell Object")
            imPath = df.loc[index, "Path_CellObjects"]
            object_im = np.squeeze(read_from_bucket.read_image(imPath))
            # a rest seems to reduce weirdness in fig creation of objects?
            time.sleep(10)
            object_im_crop = crop_single_cell_image(
                object_im, xCenter, yCenter, halfbox_size
            )

            # get value of object_im at midpoint (centered object)
            val = object_im_crop[int(box_size/2),int(box_size/2)]

            # find where in image the value is equal to val
            mask = np.where(object_im_crop == val, True, False)

            # add surrounding object outlines
            outlines = skimage.segmentation.find_boundaries(object_im_crop, mode='inner')
            maskfig = mask + outlines
            # convert boolean mask to float, invert for better look
            maskfig = maskfig.astype(int)^1
            axarr[index, c].imshow(maskfig, cmap='Greys')
            # add label matrix
            c += 1
            axarr[index, c].imshow(object_im_crop)

        if color:
            c+=1
            color_im, colorImagesList = CP_to_RGB_single(
                sc_collage_row[:, :, : len(im_channels)], im_channels, percentile
            )

            axarr[index, c].imshow(color_im, interpolation='nearest')
            if index == 0:
                axarr[index, c].set_title("composite")

            if objects:
                c+=1
                outstack = np.stack([outlines,outlines,outlines], axis=2)
                white_im = np.ones((box_size,box_size,3), dtype=np.uint8)
                white_im = white_im*outstack
                # add outline on top of color image and clip appropriately
                color_im_out = np.clip(color_im/255 + white_im, a_min=0, a_max=1)

                axarr[index, c].imshow(color_im_out, interpolation='nearest')
                if index == 0:
                    axarr[index, c].set_title("composite+outline")

        # remove ticks from all boxes
        for c in range(columns_count):
            axarr[index, c].axes.xaxis.set_ticks([])
            axarr[index, c].axes.yaxis.set_ticks([])

        if "Metadata_Label" in df.columns:
            imylabel = (
                df.loc[index, "Metadata_Label"]
                + "\n"
                + df.loc[index, "Metadata_Foci_Barcode_MatchedTo_Barcode"][0:9]
            )
        else:
            imylabel = df.loc[index, "Metadata_Foci_Barcode_MatchedTo_Barcode"][0:12]
        axarr[index, 0].set_ylabel(imylabel)
    return f