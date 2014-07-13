CellProfiler Pipeline: http://www.cellprofiler.org
Version:1
SVNRevision:11710

LoadData:[module_num:1|svn_version:\'Unknown\'|variable_revision_number:6|show_window:False|notes:\x5B\x5D]
    Input data file location:Elsewhere...\x7CC\x3A\\Users\\Sissela\\Documents\\MATLAB
    Name of the file:Only_Classification_Slides.csv
    Load images based on this data?:Yes
    Base image location:Elsewhere...\x7CC\x3A\\Users\\Sissela\\Documents\\MATLAB\\08-tophats-registered\\slideA
    Process just a range of rows?:No
    Rows to process:1,100000
    Group images by metadata?:No
    Select metadata fields for grouping:position
    Rescale intensities?:Yes

IdentifyPrimaryObjects:[module_num:2|svn_version:\'10826\'|variable_revision_number:8|show_window:False|notes:\x5B\x5D]
    Select the input image:DO
    Name the primary objects to be identified:blobs
    Typical diameter of objects, in pixel units (Min,Max):2,10
    Discard objects outside the diameter range?:Yes
    Try to merge too small objects with nearby larger objects?:Yes
    Discard objects touching the border of the image?:Yes
    Select the thresholding method:Manual
    Threshold correction factor:1
    Lower and upper bounds on threshold:0.000000,1.000000
    Approximate fraction of image covered by objects?:0.01
    Method to distinguish clumped objects:Intensity
    Method to draw dividing lines between clumped objects:Intensity
    Size of smoothing filter:1
    Suppress local maxima that are closer than this minimum allowed distance:3
    Speed up by using lower-resolution image to find local maxima?:No
    Name the outline image:BlobOutlines
    Fill holes in identified objects?:Yes
    Automatically calculate size of smoothing filter?:No
    Automatically calculate minimum allowed distance between local maxima?:No
    Manual threshold:0.05
    Select binary image:None
    Retain outlines of the identified objects?:Yes
    Automatically calculate the threshold using the Otsu method?:Yes
    Enter Laplacian of Gaussian threshold:0.5
    Two-class or three-class thresholding?:Two classes
    Minimize the weighted variance or the entropy?:Weighted variance
    Assign pixels in the middle intensity class to the foreground or the background?:Foreground
    Automatically calculate the size of objects for the Laplacian of Gaussian filter?:Yes
    Enter LoG filter diameter:5
    Handling of objects if excessive number of objects identified:Continue
    Maximum number of objects:500
    Select the measurement to threshold with:None

MeasureObjectIntensity:[module_num:3|svn_version:\'10816\'|variable_revision_number:3|show_window:False|notes:\x5B\x5D]
    Hidden:17
    Select an image to measure:1_1
    Select an image to measure:1_2
    Select an image to measure:1_3
    Select an image to measure:1_4
    Select an image to measure:2_1
    Select an image to measure:2_2
    Select an image to measure:2_3
    Select an image to measure:2_4
    Select an image to measure:3_1
    Select an image to measure:3_2
    Select an image to measure:3_3
    Select an image to measure:3_4
    Select an image to measure:4_1
    Select an image to measure:4_2
    Select an image to measure:4_3
    Select an image to measure:4_4
    Select an image to measure:DO
    Select objects to measure:blobs

ExportToSpreadsheet:[module_num:4|svn_version:\'10880\'|variable_revision_number:7|show_window:False|notes:\x5B\x5D]
    Select or enter the column delimiter:Comma (",")
    Prepend the output file name to the data file names?:No
    Add image metadata columns to your object data file?:No
    Limit output to a size that is allowed in Excel?:No
    Select the columns of measurements to export?:Yes
    Calculate the per-image mean values for object measurements?:No
    Calculate the per-image median values for object measurements?:No
    Calculate the per-image standard deviation values for object measurements?:No
    Output file location:Default Output Folder\x7CC\x3A\\\\Users\\\\User\\\\Documents\\\\MATLAB\\\\cells\\\\figs
    Create a GenePattern GCT file?:No
    Select source of sample row name:Metadata
    Select the image to use as the identifier:None
    Select the metadata to use as the identifier:None
    Export all measurements, using default file names?:Yes
    Press button to select measurements to export:blobs\x7CIntensity_MaxIntensity_DO,blobs\x7CIntensity_MaxIntensity_1_1,blobs\x7CIntensity_MaxIntensity_1_3,blobs\x7CIntensity_MaxIntensity_1_2,blobs\x7CIntensity_MaxIntensity_1_4,blobs\x7CIntensity_MaxIntensity_3_1,blobs\x7CIntensity_MaxIntensity_3_3,blobs\x7CIntensity_MaxIntensity_3_2,blobs\x7CIntensity_MaxIntensity_3_4,blobs\x7CIntensity_MaxIntensity_2_1,blobs\x7CIntensity_MaxIntensity_2_3,blobs\x7CIntensity_MaxIntensity_2_2,blobs\x7CIntensity_MaxIntensity_2_4,blobs\x7CIntensity_MaxIntensity_4_1,blobs\x7CIntensity_MaxIntensity_4_3,blobs\x7CIntensity_MaxIntensity_4_2,blobs\x7CIntensity_MaxIntensity_4_4,blobs\x7CLocation_Center_Y,blobs\x7CLocation_Center_X
    Data to export:blobs
    Combine these object measurements with those of the previous object?:No
    File name:blobs.csv
    Use the object name for the file name?:No
