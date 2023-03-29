import pybedtools


def concatenate_BED_files(bed1, bed2, sort=True):
    if sort:
        bed = pybedtools.BedTool(str(bed1) + str(bed2), from_string=True).sort()
    else:
        bed = pybedtools.BedTool(str(bed1) + str(bed2), from_string=True)
    return bed
