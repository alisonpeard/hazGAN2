import os

# add new relevant variables here
long_names = {
    "GUST_10m": "gust",
    "PRMSL_msl": "apcp",
    "APCP_sfc": "prmsl"
}

def parse_input_pattern(indir, field):
    """
    Parse the input pattern for the given field.
    """
    return os.path.join(indir, long_names[field], '*.nc')


def filter_files(files, year):
    """
    Filter the files for the given year.
    """
    def year_in_file(f):
        fyear = f.split("_")[-1]
        fyear = fyear[:4]
        return str(fyear) == str(year)

    return [f for f in files if year_in_file(f)]

def clip_to_bbox(ds, xmin, xmax, ymin, ymax):
    """
    Clip the dataset to the given bounding box.
    """
    return ds.sel(longitude=slice(xmin, xmax), latitude=slice(ymin, ymax))