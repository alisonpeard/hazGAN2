import os

# add new relevant variables here
long_names = {
    "GUST_10m": "gust",
    "PRMSL_msl": "apcp",
    "APCP_sfc": "prmsl"
}

def get_input_file_pattern(indir, field):
    """Parse the input pattern for the given field."""
    return os.path.join(indir, long_names[field], '*.nc')


def filter_files(files, year, antecendent_buffer_days=None):
    """Filter the files for the given year."""
    def year_in_file(f, year):
        fyear = f.split("_")[-1]
        fyear = fyear[:4]
        return str(fyear) == str(year)
    
    files_year = [f for f in files if year_in_file(f, year)]

    if antecendent_buffer_days:
        nyears = 1 + (antecendent_buffer_days // 365)
        for dyear in range(1, nyears + 1):
            antecedent_files = [f for f in files if year_in_file(f, year-dyear)]
            files_year.extend(antecedent_files)
    
    return files_year


def clip_to_bbox(ds, xmin, xmax, ymin, ymax):
    """Clip the dataset to the given bounding box."""
    return ds.sel(longitude=slice(xmin, xmax), latitude=slice(ymin, ymax))


def preprocess(ds):
    """Extra handling to preprocess the dataset."""
    return ds