### Notes

Use `submit-ouce.sh` to run data processing jobs on the SoGE cluster. Due to GPU issues we need to use ARC for training StyleGAN. Use the `submit-arc.sh` script for that.

#### Event extraction

Currently using the storm tracks supplied by Colin, which are extract based on distance to each DNO license region. I have downloaded the [ECMWF](https://cds.climate.copernicus.eu/datasets/sis-european-wind-storm-reanalysis?tab=download) tracks and they are stored locally in `/Users/alison/Documents/dphil/data/ecmwf/storm-tracks`.

For alignment with Colin's model, I will use the original tracks for now.


#### Marginal distributions

Using all GPD fits with Murphy (2025) threshold selection instead of Bader method. Would be good to include uncertainty but that's future work.