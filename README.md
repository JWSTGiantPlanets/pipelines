# JWST reduction pipelines
The reduction and visualisation code in this repository is used to call and extend the [standard JWST reduction pipeline](https://github.com/spacetelescope/jwst/) to process, reduce and analyse data for solar system observations. 

> [Setup](#setup) | [MIRI Pipeline](#miri-mrs-pipeline) | [NIRSPEC Pipeline](#nirspec-pipeline) | [Reference](#reference)

## Setup
Requirements:
- Python 3.10 or above
- ~10GB of disk space (for the CRDS cache files downloaded when running the JWST data reduction)

To run the pipelines yourself, first create a local clone of this repository with:
```
git clone https://github.com/JWSTGiantPlanets/pipelines.git
```

Then install the required Python modules with:

```
cd pipelines
pip install -r requirements.txt
```

## MIRI MRS pipeline

TODO 


### Flat fields
Synthetic flat fields derived from the observations of Saturn in November 2022 can be downloaded from the [supplementary material for Fletcher et al. (2023)](https://github.com/JWSTGiantPlanets/saturn-atmosphere-miri). The Saturn observations have poor SNR at some wavelengths (particularly in channel 1B) which reduces the quality of the flats at some wavelengths. This means that data corrected with the Saturn flats will have increased noise at the wavelengths where there was poor SNR in the original Saturn observations.

Our testing suggests that the flat field appears to be slightly different for different observations, potentially due to non-linear variations in the flat with brightness or [time variability of the flat](https://blogs.nasa.gov/webb/2023/04/21/mid-infrared-instrument-operations-update-2/). If the Saturn flats described above do not adequately remove striping from your data, you may wish to generate synthetic flats from your dataset instead.

You can create flat fields yourself using the [`construct_flat_field`](https://github.com/JWSTGiantPlanets/pipelines/blob/main/construct_flat_field.py) script to create a set fo flat fields directly from your dataset. The script takes a set of `stage3` or `stage3_desaturated` dithered observations and uses these to construct a flat field which minimises variation in brightness of each location on the target between dithers. This was designed to work for extended source observations which fill the MIRI field of view and have at least 4 dithers. Observations of objects which do not fill the FOV or have fewer dithers are unlikely to produce reliable flat fields.

![MIRI pipeline summary figure](images/pipeline_summary.png)
_Figure from Fletcher et al. (2023) showing the custom MIRI pipeline for Saturn. The first column shows the output of the standard JWST pipeline, which still contains significant flat field effects (a & d), saturation (g), and partial saturation (dark pixels in g & j). The second column shows the data after the desaturation step is applied, and the third column shows the data after the flat field correction is applied._

## NIRSPEC pipeline
_The NIRSPEC pipeline is in development and will be released soon..._


## Reference
Fletcher et al. (2023). _Saturn's Atmosphere in Northern Summer Revealed by JWST/MIRI_. Manuscript in preparation.