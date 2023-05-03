# JWST reduction pipelines
The reduction and visualisation code in this repository is used to call and extend the [standard JWST reduction pipeline](https://github.com/spacetelescope/jwst/) to process, reduce and analyse data for solar system observations.


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

TODO mention getting flat fields

## NIRSPEC pipeline
_The NIRSPEC pipeline is in development and will be released soon..._


## Reference
> Fletcher et al. (2023). _Saturn's Atmosphere in Northern Summer Revealed by JWST/MIRI_. Manuscript in preparation.