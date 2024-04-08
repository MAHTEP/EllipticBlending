# EllipticBlending
This project aims at implementing the k-epsilon Lag Elliptic Blending turbulent model for OpenFOAM. 

## Software requirements
The model has been developed and tested in [OpenFOAM v2212](https://www.openfoam.com/download/release-history#v2212).


## Usage

### Model source code
In the folder ```TurbulenceModels/turbulenceModels/RAS/kEpsilonLagEB``` you can find the implementation of the model in OpenFOAM.

### Test case
In the folder ```BFS``` you can find the BFS case adapted for the lagEB model. You just need to run the Allrun script.

## References
The elliptic model was originally implemented in the commercial software STAR-CCM+, and it was developed as described in the following two publications:

[1] [Development of an elliptic-blending lag model for industrial applications](https://www.researchgate.net/publication/314229391_development_of_an_elliptic-blending_lag_model_for_industrial_applications), S. Lardeau, F. Billard, 54th AIAA Aerospace Sciences Meeting (2016), DOI: 10.2514/6.2016-1600
      
[2] [AN ELLIPTIC BLENDING LAG MODEL FOR FLOWS IN THERMAL-HYDRAULICS SYSTEMS](https://www.researchgate.net/publication/308413285_AN_ELLIPTIC_BLENDING_LAG_MODEL_FOR_FLOWS_IN_THERMAL-HYDRAULICS_SYSTEMS), R. Tunstall, S. Lardeau et al. (2016)


    
