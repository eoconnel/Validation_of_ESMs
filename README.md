# Validation_of_ESMs
Code repository to accompany master's thesis: Improved Validation of Earth System Models

## Author:  
Ethan O'Connell  
et214214@dal.ca  
Dalhousie University  
Master of Science in Statistics  

## Supervisor:  
Michael Dowd  
michael.dowd@dal.ca  
Dalhousie University  
Department of Mathematics and Statistics  

## Abstract
Earth system models are complex numerical models that simulate several aspects of the Earth system. They use coupled systems of partial differential equations to model the physical dynamics of the atmosphere and the oceans, as well as other chemical and biological processes. Due to their nonlinearity and complex geometry, they are solved numerically using a variety of computational approaches. With increases in computational power and model complexity, as well as the increasing availability of calibration datasets, there is an increased need for the assessment of numerical model performance to better inform decision making. In this thesis, we provide an overview of common statistical methods and develop new approaches for the validation of numerical models applicable to Earth system models, with a particular focus on the definition of application-specific validation metrics and diagnostics.

As a concrete focus, we apply these methods to the simulation of nearshore ocean temperatures along the Atlantic coast of Nova Scotia. This coastal region is important for marine seagrass ecosystems and provides features of statistical interest since it exhibits strong short-term variability and extreme temperatures. For this reason, we emphasize techniques for the evaluation of short-term variability using frequency-based time series decomposition. We also present various approaches suitable for the evaluation of extreme temperatures, based on characterizing marine heatwaves using extreme value theory. To evaluate the effectiveness of these approaches across various types of numerical models, we compare two iterations of a regional FVCOM model, the GLORYS12v1 reanalysis dataset, and a representative CMIP6 model against temperature observations at two nearshore study sites. Our findings highlight the importance of developing tailored validation procedures, particularly in regions with substantial variability and extremes.
