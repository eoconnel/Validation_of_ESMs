# Validation_of_ESMs
Code repository to accompany master's thesis: Validation of Earth System Models

## Author:  
Ethan O'Connell  
et214214@dal.ca  
Dalhousie University
Master of Science in Statistics  

## Supervisor:  
Michael Dowd 
michael.dowd@dal.ca 
Dalhousie University, department of Mathematics and Statistics  

## Abstract


## Contents
Most of the experimentation was conducted using jupyter notebooks to facilitate visualization and testing. The contents of each of the notebooks and python modules are explained below:  
* __network.py__: neural network class as explained in `modular_neural_network.ipynb`.  
* __network_utils.py__: utility functions to facilitate neural network experimentation.  
* __modular_neural_network.ipynb__: explanation of the `network.py` module.  
* __hyperparameter_tuning.ipynb__: a set of experiments conducted to determine an approximate learning rate and number of epochs for the binary classification test problems.  
* __binary_classification.ipynb__: an explanation of classification problems and training a series of networks over binary classification regions generated using hand-drawn black and white images.  
* __concave_classification.ipynb__: an extension of binary classification problems to more complex regions involving concavity.  
* __function_approximation.ipynb__: an introduction to function approximation problems and a demonstration of one and two dimensional cases.  
* __convolutional_networks.ipynb__: an explanation of convolutional networks implemented over the MNIST hand-written digit classification problem as well as an attempted fully convolutional network.  
* __universal_approximation.ipynb__: a constructive proof and visual demonstration of the Universal Approximation Theorem for neural networks.
