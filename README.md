# MPhys

This contains the relevant code for my MPhys project "Detection of transiting exoplanets with the TESS space mission". 


## Getting started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites
This code requires the following packages. All available via pip

```
pip install pandas
pip install numpy
pip install sklearn
pip install wotan
pip install scipy
pip install lmfit
```

As well as the Probabilistic Random Forest package available at https://github.com/ireis/PRF


## Feature extraction

To train the (Probabilistic) Random Forest and make predictions with the trained classifier, one must first undertake a stage of feature extraction. In principle this can be done simply by iterating over the light curves but the independent nature of the analysis lends itself well to parallelisation. Parallelisation code is included in this repo.

The feature extraction for the training set requires:
* The output from PHT analysis for each sector
* The output from the manual vetting stage, including stellar radius
* (If TOIs are being included in the training set) The list of TOIs from https://tev.mit.edu/data/

The feature extraction for unvetted light curves requires:
* The output from PHT for each sector
* A .csv file with the same format as that produced by the manual vetting stage. The entries for the labels are ignored if they exist.

Using the `write_csv.py` file, run the
```
write_CSV(sec_list, with_label = True)
```
function to format the input properly for feature extraction. `sec_list` is the list of sectors from which one desires to analyse light curves. `with_label` is a boolean value to notify the code if the data has labels (and so is to be used as a training set) or not (and so the classifier is to be used to predict labels).

At this point, runnning `MPI_extract.py` or `MPI_extract_no_labels.py` will run the feature extraction algorithm for the labelled and un-labelled light curves, respectively. This allows for parallelised operation.

After running both of these files, the `predict.py` code will predict whether to keep or reject each light curve.

Further details on the background and feature analysis can be found in the file `Report.pdf`.

I hope you find this code useful.

Anonymous MPhys student.


