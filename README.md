# FinderSpiker

## Calcium Imaging Signal Processing Toolbox

### Experiment Profile:
* For synthetic dyes or viral flurophores
* Calcium Transients Offline Detector
  * Automatic mode and Visual Inspection Mode.
* Multiple conditions (e.g. CTRL, +DRUG A, +DRUG B, etc).
* It can merge 2 different dyes (eg GCaMP/Syn + TdTomtato/cre) in the same field

# See [**Quick User Guide**](http://htmlpreview.github.io/?https://github.com/vladscript/FinderSpiker/blob/master/html/USER_GUIDE.html)

### Brief Description
Optimized for Acquired Data in [ImPatch](http://impatch.ifc.unam.mx/)
* **Input**
  - List of videos (e.g. CTRL01,.., CTRLN,DRUGA01,...,DRUGA0N,...)
  - List of ROIs Coordianates (CSV file from ImPatch)
* **Output**
  - Activity Matrix: **R**: Cells x Frames Activity as 1's
  - It identifies neural ensembles as coactive neurons
  - It creates functional network for [Gephi](https://gephi.org/) Visualization by wiring neurons that fired together
* **Algorithms**
  - Calcium Transients Detection:
    - Detrending (RLAOSS)
    - Denoising (Wavelet Analysis)
    - Sparse Deconvolution (LASSO regularization)
  - Neural Ensembles Detection:
    - Hierarchichal Clustering
    - Cross Validation by Naive Bayes Classificator
* **Results Visualization**
  - Boxplots and Datasets for Diverse Experiments
  - Cummulative Distribution Functions of Descriptive Features:
    - Activity Features
    - Neural Ensembles Features
    - Functional Network Features (from Gephi)
  - Rebuild denoised video and Ca Transients Detection

* Third Party Software:
  - [SpaRSA](https://www.lx.it.pt/~mtf/SpaRSA/)
  - [Rain Cloud Plots](https://github.com/RainCloudPlots/RainCloudPlots)
  - 
