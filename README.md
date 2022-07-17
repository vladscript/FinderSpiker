
# FinderSpiker

## Calcium Imaging Signal Processing Toolbox
![Logo](/figures/Logo_FinderSpiker.png)

### Experiment Profile:
* For synthetic dyes or genetically encoded indicators of calcium (neuronal) activity
* Calcium Transients Offline Detector
  * Automatic mode and Visual Inspection Mode.
* Multiple conditions (e.g. CTRL, +DRUG A, +DRUG B, etc).
* It can merge 2 different dyes (eg GCaMP/Syn + TdTomtato/cre) in the same field

# See [**Quick User Guide**](http://htmlpreview.github.io/?https://github.com/vladscript/FinderSpiker/blob/master/html/USER_GUIDE.html)

### Brief Description
Optimized for Acquired Data in [ImPatch](http://impatch.ifc.unam.mx/) and CSV exported data from [Igor WaveMetrics](https://www.wavemetrics.com/downloads)
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

### Algorithms Flow Charts
  - Signal Processing
  - Neuronal Ensembles
  - Functional Network
  - Machine Learning

### Neuronal Ensembles Algorithm demo (toy example):
Open FinderSpiker as Current Folder in MATLAB and run:

```

>>Import_FinderSpiker

>>cd Demo' Scripts'\

>>NeuronalEnsemblesToyExample % see code

%>>edit NeuronalEnsemblesToyExample;

```

### Third Party Software:
  - [SpaRSA](https://www.lx.it.pt/~mtf/SpaRSA/)
  - [Rain Cloud Plots](https://github.com/RainCloudPlots/RainCloudPlots)
  - [ReadImageJROI](https://www.github.com/DylanMuir/ReadImageJROI)
  - [Neuronal Ensembles] (https://github.com/PerezOrtegaJ/Neural_Ensemble_Analysis)

### References, cites:

Version 1 (2022)

  - [Striatal Neuronal Ensembles Reveal Differential Actions of Amantadine and Clozapine to Ameliorate Mice L-DOPA-Induced Dyskinesia, Neuroscience, 2022](https://authors.elsevier.com/a/1fO~U15hTtrK7A)
  - [Dopamine D2 and adenosine A2A receptors interaction on Ca2+ current modulation in a rodent model of Parkinsonism, ASN Neuro, 2022](https://doi.org/10.1177%2F17590914221102075)
  - [Direccionamiento de axones de neuronas dopaminérgicas humanas por semaforina 3C en cultivos organotípicos de cerebro, Tesis de Maestría en Ciencias Bioquímicas UNAM, 2020](https://tesiunam.dgb.unam.mx/F/7C2N4TM4UBA5C4BJURKYUXXXAEM71I5PGHKBURIHVGC248I2TE-26944?func=full-set-set&set_number=511707&set_entry=000001&format=999)
  
  ### Use this reference to cite
  
Vladimir M. Calderón, Aldo Luna-Leal, Alejandra Gómez-Paz, Fernanda Ramírez-López, Mario Arias-García, Esther Lara-González, Elvira Galarraga, José Bargas,
Striatal Neuronal Ensembles Reveal Differential Actions of Amantadine and Clozapine to Ameliorate Mice L-DOPA-Induced Dyskinesia,
Neuroscience,
Volume 492,
2022,
Pages 92-107,
ISSN 0306-4522,
https://doi.org/10.1016/j.neuroscience.2022.03.036.
(https://www.sciencedirect.com/science/article/pii/S0306452222001610)

Abstract: Amantadine and clozapine have proved to reduce abnormal involuntary movements (AIMs) in preclinical and clinical studies of L-DOPA-Induced Dyskinesias (LID). Even though both drugs decrease AIMs, they may have different action mechanisms by using different receptors and signaling profiles. Here we asked whether there are differences in how they modulate neuronal activity of multiple striatal neurons within the striatal microcircuit at histological level during the dose-peak of L-DOPA in ex-vivo brain slices obtained from dyskinetic mice. To answer this question, we used calcium imaging to record the activity of dozens of neurons of the dorsolateral striatum before and after drugs administration in vitro. We also developed an analysis framework to extract encoding insights from calcium imaging data by quantifying neuronal activity, identifying neuronal ensembles by linking neurons that coactivate using hierarchical cluster analysis and extracting network parameters using Graph Theory. The results show that while both drugs reduce LIDs scores behaviorally in a similar way, they have several different and specific actions on modulating the dyskinetic striatal microcircuit. The extracted features were highly accurate in separating amantadine and clozapine effects by means of principal components analysis (PCA) and support vector machine (SVM) algorithms. These results predict possible synergistic actions of amantadine and clozapine on the dyskinetic striatal microcircuit establishing a framework for a bioassay to test novel antidyskinetic drugs or treatments in vitro.
Keywords: L-DOPA-Induced Dyskinesia; dyskinetic striatal microcircuit; calcium imaging; neuronal ensembles; identification of antidyskinetic drugs actions
