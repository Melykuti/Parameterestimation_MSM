# Parameterestimation_MSM
Parameter estimation of seeding and contamination rates in a subcritical percolation model with colouring on a triangular lattice, using the method of simulated moments (MSM). Code implementation in MATLAB (The MathWorks, Inc.).

#### Reference:
Felix Beck, Bence Mélykúti.  
_Parameter estimation in a subcritical percolation model with colouring_  
arXiv:1604.08908  
[https://arxiv.org/abs/1604.08908](https://arxiv.org/abs/1604.08908).

#### Note:
An additional MATLAB file must be downloaded!

fminsearchbnd:
* **fminsearchbnd** behaves similarly to fminsearch, except you can add bound constraints.
* source and description: [http://de.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon/content/FMINSEARCHBND/fminsearchbnd.m](http://de.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon/content/FMINSEARCHBND/fminsearchbnd.m)
* author: John D'Errico

#### Authors:
Felix Beck, Bence Mélykúti (University of Freiburg, Germany).  
2015-2017

#### Contents:
* **msm.m** is the main file to be run for parameter estimation, with comments about variables and input requirements.
* **structure.txt** describes the interrelation of program files.
* **Synthetic\_datasets\_for\_estimation** folder contains our standard synthetic datasets that can be used to test parameter estimation
* **Data\_estimates** contains results of parameter estimation. These results are reported in Section 6 of the paper.
* **Mathematica** has a Mathematica (Wolfram Research, Inc.) notebook for calculations in Section 7 of the paper.
* **Densityestimation** 
