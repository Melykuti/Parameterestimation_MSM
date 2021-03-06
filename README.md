# Parameterestimation_MSM
Parameter estimation of seeding and contamination rates in a subcritical percolation model with colouring on a triangular lattice, using the method of simulated moments (MSM). Code implementation in MATLAB (The MathWorks, Inc.).

#### Authors
Felix Beck, Bence Mélykúti (University of Freiburg, Germany).  
2015-2017  

With questions or comments, please contact [Bence Mélykúti](https://github.com/Melykuti).

#### Reference
Felix Beck, Bence Mélykúti.  
_Parameter estimation in a subcritical percolation model with colouring_,  
Stochastics: An International Journal of Probability and Stochastic Processes, **91**(5), 657–694, 2019, [doi:10.1080/17442508.2018.1539089](https://doi.org/10.1080/17442508.2018.1539089) (open access),  
or [arXiv:1604.08908](https://arxiv.org/abs/1604.08908)

#### Dependencies
An additional MATLAB file must be downloaded!

fminsearchbnd:

* **fminsearchbnd.m** behaves similarly to fminsearch, except you can add bound constraints.
* source and description: [http://de.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon/content/FMINSEARCHBND/fminsearchbnd.m](http://de.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd--fminsearchcon/content/FMINSEARCHBND/fminsearchbnd.m)
* author: John D'Errico

#### Contents
* **msm.m** is the main file to be run for parameter estimation. It contains comments about variables and input requirements.
* **createsynthdata.m** is a tool to generate and save random realisations of the percolation process with colouring.
* **plot_figure.m** is a versatile plotting program to display lattices with our percolation process from saved data.
* **latex_table.m**, **latex_table5.m** and **latex_table6.m** created Tables 1-4, Table 5 and Table 6 of the paper in LaTeX format.
* **structure.txt** describes the interrelation of program files.
* The folder **Synthetic\_datasets\_for\_estimation** contains our standard synthetic datasets with known parameter values that can be used to test parameter estimation.
* The folder **Synthetic\_datasets\_for\_plotting** contains our synthetic datasets that were used to create Figures 1 and 3 of the paper.
* The folder **Data\_estimates** contains results of parameter estimation. These results are reported in Section 6 of the paper.
* The folder **Mathematica** has a Mathematica (Wolfram Research, Inc.) notebook for calculations in Section 7 of the paper.
* The folder **Densityestimation** contains:
    - **densityestimation.m**, which estimates the empirical means of the densities of coloured vertices and of adjacent pairs of coloured vertices as a function of lambda (seeding rate) and mu (percolation parameter), for a mesh of lambda and mu pairs, when the number of colours is 1.
    - **data_visualise.m** displays the result of densityestimation.m saved in **dataset_A.mat** and **dataset_B.mat** in several plots.
    - **data_visualise_publication.m** creates the figures that are presented in the Appendix of the paper (arXiv:1604.08908v3).

#### License
The software in the repository _Melykuti/Parameterestimation\_MSM_ is licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License](http://rawgit.com/Melykuti/Parameterestimation_MSM/master/LICENSE.html). For commercial licensing, please contact [Bence Mélykúti](https://github.com/Melykuti).

Copyright (c) 2017, Felix Beck, Bence Mélykúti 
