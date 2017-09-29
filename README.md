# stochastic_design

## Brief description
These codes implement the optmization algorithms and simulations presented in "Optimization-basd synthesis of stochastic biocircuits with statistical specifications" by Yuta Sakurai and Yutaka Hori.

These program codes are contributed by 
 - Yuta Sakurai (Keio University)
 - Yutaka Hori (Keio University)

## Prerequisites
- MATLAB 2016b 
- SeDuMi 1.3 for the semi-algebraic optimizations. Download SeDuMi 1.3 at https://github.com/sqlp/sedumi

## How to run the codes

1. Open MATLAB and add path to the "Add_Path" folder. 

2. Follow the instructions below to produce the figures.

### Fig. 2
- Run "ssaTrscTrsl.m" for the SSA simulations (Fig. 2(B)).
- Original data are saved in "mRNAvsProteinStatistics.mat".

### Fig. 3
- Run "ssa_RepressorFeedback.m" for the SSA simulations (Fig. 3(B)).
- Original data are saved in "ssa_RepressorFeedback.mat".

- Run "SteadyStateAnalysis_OrdervsStatistics.m" to compute the mean and the CV in Figs. 3(C) and (D).
- Original data are saved in "SteadyStateAnalysis_OrdervsStatistics.mat".

### Fig. 4
- Run "Fig_B/FeasibleRegion_FigureOutput.m" to load results and plot the optimization results in Fig. 4(B).
- To run this code, you first need to run the optimization codes in "Fig_C_D" and "Fig_E" folders (see below for the details). Then, move the output files of the optimization codes (.mat files) to the "Fig_B" folder. 

- To plot Fig. 4(C)-(D), first, run "Fig_C_D/SteadyStateAnalysis_TrslDeg_Repressor.m" to compute the mean and the CV of the repressor using the optimization program.
- Then, run "Fig_C_D/OutputALLFigures.m" to generate the plot.

- To plot Fig. 4(E), first, run "Fig_E/SteadyStateAnalysis_TrslDeg_Reporter_Ave.m", and then, run "Fig_E/SteadyStateAnalysis_TrslDeg_Reporter_CV.m" to compute the CV of the reporter.
- Finally, run "Fig_E/OutputALLFigures.m" to generate the plot.

### Fig. 5
- For Fig. 5(A), run "SteadyStateAnalysis_AvevsCV_Repressor_K1.m", "SteadyStateAnalysis_AvevsCV_Repressor_K2.m" and ""SteadyStateAnalysis_AvevsCV_Repressor_K3.m" in the "Codes" folder to compute the CV of the repressor.
- Then run "outputfigures.m" to generate the plot.

- For Fig. 5(B), run "SteadyStateAnalysis_TrslDeg_mRNA1.m" first.
- Then run "trsl_mRNA_AveCV.m" to generate the plot.

- For Fig. 5(C), run "SteadyStateAnalysis_RepressorReporter_Relation.m" first.
- Then run "RepressorReporter_Relation_FigureOutput.m" to generate the plot.


### Fig. 6
- Run "SteadyStateAnalysis_SensitivityAnalysis.m" to compute the sensitivity. 
- Run "SensitivityFigureOutput.m" to plot the result (Fig. 6).

### Fig. S.1
Run "SteadyStateAnalysis_TrscTrsl_Ellipsoid.m" to plot Fig. S1.


# Appendix
The files in the "Add_Path" folder include common library codes that are referred from other codes.
The contents of "Add Path" folder are as follows.
- Analysis_Output : Codes for solving SDP with SeDuMi in the "analysis mode"
- Common_Files : Codes to transform the optimization problems into the SDP form
- SensitivityAnalysis_Output : Codes to solve SDP with SeDuMi in the "sensitivity analysis mode"



