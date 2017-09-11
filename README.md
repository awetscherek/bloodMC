blood_MC
===========

This software package provides the necessary MATLAB scripts and data sets to reproduce Fig. 5 in the manuscript 
"Characterization of the diffusion coefficient of blood" by C. Funck, F. B. Laun and A. Wetscherek, accepted for publication in the Journal Magnetic Resonance in Medicine. 

blood_MC includes code to perform Monte-Carlo simulations of water diffusion in blood, modeled as a two-compartment model and code to calculate ADC_0, ADC_b and K_app from the obtained particle trajectories. The results of the simulations for the blood count parameters of the samples in the manuscript are provided. 

---------------------------------------------------------------------------

Contents:

  0. Copyright

  1. skript_bloodMC.m : Monte-Carlo simultation of water diffusion in Blood

  2. skript_calcADC.m : Calculation of ADC and Kurtosis based on the simulated trajectories from skript_bloodMC.m
  
  3. skript_figure5.m : Just plots the results calculated with skript_calcADC.m

---------------------------------------------------------------------------

0. Copyright
------------

This software is covered by the following BSD license:
---------------------------------------------------------------------------

Copyright (c) 2017, Andreas Wetscherek
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in
  the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

---------------------------------------------------------------------------



1.  skript_bloodMC.m : Monte-Carlo simultation of water diffusion in Blood
--------------------------------------------------------------------------

NOTE: by default the number of particle trajectories is only set to 1000, since the completion of the script for the 300000 traces, which was performed for the manuscript, takes a long time. However, to get useful results regarding a T-dependence one should simulate at least 50000 traces.

The skript employs a mex function MCKernel3D_inout.cpp, which for a time interval (specified via number of steps NSteps and diffusion distance dr travelled within each timestep), a starting position (x0, y0, z0) and a geometry (L1, L2, Lr, Lz) calculates center of mass within that time interval and final position of the particle. Different membrane permeability is modelled via the parameters k2D_e and k2D_p corresponding to the products 2*kappa_i/D_i. In addition the time intervals spend within erythrocytes and plasma are saved for each particle.

MCKernel3D_inout is based on a grid of cubic unit cells of dimension L1 x L1 x L2 in which a red blood cell is centered. The red blood cell is modeled as a cylinder of diameter 2*Lr and height 2*Lz.

Settings for the simulations are specified at the beginning of skript_bloodMC.m. HCT, MCV based on blood count and T2 are provided for each of the samples measured in the manuscript (for the latter see Fig. 4 of the manuscript). The experimental settings for diffusion time T and TE used in the manuscript can also be found.

For control purposes a simulation of a fully permeable membrane with D_p and D_e set to the free diffusion coefficient of water at body temperature can also be performed. This subsection of the script should not be executed, if a blood simulation is to be performed, since it overwrites some of the parameters.

The skript will loop over all samples (HCT values in experiment.HCT_), performing the following steps for each of them:

a) Calculate the scaled red blood cell (diameter Ld and height Lh) corresponding to the MCV value in experiment.MCV_ employing the function calcRBCSize.m

b) Calculate the size of the unit cell, such that the sample HCT value corresponds to the volume ratio between RBC and unit cell (calcUnitCell.m)

c) Calculate other parameters, like:
  - number of particles starting in/outside RBC using calcNinRBC.m
  - membrane permeability corresponding to chosen pre-exchange lifetime and surface-to-volume ratio
  - number of steps needed to match the criterion that the step size should be less than dr_rel * the shortest dimension
  
d) Calculation of interval durations such that flow-compensated and bipolar ADCs can be calculated for experiment.TE and the diffusion times in experiment.T_ using calcIntervalDurations.m

e) Monte-Carlo simulations are performed for each particle, looping over the time intervals.

The resulting output files for the settings used in the manuscript are ~2.53 GB in total and therefore not included in this distribution. The author can be contacted in someone is interested in them. Those output files contain for each particle the location at beginning and end of each interval, as well as the center of mass position of the trajectory travelled during each interval. Furthermore for each interval the time of the first membrane transit and the time spent in- and outside a red blood cell is saved.


2. skript_calcADC.m : calculate diffusion properties from MC simulations
------------------------------------------------------------
 
skript loads from the data files generated by the MC-Skript:
- starting position x0,y0,z0 of each particle
- center of mass of trajectory during each interval: xcm, ycm, zcm
- time spent inside (t_in) and outside (t_out) red blood cells
- time of first membrane transit (t_first)

afterwards some sanity checks are performed, whether the geometric HCT corresponds to the HCT calculated from the residence times.

first membrane transit times are used to calculate lifetime in plasma / red blood cell

then the ADC / K_app results are calculated for each diffusion time in experiment.T_:

a) calculating center of mass for each time interval during which the diffusion gradient was constant (2 intervals for bipolar, 4 for the used flow-compensation scheme) and overall displacement weighted with the sign of the diffusion gradient amplitude

b) calculation of gradient amplitude necessary for b=400 s/mm^2

c) displacements are weighted along many different directions to obtain a trace-weighted diffusion coefficient

d) ADC_b (Dr_mp / Dr_fc) is calculated via -ln( | sum_k(e^(i * phi_k) | ), where the accumulated phase phi_k is given by the product of larmor const, gradient amplitude, diffusion time and the weighted displacement

e) ADC_0 (Da_mp / Da_fc) is calculated as <phi_k^2>/2

f) K_app is calculated based on the 4th moment of the distribution of phases phi_k.

The results are exported again into .mat files. Those .mat files correspond to the ones provided in the 300k subfolder. 



 
 