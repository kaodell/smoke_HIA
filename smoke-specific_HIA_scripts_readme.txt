smoke-specific_HIA_scripts_readme.txt
written by Katelyn O'Dell
04.27.21

This folder contains all the python scripts used in the analysis for
"Estimated mortality and morbidity attributable to smoke plumes in the US: Not just a western US problem"

Authors: Katelyn O’Dell, Kelsey Bilsback, Bonne Ford, Sheena E. Martenies, Sheryl Magzamen, 
Emily V. Fischer, and Jeffrey R. Pierce

For submission to GeoHealth special issue, “Fires in the Earth System”
Submission Open 08 April 2020 - 31 May 2021

In order to perform the analysis presented in the paper, 
these codes are run in the following order:

Stage 1: Prep Input Datasets
(1.1) regrid_population.py
	This code was adapted from a regrid.py code from Dr. Kelsey Bilsback, and is run in python 2.
	It requires the kriging grid which is available for download at: https://doi.org/10.25675/10217/230602
	and the SEDAC gridded population of the world v4, available at: https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11

(1.2) mk_state_grid.py
	Python script to assign grid cells to states. Was done by hand, but may be faster to
	just grab a higher resolution shapefile in the future.
	Requires population output from (1.1) to check state assignments

(1.3) prep_smokePM_4HIA.py
	This code loads all years of the PM data, calculates smoke PM, and saves as a netCDF 
	of the combined total PM, smoke PM, and HMS to use in the following scripts.
	
	This code creates Figures 1 and S2 in the manuscript.
	
	Code requires: 
		krigged PM2.5 dataset available for download at: https://doi.org/10.25675/10217/230602
		a state mask which assigns grid cells to each state, created in step 1.2
		re-gridded SEDAC population for 2010, created in step 1.1

Stage 2: Acute Smoke PM2.5 HIA
(2.1) pool_asth_OR.py
	This code conducts a monte-carlo fixed-effect pooling of odds ratios for asthma ED
	visits and asthma hospitalizations from the US studies included in Borchers Arrigada et al.
	(2019) with additional recently published western US studies and eastern US studies. 
	This code creates Figure S1. 
	
	Code requires:
		csv file with odds ratios to pool. All of these values are available in the
		publications referenced in Figure S1.

(2.2) calc_acute_HIA_pool_EPAregions.py
	This code calculates the acute smoke PM2.5 HIA and sums by regions in the paper.
	Figures 2, 3, S3, S4, and S5 are created in this code.
	
	Code requires:
		Krigged PM output from step 1.3, population output from step 1.1, state grid
		from step 1.2, and acute baseline rates available via HCUP (details in paper).
		
(2.3) calc_acute_HIA_pool_states.py
	This code calculates the acute HIA by state and creates Figures S5-S8.
	
	Code requires:
		HIA output from step 2.2, population output from step 1.1, state grid
		from step 1.2, and acute baseline rates available via HCUP (details in paper).

Stage 3: Chronic Smoke PM2.5 and Smoke HAPs HIA
(3.1) calc_mortalities_DALYs_gemm.py
	This code calculates total PM and smoke PM mortalities and smoke PM DALYs. 
	It is modified from code written by Dr.s Kelsey Bilsback and Jack Kodros. 
	
	Code requires:
		Krigged PM output from step 1.3, population output from step 1.1, and the state 
		grid from step 1.2.
		
(3.2) plot_mortalities_by_country_gemm.py
	This code creates Figure 4 in the main text. It is modified from code written by 
	Dr.s Kelsey Bilsback and Jack Kodros. 
	
	Code requires:
		Mortalities output from step 3.1, population output from step 1.1, and the state
		grid from step 1.2
	
(3.3) calc_DALYs_fromHAPs.py
	This code estimates HAPs concentrations, calculates DALYs from HAPs,
	and creates Figure 5 in the main text and Figures S10 and S11 in the supplement.
	
	Code requires:
		HAPs to PM ratios from O'Dell et al. (2020) supplementary data, 
		combined damage and effect factors from Huijbregts et al. (2005) appendix 1,
		krigged PM2.5 output from step 1.3, population file from step 1.1, state grid 
		from step 1.2, and PM DALYs from step 3.1
		
and that's a wrap!