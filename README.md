# Detecting changes in generation and serial intervals under varying pathogen biology, contact patterns and outbreak response

This repository contains code and data to simulate pairwise and cluster-level transmission pairs and obtain the corresponding generation (GI) and serial interval (SI). We compare the power to detect differences in the generation interval of a reference and an alternative pathogen under varying biological characteristics, contact patterns, outbreak response and epidemic dynamics. Overall, the study aims to inform the required sample size under certain outbreak conditions to ensure that future GI and SI studies are well powered.

## Quick start guide
### Data for time varying contact pattern
Original high resolution contact data from a [community in the UK] (https://doi.org/10.1101/479154) can be found in:

> data/raw/haslemere_edge.csv<br/>
> data/raw/haslemere_time.csv

The data was processed and stratified into household and non-household contacts. As the infection process for most diseases typically occurred on time scales lasting more than three days, we extended the contact between pairs of individuals by randomly sampling their daily contacts over weekdays based on the recorded contact patterns on Thursday and Friday and fixed all weekend contacts based on Saturday. The code and processed data can be found in 

> `data/processed`<br/>
> `code/step1a_haslemere_edgelist.R`<br/>
> `code/step1b_haslemere_edgelist_extend.R`

### Setting up parameters and functions for simulation
Data loading and model run script for pairwise transmission is in `code/step2c_inf_pair_simulate.r`. Calls the following R files:
> `code/step2a_inf_pair_function.r` - simulates transmission based on incubation, infectiousness, delay from onset-to-isolation.<br/>
> `code/step2b_inf_pair_param.r` - set up parameters for pairwise transmission in different diseases and epidemic dynamics

#### For pairwise transmission:
* Set `param = param.pair()` or `param = param.pair.disease`
* Use `code/step2d_inf_pair_power.r` to compute the power 

#### For pairwise transmission under varying epidemic dynamics:
* Set `param = param.pair.epi.dyn.unadj()` or `param = param.pair.epi.dyn.pow()`
* Use `code/step3d_inf_epi_dyn_power.r` to compute the power 

Data loading and model run script for cluster-level transmission is in `code/step4b_inf_cluster_simulate.r`. Calls the following R files:
> `code/step4a_inf_cluster_function.r` - simulates transmission pairs in household cluster of 3 persons with one index case<br/>
> `code/step2b_inf_pair_param.r` - set up parameters for cluster-level transmission 

#### For cluster-level transmission 
* Set `param = param.cluster()` 
* Use `code/step4c_inf_cluster_power.r` to compute the power 

### Plots
R scripts to perform the plots are:
> `code/step1c_haslemere_plot.R`<br/>
> `code/step2e_inf_pair_plot.R`<br/>
> `code/step3e_inf_epi_dyn_plot.R`<br/>
> `code/step4d_inf_cluster_plot.R`<br/>

If you plan to build on or cite this preliminary analysis for an academic publication, please ensure that you credit the underlying data sources above appropriately.
