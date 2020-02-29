wrf_crtm_radiance_simulator
Bryan M. Karpowicz NASA GMAO/GESTAR/USRA

This is a basic simulator taking a wrf output file and running it through IASI B1 (band 1), and B2 (band 2) to come up with a smaller subset of channels.

Basically, I'm putting this here to make something cleaner in the future with examples of how to use pycrtm. 
Note, to get all you need for submodules you'll need:

```shell
git clone --recursive https://github.com/karpob/wrf_crtm_radiance_simulator
```
Also, you'll need to install:

https://wrf-python.readthedocs.io/en/latest/

Along with a link to rttov, to get its hdf5 profiles where I take ozone to fill in a US standard atmosphere, since the WRF runs I was working with didn't have any ozone fields populated. You could easily pull something from MERRA-2, I suppose, but I was doing this quickly.
The scripts in this directory (in the order which I would run for large WRF domains):

wrfSplit.py --> use this to split a wrf run at a single time step into spatial chunks in the North-South direction.

wrfSim_flat_B2B3.py  --> this script runs CRTM over a given time step, follow help args (or see the run.j for an example).

catBt.py  --> assembles the output of serveral wrfSim_flat_B2B3.py runs (basically cats together runs when using wrfSplit as a fisrt step for large WRF domains).

plot.py  --> plots the output from wrfSim_flat_B2B3.py or catBt.py

run.j --> example slurm job to run on discover. Obviously, I'm selfish and have things hard coded like groups I'm a member of...etc.
 
