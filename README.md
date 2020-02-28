wrf_crtm_radiance_simulator
Bryan M. Karpowicz NASA GMAO/GESTAR/USRA

This is a basic simulator taking a wrf output file and running it through IASI B1, and B2 to come up with a smaller subset of channels.

Basically, I'm putting this here to make something cleaner in the future with examples of how to use pycrtm. 
Note, to get all you need for submodules you'll need:

```shell
git clone --recursive https://github.com/karpob/wrf_crtm_radiance_simulator
```
Also, you'll need to install

https://wrf-python.readthedocs.io/en/latest/

Along with a link to rttov, to get it's hdf profiles where I take ozone to fill in a US Std atmosphere, since the WRF runs I was working with didn't have any ozone fields populated. You could easily pull something from MERRA-2, I suppose, but I was doing this quickly.

 
