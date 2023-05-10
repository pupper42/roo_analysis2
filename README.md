# ROO Analysis (ver. 2)

Version 2 uses Skyfield instead of Astropy for the transformations.

To get this working you need to do some modifications:

* Georinex Package

    * Find where georinex is installed
        * If you use anaconda on Windows it should be `C:/Users/<name>/AppData/Local/anaconda3/envs/<env_name>/Lib/site-packages/georinex`
    * Comment out the line `ds.to_netcdf(outfn, mode="w", encoding=enc)` (which should be line 113)

* Numerical methods

    * Change `Y = Y[row0-nn:row0+nn+1, :]` in line 68 to `Y = Y[row0-nn:row0+nn+1]`
 
