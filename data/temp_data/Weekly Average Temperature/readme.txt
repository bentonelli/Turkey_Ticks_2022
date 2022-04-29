These site-specific weekly temperature files were cleaned and calculated from all available temperature data of the 9 sites I could found from NEON website. 

For each site, in order to minimize the affect of missing values, I merged all files with 30min temperature logs. All the integrating functions I used automatically ignored NA entries when calculating stuffs. However, there were still NA's due to some missing data over a long period of time. To deal with this, I tried the auto completion method from the MICE package and it worked out pretty well.

I plotted the temperature mean v.s. week plots to check the data accuracy. They all look like sine waves which conform to my common sense.
