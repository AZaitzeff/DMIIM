## DMIIM
This repository contains the code for the paper "On the Voronoi Implicit Interface Method" by Alexander Zaitzeff, Selim Esedoglu and Krishnakumar Garikipati.

# Code
The implementation of the parameterized curve code in sections four through six can be found in the para_curves folder.
The level set implementation of DMIIM is in the level_set_DMIIM folder. 

Not included is redisiting code on a grid (only needed for the level set implementation). We used https://github.com/mels630/fast-accurate-redistancing, included in utility/redist/ are the files we used to wrap the aforementioned redisting code into mex, including the mex command. The redisting code requires a linux operating system with gcc-4.9.x. To use, unzip the redisiting code into a folder called redist in utility/, then put mexmake.m, mexmake3.m, redistz.cpp, redistz3.cpp in the redist folder. Then in the test folder run 'setupredist'.

If you would like to use you own redisting code I mark the lines in dictmapping3.m and dictmapping2.m. In the tests files DMIIM2, DMIIM3, runtestgrimreaperlevelset.m we also used redisting for initialization (though it is not required).

# Simulations from the paper  
The test folder contains code to generate the simulations from the paper. The parameterized tests specify what table number they correspond to.
