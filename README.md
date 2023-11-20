# wavecode
Wave-like propagation of warps in accretion discs

## About
This code solves for the propagation of a warp wave through a Keplerian
disc. It uses the warp equations as given in Lubow, Ogilvie &
Pringle (2002) MNRAS 337, 706, Section 4, using the variables A and
D. For a Keplerian disc we set 'zeta' and 'eta' to be zero.

We can also add dissipation in the form of a decay of 'A' the radial
velocity component, using 'alpha'

Adapted from an original FORTRAN77 code by Jim Pringle
Modernised, converted to Fortran 90 by Daniel Price, April 2013

## Install
You will need a Fortran compiler, e.g. gfortran. On a Mac it's easiest to use homebrew:
```
brew install gfortran
```
then simply enter the code directory and type "make"
```
cd wavecode
make
```
this should create a binary called ``wave'' in the top-level directory

## Run
```
./wave
```
the code asks for a scenario to create a parameter file:
``` 
Please pick a mode:
  1 = blackhole
  2 = binary
  3 = binary-alpha0
```
following which you should find a setup.in file in the working directory:
```
% ls
Makefile	build/		setup.in	src/		wave*
```
Edit this file in your favourite text editor (e.g. nano, vi, emacs), and then rerun the binary with no arguments to run the code:
```
./wave
```
This should produce a series of output files as follows
```
% ./wave         
reading setup options from setup.in
 opening database from setup.in with 12 entries

 have set zi =                (0.0000000000000000,1.0000000000000000)

--- Frequencies ---------------------
 eta0    =   -4.4776119402985072E-002
 zeta0   =    4.4776119402985072E-002
 omega0  =   0.12789613706050704     

--- Disc ----------------------------
 H/R     =    5.0000000000000003E-002
 alphaSS =    2.0000000000000000E-002

--- Grid ----------------------------
 N     =                   300
 Rin   =    4.0000000000000000     
 Rout  =    160.00000000000000     
 rstep =    20.000000000000000     
 wstep =    2.0000000000000000     

--- Output --------------------------
 tstop    =    5000.0000000000000     
 toutfile =    125.00000000000000     
 ctime    =    2.9999999999999999E-002

START:
 nstep          0 time =   0.0000E+00 writing wave_00000

---------------------------------------------------------------------------
 Calculating timestep:
 timestep due eta at gridpt      2
 timestep dt =   5.4183E-03
---------------------------------------------------------------------------

 nstep      23071 time =   1.2501E+02 writing wave_00001
 cpu time since last dump =   0.14

 nstep      46142 time =   2.5001E+02 writing wave_00002
 cpu time since last dump =   0.13
 ...
```
## Plot
The snapshots are just ascii files with a header and data in each column. I used splash to visualise them:
```
brew tap danieljprice/all
brew install splash
splash wave_0* -y 4 -x 1
```
giving a result something like:
![image](https://github.com/danieljprice/wavecode/assets/12252103/3cab26ba-ce39-49fa-a6be-924c677c89d3)
