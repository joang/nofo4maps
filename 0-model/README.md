# Model

It contains the raw equations of the Hénon-Heiles system, those with extension `.eq`.

 - `hh.eq` contains the ODE system with the expression to isolate the variable `py` given an energy level (external variable) and assuming `x=0`.
 - `hh_vars.eq` same as `hh.eq` but adding the line for jets with symbols and degree 1.
 - `hh_s2.eq` same as `hh_vars.eq` but arbitrary degree.
 - `hh_s1.eq` same as `hh_s2.eq` but 1 symbol.
 
In each of these files, we added the commented lines indicating what command to run with [taylor v2.2](https://github.com/joang/taylor2-dist).

From those commands, the other files in the folder are obtained. In particular, the files with patter `_cmplx` refer to complex number arithmetic.



