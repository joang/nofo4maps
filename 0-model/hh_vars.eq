x'=px;
y'=py;
px'=-x-2*x*y;
py'=-y-(x*x-y*y);

expr pxval=y*y*(2.0*y/3.0-1) - py*py;

jet all symbols 2 deg 1;


/*
taylor -o hh_vars.h -name ode -expression -jlib n_1 -header hh_vars.eq
taylor -o hh_vars.c -name ode -expression -jlib n_1 -jet -step -headername ode.h hh_vars.eq
 */
