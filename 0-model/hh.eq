x'=px;
y'=py;
px'=-x-2*x*y;
py'=-y-(x*x-y*y);

expr pxval=y*y*(2.0*y/3.0-1) - py*py;


/*
taylor -o hh.h -name ode -expression -header hh.eq
taylor -o hh.c -name ode -expression -jet -step -headername ode.h hh.eq
 */
