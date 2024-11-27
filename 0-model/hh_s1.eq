x'=px;
y'=py;
px'=-x-2*x*y;
py'=-y-(x*x-y*y);

expr pxval=y*y*(2.0*y/3.0-1) - py*py;

jet all symbols 1 deg 8;


/*
taylor -o hh_s1.h -name ode -expression -header hh_s1.eq
taylor -o hh_s1.c -name ode -expression -jet -step -headername ode.h hh_s1.eq


taylor -o hh_s1_cmplx.h -name ode -complex -expression -header hh_s1.eq
taylor -o hh_s1_cmplx.c -name ode -complex -expression -jet -step -headername ode.h hh_s1.eq
 */
