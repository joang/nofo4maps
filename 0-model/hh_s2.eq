x'=px;
y'=py;
px'=-x-2*x*y;
py'=-y-(x*x-y*y);

expr pxval=y*y*(2.0*y/3.0-1) - py*py;

jet all symbols 2 deg 15;


/*
taylor -o hh_s2.h -name ode -expression -header hh_s2.eq
taylor -o hh_s2.c -name ode -expression -jet -step -headername ode.h hh_s2.eq

taylor -o hh_s2_cmplx.h -name ode -complex -expression -header hh_s2.eq
taylor -o hh_s2_cmplx.c -name ode -complex -expression -jet -step -headername ode.h hh_s2.eq
 */
