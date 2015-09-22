#! /usr/bin/yorick -batch
//#! /home/thierry/apps/yorick/Linux-x86_64/bin/yorick -batch
#include "$HOME/Yorick/Chris/utls.i"
#include "$HOME/work/morse/poisson/PoissonPersistence.i"

POWER_INDEX = 0;

args=get_argv(); // waits for arg
if (numberof(args)>=2) N=int(str2long(args(2)));
if (numberof(args)>=3) NGEN=int(str2long(args(3)));
if (numberof(args)>=4) POWER_INDEX=double(str2double(args(4)));
if (numberof(args)>=5) NDIMS=int(str2long(args(5)));
if (numberof(args)>=6) error,"error: 4 arguments or less: N, NGEN, POWER_INDEX and NDIMS";

NPART = N^NDIMS;
write,swrite(format="will generate %d %d^%d = %d realisations, index = %.2g",NGEN,N,NDIMS,NPART,POWER_INDEX);

func __PS(s) {if (POWER_INDEX==0) return 1; else return 1./(1e-12+s^POWER_INDEX);}

generate();
