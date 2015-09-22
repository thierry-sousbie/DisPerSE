#include "~/prog/yhealpix/yhealpix/healpix.i"

func skl2hp(skl,nside=, ext=)
{
  if (is_void(nside)) nside=512;
  npix=nside2npix(nside);

  map=array(char(0),npix);

  if (is_void(ext)) {
    p=*skl.segpos;
    pa=p(,avg,);  
    for (i=1;i<=skl.nsegs;i++) {
      map(1+ang2pix_nest(nside,p(1,1,i),p(2,1,i))) = 1;
      map(1+ang2pix_nest(nside,p(1,2,i),p(2,2,i))) = 1;
      map(1+ang2pix_nest(nside,pa(1,i),pa(2,i))) = 1;
    }
  } else {
    p=*skl.nodepos;
    for (i=1;i<=skl.nnodes;i++)
      map(1+ang2pix_nest(nside,p(1,i),p(2,i))) = 1;
  }
  
  return map;
}
