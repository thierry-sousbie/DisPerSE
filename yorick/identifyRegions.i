#include "~/prog/octree2/yorick/octree.i"

/*
  name = "fof_100Mpc_512_z0_b0p2_y";
  disp_name = "snap100Mpc128_z0.NDnet_s6";
  skl=NDskel_read(disp_name+".up.NDskl.BRK.NDskl");
  wall=NDnetwork_read(disp_name+"_manifolds_J1a.NDnet");
  pos=read_ascii(name)(2:4,)*1000.;
  dlim=2000.0;
  net=identifyRegion(pos,skl,dlim,wall=wall,boxMin=0.0,boxMax=100000.0);
  NDnetwork_write(net,name+"_webID");

  pos=read_ascii(name);
  pp=transpose(pos);grow,pp,double(*((*net.data)(1).data));
  smwrite(name+"_webID",pp,head="      x     y    z h^-1 Mpc         velocity[km/s]      lkl   mass h^-1 M_s  sigma  sigma_v  r_sph  delta  spin p. webID");
*/

func treatFiles(name,disp_name,dlim=,boxMin=,boxMax=,scaleInput=,outName=)
{
  if (is_void(dlim)) error,"Need dlim param !";
  if (is_void(boxMax)) error,"Need boxMax param !";
  if (is_void(boxMin)) boxMin=0;
  
  boxMin=double(boxMin);boxMax=double(boxMax);
  
  skl=NDskel_read(disp_name+".up.NDskl.BRK.NDskl");
  wall=NDnetwork_read(disp_name+"_manifolds_J1a.NDnet");
  pos=read_ascii(name)(2:4,);
  if (!is_void(scaleInput)) pos*=scaleInput;
  
  net=identifyRegion(pos,skl,dlim,wall=wall,boxMin=boxMin,boxMax=boxMax);
  NDnetwork_write(net,name+"_webID",toVtu=1);//_webID

  if (is_void(outName)) outName=name+"_webID";
  
  pos=read_ascii(name);
  pp=transpose(pos);grow,pp,double(*((*net.data)(1).data));
  smwrite("FOF_"+outName,pp,head="      x     y    z h^-1 Mpc         velocity[km/s]      lkl   mass h^-1 M_s  sigma  sigma_v  r_sph  delta  spin p. webID",fmt="%g\t");
}

func identifyRegion(pos,skl,dlim,wall=,boxMin=,boxMax=)
{  
  id=array(0,dimsof(pos)(0));

  if (!is_void(wall))
    {
      segp = (*wall.v_coord);
      if (!is_void(boxMin)) {
        w=where(segp(*)<boxMin);
        if(numberof(w)) segp(*)(w)+=(boxMax-boxMin);
      }
      if (!is_void(boxMax)) {
        w=where(segp(*)>=boxMax);
        if(numberof(w)) segp(*)(w)-=(boxMax-boxMin);
      }
      
      tree=BuildOctreeFromPos(segp,Min=boxMin,Max=boxMax);
      dist=[];
      nei=FindNeighbours(tree,pos,1,dist,dist=1,periodic=1);     
      id=array(0,dimsof(pos)(0));
      w=where(dist<dlim);
      id(w)=1;
    }
  
  segp=(*skl.segpos)(,avg,);
  if (!is_void(boxMin)) {
    w=where(segp(*)<boxMin);
    if(numberof(w)) segp(*)(w)+=(boxMax-boxMin);
  }
  if (!is_void(boxMax)) {
    w=where(segp(*)>=boxMax);
    if(numberof(w)) segp(*)(w)-=(boxMax-boxMin);
  }
  
  tree=BuildOctreeFromPos(segp,Min=boxMin,Max=boxMax);
  dist=[];
  nei=FindNeighbours(tree,pos,1,dist,dist=1,periodic=1);
  FreeOctree(tree);
  w=where(dist<dlim);
  id(w)=2;
  
  net=NDnetwork();
  net.comment="vertex regions";
  net.periodicity=(1<<4)-1;
  net.ndims=dimsof(pos)(2);
  net.ndims_net=dimsof(pos)(2);
  net.f_areSimplex=1;
  net.x0=&array(double(boxMin),net.ndims);
  net.delta=&array(double(boxMax-boxMin),net.ndims);
  net.nvertex=dimsof(pos(,*))(0);
  net.v_coord=&pos;
  data=array(NDnetworkData(),1);
  data(1).type=0;
  data(1).name="web_ID";
  data(1).data=&double(id);
  net.data=&data;
  nv=numberof(where(id==0));
  nw=numberof(where(id==1));
  nf=numberof(where(id==2));
  nt=numberof(id);
  write,format="Counted %ld(%2.1f%%) voids, %ld(%2.1f%%) wall and %ld(%.1f%%) fil over %ld.\n",nv,100.0*double(nv)/nt,nw,100.0*double(nw)/nt,nf,100.0*double(nf)/nt,nt;
  return net;
}

func treatAllFiles
{
  treatFiles("fof_100Mpc_512_z0_b0p2_y","snap100Mpc128_z0.NDnet_s6",dlim=2000.0,boxMax=100000,scaleInput=1000.0,outName="100_512_sousbie.txt");
  treatFiles("fof_100Mpc_256_z0_b0p2_y","snap100Mpc128_z0.NDnet_s6",dlim=2000.0,boxMax=100000,scaleInput=1000.0,outName="100_256_sousbie.txt");
  treatFiles("fof_100Mpc_128_z0_b0p2_y","snap100Mpc128_z0.NDnet_s6",dlim=2000.0,boxMax=100000,scaleInput=1000.0,outName="100_128_sousbie.txt");

  treatFiles("fof_200Mpc_512_z0_b0p2_y","snap200Mpc128_z0.NDnet_s6",dlim=4000.0,boxMax=200000,scaleInput=1000.0,outName="200_512_sousbie.txt");
  treatFiles("fof_200Mpc_256_z0_b0p2_y","snap200Mpc128_z0.NDnet_s6",dlim=4000.0,boxMax=200000,scaleInput=1000.0,outName="200_256_sousbie.txt");
  treatFiles("fof_200Mpc_128_z0_b0p2_y","snap200Mpc128_z0.NDnet_s6",dlim=4000.0,boxMax=200000,scaleInput=1000.0,outName="200_128_sousbie.txt");
  
}
