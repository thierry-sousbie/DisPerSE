#include "$HOME/prog/morse/yorick/network.i"
#include "$HOME/Yorick/Chris/skel.i"
#include "$HOME/prog/morse/yorick/skel.i"
#include "$HOME/Yorick/Chris/gadget.i"
#include "yeti.i"
#include "$HOME/Yorick/Chris/randfield.i"
#include "$HOME/Yorick/Eric/histo.i"

NGEN = int(-1);
NDIMS = 3;
N = 64;
NPART = N^NDIMS;
NBINS=200;
POWER_INDEX = 2;

/*
NGEN = 1;
NDIMS = 2;
N = 128;
NPART = N^NDIMS;
NBINS=200;
*/

func __PS(s) {if (POWER_INDEX==0) return 1; else return 1./(1e-12+s^POWER_INDEX);}

func generate(Ngen)
{
    if (is_void(Ngen)) Ngen = NGEN;
    seed = randomize(void);
    seedI = long(seed*(1<<31));
    dms= [NDIMS,NDIMS,NPART];
    while (numberof(dms)!=NDIMS+1) {grow,dms,1;}
    name_list=[];
    
    nreal=int(1);
    for (nreal=int(1);nreal!=1+Ngen;nreal++) {
        //d=random(dms);
        d=array(float,dms);
        if (NDIMS==3) d(,,*)=transpose(Dgenrandfield3D(NPART,pp=N,fun=__PS));
        else if (NDIMS==2) d(,,*)=transpose(Dgenrandfield2D(NPART,pp=N,fun=__PS));
        //nreal=int(1);
        nreal;
        fname = swrite(format="P%d_%d_%1.2f_S%9d-%3.3d.ND",NDIMS,N,POWER_INDEX,seedI,nreal);
        
        NDfield_write(d,fname);
        if (NDIMS==2) cmd = "delaunay_2D -margin 0.4 -periodic "+fname;
        else cmd = "delaunay_3Dp "+fname;
        write,"\nExecuting :",cmd;
        out=exec(cmd);
        write,"Output:";write,out;
        net_fname = fname+".NDnet";
        
        cmd = "mse -ppairs -nsig 0 -no_manifolds -no_arcsGeom -no_skl "+net_fname;
        write,"\nExecuting :",cmd;
        out=exec(cmd);
        write,"Output:";write,out;
        
        cmd="rm "+fname+" "+net_fname;
        write,"\nExecuting :",cmd;
        out=exec(cmd);
        fname = net_fname+".ppairs.NDnet";
        grow,name_list,pwd()+fname;
        //write,"Output:";write,out;
    }

    return name_list;
}

func loadPairs(fnames,dir=)
{
    if (is_void(dir)) dir="./";
    else dir+="/";
    
    if (is_void(fnames)) {
        fnames = exec("ls "+dir+"*.ppairs.NDnet");
    }
    
    local net;
    
    net=[];
    
    for (i=1;i<=numberof(fnames);i++) {
        write,"Loading ",fnames(i);
        grow,net,NDnetwork_read(fnames(i));
    }

    return net;        
}

func dataFromNet(&net, index) {
    return (*net.data)(index).data;
}

func binnedIndiceInt(val, &index2val, uniqueVal=)
{
    if ((typeof(val)!="long")&&(typeof(val)!="int"))
    {
        error,"in binnedIndiceInt, parameter must be an int or a long";
    }

    if (is_void(uniqueVal)) uv=unique(val);
    else uv=uniqueVal;

    delta=min(uv);
    if (delta>=1) delta =0;
    else delta = -(delta-1);
    
    val2index = array(int,max(uv)+1-min(uv));
    index2val = array(int,numberof(uv));
    for (i=1;i<=numberof(uv);i++)
    {
        val2index(delta+uv(i)) = i;
        index2val(i)=uv(i);
    }
    
    nin=array(int(0),numberof(uv));
    
    for (i=1;i<=numberof(val);i++) {
        index = val2index(delta+val(i));
        nin(index)++;
    }

    result=array(pointer,numberof(uv));
    for (i=1;i<=numberof(result);i++)
        result(i)=&array(int(),nin(i));
    
    nin=array(int(0),numberof(uv));
    
    for (i=1;i<=numberof(val);i++) {
        index = val2index(delta+val(i));
        (*result(index))(nin(index)) = i;
        nin(index)++;
    }

    return result;
}

struct pdf {
    pointer ph, rh, dh, rdh;
    pointer px, rx, dx, sdx;
}

func analyse( npart, dir, fnames=, dw=, nopreload=, ravg=)
{
     
    if (is_void(nopreload)) nopreload=0;
    if (is_void(dir)) dir="./pairs";
    
    if (is_void(fnames)) {
        lst = swrite(format="%s/*_%d_*ppairs.NDnet",dir,npart);
        fnames = exec("ls "+lst);
    }
    
    if (!nopreload) {
        all_net = loadPairs(fnames);
        net=all_net(1);
    }
    else net = NDnetwork_read(fnames(1));
    
    if (is_void(ravg))
        ravg = (npart^NDIMS);
    
    vti = NDNetDataIndex(net(1),0,"type");
    vdi = NDNetDataIndex(net(1),0,"density");
    vpi = NDNetDataIndex(net(1),0,"persistence");
    vri = NDNetDataIndex(net(1),0,"persistenceR");

    sti = NDNetDataIndex(net(1),1,"type");
    spi = NDNetDataIndex(net(1),1,"persistence");
    sri = NDNetDataIndex(net(1),1,"persistence Ratio");
    siui = NDNetDataIndex(net(1),1,"index up");
    sidi = NDNetDataIndex(net(1),1,"index down");
    //dataFromNet(net(1), vti);

    pstat=array(double,[3,numberof(fnames),NDIMS,4]);
    rstat=array(double,[3,numberof(fnames),NDIMS,4]);
    dstat=array(double,[3,numberof(fnames),NDIMS+1,4]);   
    sdstat = array(double,[3,numberof(fnames),NDIMS,4]);   
    
    for (i=1;i<=numberof(fnames);i++)
    {
        if (!nopreload) net = all_net(i);
        else net = NDnetwork_read(fnames(i));
        
        st = int(*dataFromNet(net(1), sti));
        vt = int(*dataFromNet(net(1), vti));
        sp = (*dataFromNet(net(1), spi))/ravg;
        sr = (*dataFromNet(net(1), sri));
        vd = (*dataFromNet(net(1), vdi));vd=(vd-ravg)/ravg +1;
        siu = int(1+(*dataFromNet(net(1), siui)));
        sid = int(1+(*dataFromNet(net(1), sidi)));
        write,"Realisation : ",i;
        
        if (NDIMS==3) stb = binnedIndiceInt(st,typeIndex,uniqueVal=int([1,2,3]));
        else stb = binnedIndiceInt(st,typeIndex,uniqueVal=int([1,2]));
        for (j=1;j<=numberof(stb);j++) {
            //write,"t:",type2index(j);
            cur = *stb(j);
            
            pstat(i,j,)=[min(sp(cur)),max(sp(cur)),avg(sp(cur)),RMS(sp(cur))];
            rstat(i,j,)=[min(sr(cur)),max(sr(cur)),avg(sr(cur)),RMS(sr(cur))];
            
            sdstat(i,j,)=[min(vd(siu(cur))),max(vd(siu(cur))),avg(vd(siu(cur))),RMS(vd(siu(cur)))];  
        }

        if (NDIMS==3) vtb = binnedIndiceInt(vt,typeIndex,uniqueVal=int([0,1,2,3]));
        else vtb = binnedIndiceInt(vt,typeIndex,uniqueVal=int([0,1,2]));
        for (j=1;j<=numberof(vtb);j++) {
            //write,"t:",type2index(j);
            cur = *vtb(j);
            dstat(i,j,)=[min(vd(cur)),max(vd(cur)),avg(vd(cur)),RMS(vd(cur))];
            
        }
            
    }
    
    px=spanl(min(pstat(,,1)(*)),max(pstat(,,2)(*)),NBINS+1);
    rx=spanl(min(rstat(,,1)(*)),max(rstat(,,2)(*)),NBINS+1);
    dx=spanl(min(dstat(,,1)(*)),max(dstat(,,2)(*)),NBINS+1);
    sdx=spanl(min(sdstat(,,1)(*)),max(sdstat(,,2)(*)),NBINS+1);
    
    rh=array(double,[3,numberof(fnames),NDIMS,NBINS]);
    ph=array(double,[3,numberof(fnames),NDIMS,NBINS]);
    dh=array(double,[3,numberof(fnames),NDIMS+1,NBINS]);
    rdh=array(double,[4,numberof(fnames),NDIMS+1,NBINS,NBINS]);
    
    for (i=1;i<=numberof(fnames);i++)
    {
        if (!nopreload) net = all_net(i);
        else net = NDnetwork_read(fnames(i));

        st = int(*dataFromNet(net(1), sti));
        vt = int(*dataFromNet(net(1), vti));
        sp = (*dataFromNet(net(1), spi))/ravg;
        sr = (*dataFromNet(net(1), sri));
        vd = (*dataFromNet(net(1), vdi));vd=(vd-ravg)/ravg +1;
        siu = int(1+(*dataFromNet(net(1), siui)));
        sid = int(1+(*dataFromNet(net(1), sidi)));
        
        write,"Realisation : ",i;
        if (NDIMS==3) stb = binnedIndiceInt(st,typeIndex,uniqueVal=int([1,2,3]));
        else stb = binnedIndiceInt(st,typeIndex,uniqueVal=int([1,2]));
        
        
        for (j=1;j<=numberof(stb);j++) {
            //write,"t:",type2index(j);
            cur = *stb(j);

            ph(i,j,) = double(histo1d(sp(cur),px));
            //ph(i,j,)/=ph(i,j,sum);
            
            rh(i,j,) = double(histo1d(sr(cur),rx));
            //rh(i,j,)/=rh(i,j,sum);

            rdh(i,j,,) = double(histo2d([sr(cur),vd(siu(cur))],rx,sdx));
            //rdh(i,j,,) /= rdh(i,j,*)(sum);
        }

        if (NDIMS==3) vtb = binnedIndiceInt(vt,typeIndex,uniqueVal=int([0,1,2,3]));
        else vtb = binnedIndiceInt(vt,typeIndex,uniqueVal=int([0,1,2]));
        
        for (j=1;j<=numberof(vtb);j++) {
            cur = *vtb(j);

            dh(i,j,) = double(histo1d(vd(cur),dx));
            //dh(i,j,)/=dh(i,j,sum);
        }
        
        //s=dh(i,sum,sum);
        //dh(i,,)/=numberof(vt);
        
    }

    result = pdf();
    result.ph=&ph;
    result.rh=&rh;
    result.dh=&dh;
    result.rdh=&rdh;
    
    result.px=&px;
    result.rx=&rx;
    result.dx=&dx;
    result.sdx=&sdx;

    return result;
}

func plotResults(result,dw=,cr=,overplot=,width=, noPDF=)
{
    if (is_void(width)) width=2;
    if (is_void(dw)) dw=0;
    if (is_void(overplot)) overplot=0;
    if (is_void(cr)) cr=1.0;
    if (is_void(noPDF)) noPDF=0;
    
    dw*=(3+NDIMS-1);
    
    ph=*result.ph;
    rh=*result.rh;
    dh=*result.dh;
    rdh=*result.rdh;
    
    px=*result.px;
    rx=*result.rx;
    dx=*result.dx;
    sdx=*result.sdx;

    nreal = double(dimsof(ph)(2));
    
    if (!overplot) WS,dw+1,dpi=150;
    else window,1;
    if (!noPDF) {
    ra=rh(avg,,);
    nrm=(ra(,zcen)*(rx(zcen)(dif))(-,))(,sum);
    
    if (NDIMS==3) {
    w=where(ra(3,) >=2./nreal);
    plg,ra(3,w)/nrm(3),rx(zcen)(w),width=width,color=char(cr*(__red)),type=2;
    }
    w=where(ra(2,) >=2./nreal);
    plg,ra(2,w)/nrm(2),rx(zcen)(w),width=width,color=char(cr*(__blue)),type=2;
    w=where(ra(1,) >=2./nreal);
    plg,ra(1,w)/nrm(1),rx(zcen)(w),width=width,color=char(cr*(__green)),type=2;
    }
    logxy,1,1;
    
    
    if (!overplot) WS,dw+2,dpi=150;
    else window,2;
    if (!noPDF) {
    pa=ph(avg,,);
    nrm=(pa(,zcen)*(px(zcen)(dif))(-,))(,sum);
    if (NDIMS==3) {
    w=where(pa(3,) >=2./nreal);
    plg,pa(3,w)/nrm(3),px(zcen)(w),width=width,color=char(cr*(__red)),type=2;
    }
    w=where(pa(2,) >=2./nreal);
    plg,pa(2,w)/nrm(2),px(zcen)(w),width=width,color=char(cr*(__blue)),type=2;
    w=where(pa(1,) >=2./nreal);
    plg,pa(1,w)/nrm(1),px(zcen)(w),width=width,color=char(cr*(__green)),type=2;
    logxy,1,1;
    }
    
    //WS,3,dpi=150;
    if (!overplot) window,dw+1;
    else window,1;
    ra=rh(avg,,);
    it=ra(,::-1)(,psum)(,::-1)
    if (NDIMS==3) {
    w=where(ra(3,) >=2./nreal);
    plg,(it(3,w)/it(3,w(1)))(zcen),rx(:-1)(w)(zcen),width=width,color=char(cr*(__red));
    }
    w=where(ra(2,w) >=2./nreal);
    plg,(it(2,w)/it(2,w(1)))(zcen),rx(:-1)(w)(zcen),width=width,color=char(cr*(__blue));
    w=where(ra(1,) >=2./nreal);
    plg,(it(1,w)/it(1,w(1)))(zcen),rx(:-1)(w)(zcen),width=width,color=char(cr*(__green));
    logxy,1,1;
    
    //WS,4,dpi=150;
    if (!overplot) window,dw+2;
    else window,2;
    pa=ph(avg,,);
    it=pa(,::-1)(,psum)(,::-1)
        //nrm=(pa(,zcen)*(px(zcen)(dif))(-,))(,sum);
    if (NDIMS==3) {
    w=where(pa(3,) >=2./nreal);
    plg,it(3,w)/it(3,w(1)),px(zcen)(w),width=width,color=char(cr*(__red));
    }
    w=where(pa(2,) >=2./nreal);
    plg,it(2,w)/it(2,w(1)),px(zcen)(w),width=width,color=char(cr*(__blue));
    w=where(pa(1,) >=2./nreal);
    plg,it(1,w)/it(1,w(1)),px(zcen)(w),width=width,color=char(cr*(__green));
    logxy,1,1;

    if (!overplot) WS,dw+3,dpi=150;
    else window,3;
    da=dh(avg,,);
    nrm=(da(,zcen)*(dx(zcen)(dif))(-,))(,sum);
    if (NDIMS==3) {
        plg,da(4,)/nrm(sum),dx(zcen),width=width,color=char(cr*(__red));
    }
    plg,da(3,)/nrm(sum),dx(zcen),width=width,color=char(cr*(__yellow));
    plg,da(2,)/nrm(sum),dx(zcen),width=width,color=char(cr*(__green));
    plg,da(1,)/nrm(sum),dx(zcen),width=width,color=char(cr*(__blue));
    logxy,1,0;


    if (NDIMS==3) id=3;
    else id=2;

    for (;id>0;id--) {
        WS,dw+4+id-1,dpi=150;
    
    
        rda=rdh(avg,,,);
        nrm=rda(,sum,sum);
        w=where(rda(id,,)(*) != 0);
        minimum = min(rda(id,,)(*)(w));
        rda(id,,) += minimum/10;
    
        plcont,log(rda(id,,)/nrm(id)),sdx(-::NBINS-1,zcen),rx(zcen,-::NBINS-1);
        logxy,1,1;
        c_c;
    }
}

func __fitmode_t1(x, a, &grad, deriv=, getorder=, getsuggest=)
{
  if (!is_void(getorder)) return int(3);
  if (!is_void(getsuggest)) {
      if (NDIMS==2) return [2,0.1,2];
      else return [3.,1.,2.];
  }
  
  ff=exp(-a(1)*(x-1)-a(2)*pow((x-1),a(3)));
  if (deriv) grad=ff*[-(x-1),-pow((x-1),a(3)),-a(2)*log(x-1)*pow(x-1,a(3))];
    
  return (ff);
}
func __fitmode_t1_log(x, a, &grad, deriv=, getorder=, getsuggest=)
{
    return log(__fitmode_t1(x, a, grad, deriv=deriv, getorder=getorder, getsuggest=getsuggest));
}

func __fitmode_t2(x, a, &grad, deriv=, getorder=, getsuggest=)
{
  if (!is_void(getorder)) return int(2);
  
  if (!is_void(getsuggest)) return [2.55316];

  u=x-1;
  a2=4;
  a3=9.;
  a4=1.785;
  ff1 = exp(-a(1)*u);
  ff2 = a2*x^(-a3);
  
  t = 1./(1.+(a4/u)^14);

  ff = (ff2*t) + (1.-t)*ff1;
  if (deriv)
  {
      grad=[(1.-t)*(-u*exp(-a(1)*u)),
            x^(-a3)*t];
  }
  
  return (ff);
}
func __fitmode_t2_log(x, a, &grad, deriv=, getorder=, getsuggest=)
{
    return log(__fitmode_t2(x, a, grad, deriv=deriv, getorder=getorder, getsuggest=getsuggest));
}

func __fitmode_t3(x, a, &grad, deriv=, getorder=, getsuggest=)
{
   
  if (!is_void(getorder)) return int(2);
  if (!is_void(getsuggest)) return [0.5,3.];
   
  u=(x-1);
  
  ff=pow(1+a(1)*u,-a(2));
  /*
  if (deriv) grad=[exp(-a(2)*u),
                   -a(1)*u*exp(-a(2)*u),
                   pow(u,-a(4))*exp(-a(5)/u),
                   -a(4)*a(3)*pow(u,-a(4)-1)*exp(-a(5)/u),
                   a(3)*pow(u,-a(4))*(-1/u)*exp(-a(5)/u)];
  */
  return (ff);
}
func __fitmode_t3_log(x, a, &grad, deriv=, getorder=, getsuggest=)
{
    return log(__fitmode_t3(x, a, grad, deriv=deriv, getorder=getorder, getsuggest=getsuggest));
}


func doFit(y,x,type,wei=,verb=,plot=, firstGuess=, uselog=)
{
    if (is_void(uselog)) uselog=0;
    if (is_void(type)) type=1;
    //if (type>3) type=3;
    uselog=0;
    if (type==2*1) {
        extern __fitmode_t1;
        fitfunc=__fitmode_t1;
        if (uselog) fitfunc=__fitmode_t1_log;
        
    } else if ((type==2*2)&&(NDIMS==3)) {
        extern __fitmode_t2;
        fitfunc=__fitmode_t2;
        if (uselog) fitfunc=__fitmode_t2_log;
        
    } else if (((type==2*3)&&(NDIMS==3))||((type==2*2)&&(NDIMS==2))) {
        extern __fitmode_t3;
        fitfunc=__fitmode_t3;
        //uselog=0;
        if (uselog) fitfunc=__fitmode_t3_log;
    }
    
    nparms = fitfunc(getorder=1);
    a=fitfunc(getsuggest=1);
    
    if (is_void(wei))
        wei=y*0.+1;
    wei(:-1) = (x(dif)/(x(dif)(max)));
    //wei(int(numberof(wei)*0.5):) = 0;
    
    d=[];
    fitfunc,x,a,d,deriv=1;
    if (numberof(d) == nparms)
    {
        if (uselog)
            r=lmfit(fitfunc,x,a,log(y),wei,deriv=1,  stdev=1, monte_carlo=500, correl=1);
        else
            r=lmfit(fitfunc,x,a,y,wei,deriv=1,  stdev=1, monte_carlo=500, correl=1);
    }
    else
    {
        if (uselog)
            r=lmfit(fitfunc,x,a,log(y),wei);
        else
            r=lmfit(fitfunc,x,a,y,wei);
    }
    
    if(verb) {
        if (!is_void(plot))
        {
            //r;
            //plg,y,x;
        }
        r;
        a;
        w=where(wei!=0);
        plg,fitfunc(x,a),x,color=__purple,type=2,width=2;
    }
    return a;
}

func fitProba(result,noplot=,uselog=)
{
    if (is_void(noplot)) noplot=0;
    
    ph=*result.ph;
    rh=*result.rh;
    dh=*result.dh;
    rdh=*result.rdh;
    
    px=*result.px;
    rx=*result.rx;
    dx=*result.dx;
    sdx=*result.sdx;
    parms=[];
    type=3;
    for (type=1;type<=dimsof(rh)(3);type++) {
        
        write,"Fitting type ",type;
        ra=rh(avg,,);
        it=ra(,::-1)(,psum)(,::-1);
        for (i=1;i<=dimsof(it)(2);i++) {
            w=where(ra(i,) != 0);
            it(i,)/=it(i,w(1));        
        }
        w=where(ra(type,) != 0);
        y=it(type,w)(zcen);
        x=rx(:-1)(w)(zcen);
        //y;x;
        if (noplot)
            a=doFit((y),x,type*2,uselog=uselog);
        else a=doFit((y),x,type*2,verb=1,uselog=uselog);
        grow,parms,&a;
    }

    return parms;
}

func publiPlot(results,win=,save=,overplot=)
{
    if (is_void(win)) win=window();
    
    if (is_void(overplot)) WS,win,dpi=150;
    width=numberof(results)+2;
    cr=1.00-0.2*(numberof(results)-1);

    plg,[erfc(1/sqrt(2)),erfc(1/sqrt(2))],[1,500],type=2;
    plg,[erfc(2/sqrt(2)),erfc(2/sqrt(2))],[1,500],type=2;
    plg,[erfc(3/sqrt(2)),erfc(3/sqrt(2))],[1,500],type=2;
    plg,[erfc(4/sqrt(2)),erfc(4/sqrt(2))],[1,500],type=2;
    plg,[erfc(5/sqrt(2)),erfc(5/sqrt(2))],[1,500],type=2;
    
    for (i=1;i<=numberof(results);i++) {
        result = results(i);
        
        rh=*result.rh;
        rx=*result.rx;
        nreal = double(dimsof(rh)(2));
        
        ra=rh(avg,,);
        it=ra(,::-1)(,psum)(,::-1)
            if (NDIMS==3) {
                w=where(ra(3,) >=2./nreal);
                plg,(it(3,w)/it(3,w(1)))(zcen),rx(:-1)(w)(zcen),width=width,color=char(cr*(__red));
            }
        w=where(ra(2,w) >=2./nreal);
        plg,(it(2,w)/it(2,w(1)))(zcen),rx(:-1)(w)(zcen),width=width,color=char(cr*(__blue));
        w=where(ra(1,) >=2./nreal);
        plg,(it(1,w)/it(1,w(1)))(zcen),rx(:-1)(w)(zcen),width=width,color=char(cr*(__green));
        logxy,1,1;
        width-=1;
        cr+=0.2;      
    }
    
    //a=fitProba(results(0),noplot=1,uselog=1);
    a=fitProba(results(0),noplot=1);
    
    w=where(ra(1,)>2./nreal);
    x=spanl(min(rx(w)),max(rx(w))+(max(rx(w))-min(rx(w)))/10.,100);
    plg,__fitmode_t1(x, *a(1)),x,type=2,color=__black,width=2;
    *a(1);
    if (NDIMS==3)
    {
        w=where(ra(2,) >2./nreal);
        x=spanl(min(rx(w)),max(rx(w))+(max(rx(w))-min(rx(w)))/10.,100);
        plg,(__fitmode_t2(x, *a(2))),x,type=2,color=__black,width=2;
        *a(2);
    
        w=where(ra(3,) >2./nreal);
        x=spanl(min(rx(w)),max(rx(w))+(max(rx(w))-min(rx(w)))/10.,100);
        plg,(__fitmode_t3(x, *a(3))),x,type=2,color=__black,width=2;
        *a(3);
    }
    else if (NDIMS==2)
    {
        w=where(ra(2,) >2./nreal);
        x=spanl(min(rx(w)),max(rx(w))+(max(rx(w))-min(rx(w)))/10.,100);
        plg,(__fitmode_t3(x, *a(2))),x,type=2,color=__black,width=2;
        *a(2);
    }

    xytitles,"persistence ratio r","P(r)";

    if (!is_void(save)) eps,"Pr_PDF.ps";
}

func plotPDF(d,N)
{
//d=read_ascii("denisty_64.dat");
    e=double(d);
    davg=double(N^3);e=(e-davg)/davg+1;
    x=spanl(min(e),max(e),200);
    h=double(histo1d(e,x));
    w=where(h!=0);
    v=h(w)/(h(w)*x(dif)(w))(sum);
    u=x(zcen)(w);
    plg,v,u,width=2;
}

func plotPS(xyz_p,bsize=,color=)
{
    if (dimsof(xyz_p)(0)>dimsof(xyz_p)(-1))
        xyz=transpose(xyz_p);
    else xyz=xyz_p;
    if (is_void(color)) color=__red;
    if (is_void(bsize)) bsize=1.;

    pp=1+int(pow(double(dimsof(xyz)(-1)),1./dimsof(xyz)(0)));
    if (pow(pp,dimsof(xyz)(0)) != dimsof(xyz)(-1)) pp++;
    if (pow(pp,dimsof(xyz)(0)) != dimsof(xyz)(-1)) error,"could not determine box size.";
    h=bsize * double(indgen(pp))/pp;
    h;
    u=histo3d((xyz),h,h,h);
    u /= double(sum(u));
    //stat,tt;tt(*)(sum);
    z=fft(u,1)/numberof(u);
    ns=dimsof(u)(2);
    z2=float(fft(z*conj(z),[-1,-1,-1]));
    x1=(indgen(ns)-1)*2*pi/float(ns);
    y1=(indgen(ns)-1)*2*pi/float(ns);
    z1=(indgen(ns)-1)*2*pi/float(ns);
    r1= abs(x1, y1(-,),z1(-,-,));
    zc=array(0.,ns,ns,ns);z=fft(u,1)/numberof(u); z2=float(z*conj(z));zc+=z2;
    py= histo2(r1, px, weight=zc, average=1, interp=1,binsize=r1(dif)(min));
    plg, py(2:), px(2:), color=color; //azimuthal average
    //plg, py2,px(2:),color="blue",width=4;
    logxy,1,1;
}
/*

  cd,"~/work/morse/poisson/";
  #include "PoissonPersistence.i"
  
  //r32=analyse(32);
  //r64=analyse(64);
  //r128=analyse(128,nopreload=1);


upload,"poisson_results.pdb"

  
  //plotResults(r128,width=3,cr=0.5);
  //plotResults(r64,width=2,overplot=1,cr=0.75,dw=1);
  //plotResults(r32,overplot=1,width=1,dw=2);

  plotResults(r32,width=4,noPDF=1);
  plotResults(r64,overplot=1,cr=0.75,dw=1,width=3,noPDF=1);
  plotResults(r128,overplot=1,cr=[0.5,0.5,0.5],dw=2,width=2,noPDF=1);
  window,1;
  fitProba,r128;
 */


/*
  
//r0_128=analyse( 128, "./pairs_0/", nopreload=1);
//r1_128=analyse( 128, "./pairs_1/", nopreload=1);
//r2_128=analyse( 128, "./pairs_2/", nopreload=1);
//save,createb("poisson_k_results.pdb"),r0_128,r1_128,r2_128;

upload,"poisson_k_results.pdb"

plotResults(r0_128,width=4,noPDF=1);
plotResults(r1_128,overplot=1,cr=0.75,dw=1,width=3,noPDF=1);
plotResults(r2_128,overplot=1,cr=[0.5,0.5,0.5],dw=2,width=2,noPDF=1);

 */
