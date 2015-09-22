/*
 # Copyright (C) 2009-2013 Thierry Sousbie
 # University of Tokyo / CNRS, 2009
 #
 # This file is part of porject DisPerSE
 # 
 #  Author          : Thierry Sousbie
 #  Contact         : tsousbie@gmail.com	
 #
 #  Licenses        : This file is 'dual-licensed', you have to choose one
 #                    of the two licenses below to apply.
 #
 #                    CeCILL-C
 #                    The CeCILL-C license is close to the GNU LGPL.
 #                    ( http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html )
 #
 #                or  CeCILL v2.0
 #                    The CeCILL license is compatible with the GNU GPL.
 #                    ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
 #
 #  This software is governed either by the CeCILL or the CeCILL-C license
 #  under French law and abiding by the rules of distribution of free software.
 #  You can  use, modify and or redistribute the software under the terms of
 #  the CeCILL or CeCILL-C licenses as circulated by CEA, CNRS and INRIA
 #  at the following URL : "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and  rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL and CeCILL-C licenses and that you accept its terms.
*/
#ifndef __PD_DRAW_HXX__
#define __PD_DRAW_HXX__

#include <string>

#include <mgl/mgl.h>
#include <mgl/mgl_data.h>

#include <algorithm>
#include "PD_params.hxx"

class PD_draw //: public mglDraw
{
private:
  
  
public:
  static int Draw(mglGraph *gr, void *par)
  {
    static std::string markOpt[4] = {"b o#", "g s#", "r d#", "p v#"};
    static std::string markOutOpt[4] = {"k o", "k s", "k d", "k v"};

    PD_params *p = static_cast<PD_params *>(par);
    //NDnetwork *net = p->net;
    long i;
    //char prm[255];

    if (p->x_data==NULL) return gr->GetNumFrame();  
    if (p->y_data==NULL) return gr->GetNumFrame();  
    
    gr->ResetFrames();
    gr->NewFrame(); 

    //gr->Alpha(true);
    //gr->SetAlphaDef(0.5);

    gr->Rotate(0,0);
    gr->Axis(mglPoint(p->x1,p->y1,1),mglPoint(p->x2,p->y2,1));
    gr->SetMarkSize(0.01);
    gr->Box("",false);
    gr->Label('x',p->x_label.c_str(),0);
    gr->Label('y',p->y_label.c_str(),0);

 
    if (p->showHisto)
      {
	//gr->Dens(p->histo());
	//printf (" [%ld %ld %ld/%ld]\n",p->Hx.nx,p->Hy.nx,p->H.nx,p->H.ny);
	//gr->SetPalColor(0,1.0,1.0,1.0);

	gr->CRange(p->H);
	gr->Dens(p->Hx,p->Hy,p->H,"",-1);	
	gr->Colorbar();
	float level=(p->logH)?(log10(2)/30):(2./30);
	mglData levs;
	levs.Set(&level,1);
	gr->Cont(levs,p->Hx,p->Hy,p->H,"k",-0.99);
      }
    /*
    mglData a(50,40);
    a.Modify("0.6*sin(2*pi*x)*sin(3*pi*y) + 0.4*cos(3*pi*(x*y))");
    gr->Dens(a);
    */
    if (p->logX) 
      {
	if (p->logY) gr->SetFunc("lg(x)", "lg(y)", 0);
	else gr->SetFunc("lg(x)", 0, 0);
      }
    else if (p->logY) gr->SetFunc(0, "lg(y)", 0);

    gr->Axis("xy",true);


    
    
    //gr->Box();
    int n=0;
    
    for (std::list<int>::iterator it=p->shownP.begin();it!=p->shownP.end();it++)
      //for (i=0;i<4;i++)
      {
	if (p->hideP) continue;
	i=*it;
	//if (!p->showP[i]) continue;
	if (p->x_data[i].nx<=1) {/*printf("skipped %ld\n",i);*/continue;}
	if (p->y_data[i].nx<=1) {/*printf("skipped %ld\n",i);*/continue;}

	mglData z(p->x_data[i].nx);
	z.Fill(2*n*0.1,(2*n+1)*0.1);
	
	gr->SetMarkSize(0.02);
	gr->Plot(p->x_data[i],p->y_data[i],z,markOutOpt[i].c_str());//,0+i*0.1);

	//z+=0.1;

	gr->SetMarkSize(0.02);
	gr->Plot(p->x_data[i],p->y_data[i],z,markOpt[i].c_str());//,0+i*0.1);

	//z+=0.1;


	n++;
      }
    
    gr->EndFrame();             // end of the second frame
      
    return gr->GetNumFrame();       // returns the frame number
  }


};


#endif
