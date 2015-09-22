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
#ifndef __PD_WIDGET_HXX__
#define __PD_WIDGET_HXX__

#include <QMouseEvent>
#include <QWidget>
#include <QToolBar>
#include <QPainter>
#include <QPixmap>

#include <algorithm>
#include <limits>

#include <mgl/mgl_qt.h>
#include <mgl/mgl_data.h>
#include <mgl/mgl.h>

#include "PD_params.hxx"

#include "NDnetwork.h"
#include "NDnetwork_tags.h"
#include "global.h"
#include "probas.h"

class PD_widget : public QMathGL
{
  Q_OBJECT

  private: 

  double mglDataMin(mglData &dat, bool positive=false, double ret=std::numeric_limits<double>::max())
  {
    long i,N;
    double val;

    if (!positive) val=dat.Minimal();
    else {
      val=dat.Maximal();
      if (val<=0) return ret;

      N=dat.nx*dat.ny*dat.nz;
      for (i=0;i<N;i++)
	{
	  if ((dat.a[i]>0)&&(dat.a[i]<val))
	    {
	      val=dat.a[i];
	    }
	}
    }

    return val;
    
  }

  double mglDataMax(mglData &dat, bool positive=false, double ret=std::numeric_limits<double>::epsilon())
  {
    double val;

    if (!positive) val=dat.Maximal();
    else {
      val=dat.Maximal();
      if (val<=0) val=ret;
    }

    return val;  
  }

  mglData histo(mglData &x, mglData &y, int _Nx=-1, int _Ny=-1)
  {  
    
    int Nx=_Nx;
    int Ny=_Ny;

    if (Nx<=0) Nx=512;
    if (Ny<=0) Ny=(int)((double)height()/(double)width() * Nx);

    long i,j;
    std::vector<double> h(Nx*Ny,0);    
    double x1=(p->logX)?log10(p->x1):p->x1;
    double x2=(p->logX)?log10(p->x2):p->x2;
    double y1=(p->logY)?log10(p->y1):p->y1;
    double y2=(p->logY)?log10(p->y2):p->y2;
    double dx = (x2-x1)/(Nx-1);
    double dy = (y2-y1)/(Ny-1);

    x.Create(Nx);
    y.Create(Ny);

    for (i=0;i<Nx;i++) x.a[i]=x1+i*dx;
    for (i=0;i<Ny;i++) y.a[i]=y1+i*dy;

    for (i=0;i<Nx;i++) x.a[i]=p->x1+i*(p->x2-p->x1)/(Nx-1);
    for (i=0;i<Ny;i++) y.a[i]=p->y1+i*(p->y2-p->y1)/(Ny-1);

    /* because of a bug in mathGL with log coord ... */
    //if (p->logX) for (i=0;i<Nx;i++) x.a[i]=pow(10,x.a[i]);
    //if (p->logY) for (i=0;i<Ny;i++) y.a[i]=pow(10,y.a[i]);

    for (i=0;i<4;i++)
      {
	if (!p->showP[i]) continue;
	
	for (j=0;j<p->x_data[i].nx;j++)
	  {
	    double u=(p->logX)?log10(p->x_data[i].a[j]):p->x_data[i].a[j];
	    double v=(p->logY)?log10(p->y_data[i].a[j]):p->y_data[i].a[j];

	    if ((u!=u)||(v!=v)) continue;

	    double ui = (u-x1)/dx;
	    double vi = (v-y1)/dy;

	    if ((ui<0)||(ui>=Nx)) continue;
	    if ((vi<0)||(vi>=Ny)) continue;

	    h[(long)ui + ((long)vi) * Nx]+=1;
	  }

      }

    mglData ret;
    if (p->logH) for (i=0;i<(long)h.size();i++) h[i]=log10(1.+h[i]);

    //double mn = *std::min_element(h.begin(),h.end());
    //double mx = *std::max_element(h.begin(),h.end());
    


    //for (i=0;i<(long)h.size();i++) h[i]=2.*(h[i]-mn)/(mx-mn)-1;

    ret.Set(&h[0],Nx,Ny);
    ret.Smooth("xy");ret.Smooth("xy");
    for (i=0;i<(long)h.size();i++) if (ret.a[i]<0) ret.a[i]=0;
    
    
    //ret.Set(&h[0],Nx,Ny);
    //printf("ret.Nx=%ld %ld\n",ret.nx,ret.ny);
    return ret;
  }

public:

  bool autoZoom()
  {
    if (p->net==NULL) return false;
    if (p->x_data==NULL) return false;
    if (p->y_data==NULL) return false;

    long i;
    for (i=0;i<4;i++) if ((p->showP[i])&&(p->x_data[i].nx>1)) break;
    long imin=i;
    if (imin==4) return false;

    p->x1=mglDataMin(p->x_data[imin],p->logX,std::numeric_limits<double>::max());
    p->x2=mglDataMax(p->x_data[imin],p->logX,std::numeric_limits<double>::epsilon());
    for (i=imin+1;i<4;i++) 
      {
	if (!p->showP[i]) continue;
	p->x1=std::min(p->x1,mglDataMin(p->x_data[i],p->logX,p->x1));
	p->x2=std::max(p->x2,mglDataMax(p->x_data[i],p->logX,p->x2));
      }
    p->y1=mglDataMin(p->y_data[imin],p->logY,std::numeric_limits<double>::max());
    p->y2=mglDataMax(p->y_data[imin],p->logY,std::numeric_limits<double>::epsilon());
    for (i=imin+1;i<4;i++) 
      {
	if (!p->showP[i]) continue;
	p->y1=std::min(p->y1,mglDataMin(p->y_data[i],p->logY,p->y1));
	p->y2=std::max(p->y2,mglDataMax(p->y_data[i],p->logY,p->y2));
      }
    
    if (p->x1 == std::numeric_limits<double>::max())
      {
	p->x1=0.1;p->x2=100;
	fprintf (stderr,"WARNING: trying to plot strictly negative x_values in log scale !");
      }
    if (p->y1 == std::numeric_limits<double>::max())
      {
	p->y1=0.1;p->y2=100;
	fprintf (stderr,"WARNING: trying to plot strictly negative y_values in log scale !");
      }

    return true;
  }

  bool zoom(double zoomF, mglPoint center)
  {
    if (zoomF==1.) return false;
    if (zoomF<=0) return autoZoom();

    double x=center.x;
    double x1=p->x1;
    double x2=p->x2;
    if (p->logX) 
      {
	if (x<0) x=1;
	if (x1<0) x1=x/10;
	if (x2<0) x2=x1*100;
	x=log10(x);x1=log10(x1);x2=log10(x2);
      }

    if ((x>x1)&&(x<x2))
      {
	double dx = zoomF*(x2-x1);
	double cx = (x - x1) / (x2-x1);	    
	p->x1 = x-cx*dx;
	p->x2 = p->x1 + dx;
	if (p->logX) {p->x1=pow(10,p->x1);p->x2=pow(10,p->x2);}
      }	

    double y=center.y;
    double y1=p->y1;
    double y2=p->y2;
    if (p->logY) 
      {
	if (y<0) y=1;
	if (y1<0) y1=y/10;
	if (y2<0) y2=y1*100;
	y=log10(y);y1=log10(y1);y2=log10(y2);
      }
    if ((y>y1)&&(y<y2))
      {
	double dy=zoomF*(y2-y1);
	double cy = (y - y1) / (y2-y1);
	p->y1 = y-cy*dy;
	p->y2 = p->y1 + dy;
	if (p->logY) {p->y1=pow(10,p->y1);p->y2=pow(10,p->y2);}
      }
    return true;
  }

  bool zoom(double zoomF, QPoint center)
  {
    return zoom(zoomF,XY(center.x(),center.y()));
  }

  public:
  PD_params *p;
 
  PD_widget(QWidget *parent = 0, Qt::WindowFlags f = 0):
    QMathGL(parent,f),
    mouseMode(0),
    pLevel(4,-1),
    sLevel(-1),
    noUpdates(false)
  {
    p=new PD_params();
  }

  ~PD_widget()
  {
    delete p;
  }

  std::pair<bool,bool>  setNetwork(NDnetwork *net)
  {
    std::pair<bool,bool> ret=p->setNetwork(net);
    if (ret.first) 
      {
	plotRatio_slot(false,false);
	autoZoom();  
      }
    //myUpdate();

    return ret;
  }

  QPointF pix2coords(QPoint p)
  {
    mglPoint pt = XY(p.x(),p.y());
    return QPointF(pt.x,pt.y);
  }
  
private:
  //int mouseButtons;
  int old_mouseButtons;
  //QPoint mousePos;
  QPoint old_mousePos;
  int mouseMode;
  
  
protected:

  mglPoint XY(int x, int y)
  {
    mglPoint pt = graph->CalcXYZ(x,y);
    if ((!p->logX)&&(!p->logY)) return pt;

    y=height()-y; // Because Qt and mathGL have different conventions

    mglPoint p1 = graph->CalcScr(mglPoint(p->x1,p->y1));
    mglPoint p2 = graph->CalcScr(mglPoint(p->x2,p->y2));
    double px,py;

    if (p->logX)
      {
	double delta = (x-p1.x);
	delta/=(double)(p2.x-p1.x);
	px=pow(10, log10(p->x1) + delta*(log10(p->x2)-log10(p->x1)));
      }
    else px=pt.x;

    if (p->logY)
       {
	double delta = (y-p1.y);
	delta/=(double)(p2.y-p1.y);
	py=pow(10, log10(p->y1) + delta*(log10(p->y2)-log10(p->y1)));
      }
    else py=pt.y;

    return mglPoint(px,py);
  }

  void mousePressEvent ( QMouseEvent * event )
  {
    if (noUpdates) return;

    int mouseButtons = event->buttons();
    QPoint mousePos = event->pos();

    mglPoint pt = XY(mousePos.x(),mousePos.y());
    emit mousePosChanged(pt.x,pt.y);
    //bool needUpdate=false;

    if (event->modifiers() & Qt::ControlModifier) 
      {
	if (!mouseMode) mouseMode=2;
      }

    if (mouseMode)
      {
	//needUpdate=false;
	old_mousePos=mousePos;
	myUpdate(false);
      }

    old_mousePos=mousePos;
    old_mouseButtons=mouseButtons;
  }

  void mouseReleaseEvent( QMouseEvent * event )
  {
    if (noUpdates) return;

    QPoint mousePos = event->pos();
    int mouseButtons = event->buttons();

    bool needUpdate=false;
    mglPoint pt = XY(mousePos.x(),mousePos.y());

    bool old_L=old_mouseButtons & Qt::LeftButton;
    bool old_R=old_mouseButtons & Qt::RightButton;
    bool old_M=old_mouseButtons & Qt::MidButton;
    bool old_LR=old_L&&old_R;

    bool L=mouseButtons & Qt::LeftButton;
    bool R=mouseButtons & Qt::RightButton;
    bool M=mouseButtons & Qt::MidButton;
    bool LR=L&&R;

    if (mouseMode)
      {
	if ((!L)&&(!R)&&(!M)&&(mouseMode==2)) 
	  mouseMode=0;
	
	needUpdate=false;
	old_mousePos=mousePos;
	//this->repaint();
	myUpdate(false);
      }
    else
      {
	if (old_LR&&(!LR))
	  {
	    autoZoom();
	    needUpdate=true;
	    mouseButtons&=~Qt::LeftButton;
	    mouseButtons&=~Qt::RightButton;
	  }
	else if (old_L&&(!L)) 
	  {
	    mglPoint p1=XY(mousePos.x(),mousePos.y());
	    mglPoint p2=XY(old_mousePos.x(),old_mousePos.y());

	    if (p1==p2)
	      {
		zoom(0.75,pt);
		needUpdate=true;
	      }
	    else
	      {
		p->x1=std::min(p1.x,p2.x);
		p->y1=std::min(p1.y,p2.y);
		p->x2=std::max(p1.x,p2.x);
		p->y2=std::max(p1.y,p2.y);
		needUpdate=true;
	      }
	  }
	else if (old_R&&(!R)) 
	  {zoom(1./0.75,pt);needUpdate=true;}
	
	if (old_M&&(!M)) 
	  {
	    mglPoint p1=XY(mousePos.x(),mousePos.y());
	    mglPoint p2=XY(old_mousePos.x(),old_mousePos.y());
	    
	    if (p->logX)
	      {
		p->x1=pow(10,log10(p->x1) + log10(p2.x)-log10(p1.x));
		p->x2=pow(10,log10(p->x2) + log10(p2.x)-log10(p1.x));
	      }
	    else
	      {
		p->x1-=(p1-p2).x;
		p->x2-=(p1-p2).x;
	      }
	    
	    if (p->logY)
	      {
		p->y1=pow(10,log10(p->y1) + log10(p2.y)-log10(p1.y));
		p->y2=pow(10,log10(p->y2) + log10(p2.y)-log10(p1.y));
	      }
	    else
	      {
		p->y1-=(p1-p2).y;
		p->y2-=(p1-p2).y;
	      }
	    
	    needUpdate=true;
	  }
      }

    old_mousePos=mousePos;
    old_mouseButtons=mouseButtons;
    
    if (needUpdate) 
      {
	//if ((p->logX)||(p->logY)) 
	{
	    if (p->showHisto) p->H=histo(p->Hx,p->Hy);
	}
	myUpdate();
      }
  }

  /*
  bool eventFilter(QObject *obj, QEvent *event)
  {
    if (event->type() == QEvent::MouseMove)
      {
	QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
	QPoint mousePos = mouseEvent->pos();
	QPointF pt=pix2coords(mousePos);
	emit mousePosChanged(pt.x(),pt.y());
      }
    return false;
  }
*/
  
  void mouseMoveEvent (  QMouseEvent * event  )
  {
    if (noUpdates) return;

    QPointF pt=pix2coords(event->pos());
    emit mousePosChanged(pt.x(),pt.y());
    /*
    if (mouseMode) 
      {
	old_mousePos=mousePos;
	
	//this->repaint();
      }
    */
    myUpdate(false);
    //printf("EMIT\n");  
  }
  
signals:
  void mousePosChanged(double x, double y);
  void uncheckPlotRatio();
 
public slots:
  void logX_slot(bool updt=true)
  {
    p->logX=!p->logX;
    if (updt)
      {
	autoZoom();
	if (p->showHisto) p->H=histo(p->Hx,p->Hy);
	//if (p->showHisto) p->histo=histo();
	myUpdate();
      }
  }

  void logY_slot(bool updt=true)
  {
    p->logY=!p->logY;
    if (updt)
      {
	autoZoom();
	if (p->showHisto) p->H=histo(p->Hx,p->Hy);
	//if (p->showHisto) p->histo=histo();
	myUpdate();
      }
  }

  void P0_slot(bool val)
  {
    if (val) p->shownP.push_back(0);
    else 
      {
	std::list<int>::iterator it;
	for (it=p->shownP.begin();it!=p->shownP.end();it++)
	  if (*it==0) break;
	if (it!=p->shownP.end()) p->shownP.erase(it);
      }
    p->showP[0]=val;
    if ((p->showHisto)&&(p->hideP)) p->H=histo(p->Hx,p->Hy);
    //if (p->showHisto) p->histo=histo();
    //autoZoom();
    myUpdate();
  }

  void P1_slot(bool val)
  {
    if (val) p->shownP.push_back(1);
    else 
      {
	std::list<int>::iterator it;
	for (it=p->shownP.begin();it!=p->shownP.end();it++)
	  if (*it==1) break;
	if (it!=p->shownP.end()) p->shownP.erase(it);
      }
    p->showP[1]=val;
    if ((p->showHisto)&&(p->hideP)) p->H=histo(p->Hx,p->Hy);
    //if (p->showHisto) p->histo=histo();
    //autoZoom();
    myUpdate();
  }

  void P2_slot(bool val)
  {
    if (val) p->shownP.push_back(2);
    else 
      {
	std::list<int>::iterator it;
	for (it=p->shownP.begin();it!=p->shownP.end();it++)
	  if (*it==2) break;
	if (it!=p->shownP.end()) p->shownP.erase(it);
      }
    p->showP[2]=val;
    if ((p->showHisto)&&(p->hideP)) p->H=histo(p->Hx,p->Hy);
    myUpdate();
  }

  void showH_slot(bool val,bool updt=true)
  {
    p->showHisto=val;
    if (val) 
      {
	//if ((p->H.nx<=1)||(p->logX)||(p->logY))
	  p->H=histo(p->Hx,p->Hy);
      }
    if (updt) myUpdate();
  }

  void hideP_slot(bool val,bool updt=true)
  {
    p->hideP=val;
    if (updt) myUpdate();
  }

  void logH_slot(bool val)
  {
    p->logH=val;
    p->H=histo(p->Hx,p->Hy);
    if (p->showHisto) myUpdate();
  }

  void recompH_slot()
  {
    p->H=histo(p->Hx,p->Hy);
    if (p->showHisto) myUpdate();
  }

  void setMouseMode_slot(bool val)
  {  
    mouseMode=val;
    myUpdate(false);
  }

  void plotRatio_slot(bool val, bool updt=true)
  {
    long i;
    mglData *old = p->y_data;

    pLevel.assign(4,-1);

    if (val)
      {
	for (i=0;i<4;i++) if (p->Pr[i].nx>1) break;
	if (i!=4) 
	  {
	    p->y_data=p->Pr;
	    p->y_label=std::string("Persistence ratio");
	  }
	else 
	  {
	    p->y_data=p->P;
	    p->y_label=std::string("Persistence");
	    emit uncheckPlotRatio();
	  }
      }
    else
      {
	p->y_data=p->P;
	p->y_label=std::string("Persistence");
      }

    if ((old != p->y_data)&&(updt))
      {
	autoZoom();
	if (p->showHisto) p->H=histo(p->Hx,p->Hy);   
	myUpdate();
      }
  }

private:
  std::vector<double> pLevel;
  double sLevel;
  bool noUpdates;
  
public:
  
  void blockUpdates(bool val=true)
  {
    noUpdates=val;
    if (!val) this->repaint();
  }

  void setNSig(double nsig)
  {
    sLevel=nsig;
    persistenceRatioFromSigma(nsig,p->ndims,&pLevel[0]);
  }

  void setCut(double cut)
  {
    for (int i=0;i<p->ndims;i++) pLevel[i]=cut;
  }


  double getLevels(std::vector<double> &pLevel_)
  {
    pLevel_=pLevel;
    return sLevel;
  }

  void myUpdate(bool all=true)
  {
    
    if (all) update();
    else this->repaint();
  }

  void paintEvent(QPaintEvent * event )
  {   
    static bool old_mouseMode=0;
    if (noUpdates) return;
    if (!mouseMode)
      {
	QPainter painter(this);
	QPen pen(Qt::DashDotLine);
	pen.setWidth(2);
	pen.setColor(Qt::magenta);
	painter.drawPixmap(0,0,pic);
	painter.setPen(pen);
	
	for (long i=0;i<p->ndims;i++)
	  {
	    if (pLevel[i]<0) continue;
	    if (pLevel[i]>p->y2) continue;
	    if (pLevel[i]<p->y1) continue;
	    
	    mglPoint p1 = graph->CalcScr(mglPoint(p->x1,pLevel[i]));
	    mglPoint p2 = graph->CalcScr(mglPoint(p->x2,pLevel[i]));
	    painter.drawLine(QLineF(p1.x,height()-p1.y,p2.x,height()-p2.y));

	    if (i==p->ndims-1)
	      {
		p1.x+=20;p1.y=height()-p1.y-5;
		if (sLevel>=0)
		  painter.drawText(p1.x,p1.y,QString("p=%1-sigma").arg(sLevel));
		else 
		  painter.drawText(p1.x,p1.y,QString("p=%1").arg(pLevel[0]));
	      }	
	  }
	
	if (old_mouseButtons&Qt::LeftButton)
	  {
	    QPointF p1(old_mousePos);
	    QPointF p2(this->mapFromGlobal(QCursor::pos()));
	    //printf("(%g %g) (%g %g)\n",p1.x(),p1.y(),p2.x(),p2.y());
	    pen.setWidth(0);
	    pen.setStyle(Qt::DashLine);
	    pen.setColor(Qt::black);
	    painter.setPen(pen);
	    painter.drawRect(std::min(p1.x(),p2.x()),std::min(p1.y(),p2.y()),
			     fabs(p1.x()-p2.x()),fabs(p1.y()-p2.y()));
	  }
      }
    else 
      {	
	QPointF mousePos(this->mapFromGlobal(QCursor::pos()));

	if ((!old_mouseMode)&&(mouseMode!=2))
	  {
	    //printf("HI\n");
	    if (pLevel[p->ndims-1]<0) 
	      mousePos.setY(-1);
	    else
	      mousePos.setY(height()-graph->CalcScr(mglPoint(p->x1,pLevel[p->ndims-1])).y);
	  }

	QPointF pt(pix2coords(old_mousePos));
	mglPoint p1 = graph->CalcScr(mglPoint(p->x1,p->y1));
	mglPoint p2 = graph->CalcScr(mglPoint(p->x2,p->y2));
	double y=mousePos.y();	

	if (y<p1.y) y=p1.y;
	if (y>p2.y) y=p2.y;

	double level=XY(p1.x,y).y;
	double nsig=-1;
	if (p->y_label == std::string("Persistence ratio"))
	  {
	    if (mousePos.y()>=0)
	      {
		if (level>=1)
		  nsig=sigmaFromPersistenceRatio(level,p->ndims,p->ndims-1);
		else nsig=0;
		if (nsig<0) nsig=0;
		persistenceRatioFromSigma(nsig,p->ndims,&pLevel[0]);
		sLevel=nsig;
	      }
	  }
	else
	  {
	    if (mousePos.y()>=0)
	      {
		if (level<0) level=0;
		for (long i=0;i<p->ndims;i++)  
		  pLevel[i]=level;
		sLevel=-1;
	      }
	  }

	
	//printf("level(%g-sig) :%g %g %g\n",nsig,pLevel[0],pLevel[1],pLevel[2]);

	QPen pen;
	pen.setWidth(2);
	pen.setColor(Qt::magenta);

	QPainter painter(this);
	painter.setPen(pen);
	painter.drawPixmap(0,0,pic);
	//painter.drawLine(QLineF(p1.x,y,p2.x,y));
	for (long i=0;i<p->ndims;i++)
	  {
	    if (pLevel[i]<0) continue;
	    if (pLevel[i]>p->y2) continue;
	    if (pLevel[i]<p->y1) continue;
	    
	    mglPoint pa = graph->CalcScr(mglPoint(p->x1,pLevel[i]));
	    mglPoint pb = graph->CalcScr(mglPoint(p->x2,pLevel[i]));
	    painter.drawLine(QLineF(pa.x,height()-pa.y,pb.x,height()-pb.y));
	  }

	if (pLevel[0]>=0)
	  {
	    mglPoint pa = graph->CalcScr(mglPoint(p->x1,pLevel[p->ndims-1]));
	    pa.x+=20;pa.y=height()-pa.y-5;

	    if (nsig>=0)
	      painter.drawText(pa.x,pa.y,QString("p=%1-sigma").arg(nsig));
	    else 
	      painter.drawText(pa.x,pa.y,QString("p=%1").arg(pLevel[0]));
	  }	
      }
  
     old_mouseMode=mouseMode;
  }
  
};

#endif
