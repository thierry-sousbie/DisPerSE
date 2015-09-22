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
#ifndef __PDV_INTERFACE_HXX__
#define __PDV_INTERFACE_HXX__

#include <QMouseEvent>
#include <QApplication>
#include <QWidget>
#include <QLayout>
#include <QScrollArea>
#include <QToolBar>
#include <QAction>
#include <QToolButton>
#include <QPushButton>
#include <QSpacerItem>
#include <QLabel>
#include <QMainWindow>
#include <QDialog>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGraphicsTextItem>
#include <QTextStream>
#include <QFile>
#include <QString>
#include <QIODevice>

#include "PD_widget.hxx"
#include "PD_draw.hxx"

#include "NDnetwork.h"
#include "NDnet_interface.hxx"

class PDV_interface : public QWidget
{
  Q_OBJECT
  
  private:

  class myHelpDialog : public QDialog
  {
  public:
    myHelpDialog(QWidget * parent = 0, Qt::WindowFlags f = 0):QDialog(parent,f)
    {
      QPushButton *b=new QPushButton(this);     
      QGraphicsScene *s = new QGraphicsScene(this);
      QGraphicsView *v =new QGraphicsView(s,this);
      QVBoxLayout *L=new QVBoxLayout(this);

      b->setText("close");
      L->addWidget(v);
      L->addWidget(b);
      connect(b,SIGNAL(pressed()),SLOT(hide()));

      QString str;
      QFile data(":/help.html");
      if (!data.open(QIODevice::ReadOnly | QIODevice::Text))
        str=QString("ERROR READING RESSOURCE FILE !!!");
      else
	{
	  str=QTextStream(&data).readAll();
	}
    
      QGraphicsTextItem *t=new QGraphicsTextItem();
      /*
      QString str("<p><b><u>PDview v1.0</u></p></b>");
      
      str+=QString("Left click:\t Zoom in\n");
      str+=QString("Right click:\t Zoom out\n");
      */
      t->setHtml(str);
      s->addItem(t);
      
    }

  };

  const QString progName;
  PD_widget *graphW;

  QLabel *mousePosLabel;
  QAction *plotRatio;
  QAction *hideP;
  QAction *showH;
  QAction *logX;
  QAction *logY;
  QAction *help;
  
  myHelpDialog *helpDialog;
 
  public:

  bool loadNetwork(const char *fname, bool useRatio=false)
  {
    NDnetwork *net=ndnet::IO::load(std::string(fname));
    if (net==NULL) return false;
    
    graphW->blockUpdates();

    std::pair<bool,bool> ret=graphW->setNetwork(net);
    if (ret.first)
      {
	bool noUpdt=false;
	setWindowTitle(progName+QString(" (%1 pairs in \"").arg(net->nfaces[1])+QString(fname)+QString("\")"));
	if (net->nfaces[1]>20000)
	  {
	    hideP->setChecked(true);	    
	    showH->setChecked(true);
	    noUpdt=true;
	  }
	if (useRatio) 
	  {
	    plotRatio->setChecked(useRatio);
	    graphW->logY_slot();
	    graphW->logX_slot();
	    noUpdt=true;
	  }
	if (!noUpdt) graphW->myUpdate();
      }
    else graphW->myUpdate();

    graphW->blockUpdates(false);

    return true;
  }

  void setNSig(double nsig)
  {
    graphW->setNSig(nsig);
  }

  void setCut(double cut)
  {
    graphW->setCut(cut);
  }

  PDV_interface(QWidget *parent=0): QWidget(parent), 
				    progName("Persistence Diagram Viewer")
  {
    regularQuit=false;
    graphW=new PD_widget(this);
    setWindowTitle(progName);
    
    QToolBar *toolbar = new QToolBar(this);
    //QToolButton *button = new QToolButton(toolbar);
    //button->setIcon(QIcon("/usr/share/easycrypt/icons/create.png"));
    //QToolBar *toolbar = addToolBar("main toolbar");
    
    
    logX =toolbar->addAction("lg(X)");
    logX->setToolTip("Linear / logarithmic coordinates for X axis");
    logY =toolbar->addAction("lg(Y)");
    logY->setToolTip("Linear / logarithmic coordinates for Y axis");
    
    toolbar->addSeparator();
    
    hideP =toolbar->addAction("Hide");
    hideP->setToolTip("Show/hide marks");
    hideP->setCheckable(true);
    hideP->setChecked(false);   
    QAction *P0=toolbar->addAction("0");
    P0->setToolTip("Show/hide type 0 persistence pairs (minima - saddle)");
    P0->setCheckable(true);
    P0->setChecked(true);
    QAction *P1=toolbar->addAction("1");
    P1->setToolTip("Show/hide type 1 persistence pairs");
    P1->setCheckable(true);
    P1->setChecked(true);
    QAction *P2=toolbar->addAction("2");
    P2->setToolTip("Show/hide type 2 persistence pairs");
    P2->setCheckable(true);
    P2->setChecked(true);

    toolbar->addSeparator();

    showH=toolbar->addAction("Histo");
    showH->setToolTip("Display the persistence diagram as an histogram");
    showH->setCheckable(true);
    showH->setChecked(false);
    QAction *logH=toolbar->addAction("Log");
    logH->setToolTip("Enable/disable log scale for the histogram color scale");
    logH->setCheckable(true);
    logH->setChecked(true);

    
    QAction *recompH=toolbar->addAction("Updt");
    recompH->setToolTip("Force updating the histogram");    
    

    toolbar->addSeparator();

    plotRatio=toolbar->addAction("Ratio");
    plotRatio->setToolTip("Use this mode for delaunay tesselations with DTFE density.\ne.g. if network was computed from a point sample with 'delaunay_nD'.\nPlots persistence ratios **if available**.\nAlso switches to 'nsig' type thresholds.");
    plotRatio->setCheckable(true);
    plotRatio->setChecked(false);
    /*
    QAction *unlink=toolbar->addAction("Unlink");
    unlink->setToolTip("Link/unlink persistence thresholds.\nWhen unlinked, different thresholds can be set for different pair types.\nnb: it is usually safer to keep them linked");
    unlink->setCheckable(true);
    unlink->setChecked(false);
    */
    QAction *setPer=toolbar->addAction("Set");
    setPer->setToolTip("Shortcut : Ctrl + click\n Enable setting persistence threshold using the mouse.\n nb: threshold is expressed in 'number of sigmas' when persistence ratio is used");// (leftB: 0/midB: 1/rightB: 2)
    setPer->setCheckable(true);
    setPer->setChecked(false);

    toolbar->addSeparator();

    QAction *allDone=toolbar->addAction("DONE");
    allDone->setToolTip("Quits the application returning persistence thresholds");
    QFont ft=allDone->font();
    ft.setBold(true);
    allDone->setFont(ft);
    toolbar->addSeparator();
    help=toolbar->addAction("?");
    

    helpDialog = new myHelpDialog(this);
       
    //toolbar->addSeparator();
    
    //QAction *help=toolbar->addAction("?");
    
    connect(help, SIGNAL(triggered()), SLOT(helpPopup()));
    
    connect(hideP, SIGNAL(toggled(bool)), graphW, SLOT(hideP_slot(bool)));

    connect(logX, SIGNAL(triggered()), graphW, SLOT(logX_slot()));
    connect(logY, SIGNAL(triggered()), graphW, SLOT(logY_slot()));

    connect(P0, SIGNAL(toggled(bool)), graphW, SLOT(P0_slot(bool)));
    connect(P1, SIGNAL(toggled(bool)), graphW, SLOT(P1_slot(bool)));
    connect(P2, SIGNAL(toggled(bool)), graphW, SLOT(P2_slot(bool)));

    connect(showH, SIGNAL(toggled(bool)), graphW, SLOT(showH_slot(bool)));
    connect(logH, SIGNAL(toggled(bool)), graphW, SLOT(logH_slot(bool)));
    connect(recompH, SIGNAL(triggered()), graphW, SLOT(recompH_slot()));

    connect(plotRatio,SIGNAL(toggled(bool)), graphW, SLOT(plotRatio_slot(bool)));
    connect(setPer,SIGNAL(toggled(bool)), graphW, SLOT(setMouseMode_slot(bool)));


    connect(graphW,SIGNAL(uncheckPlotRatio()), this, SLOT(uncheckPlotRatio()));
    connect(graphW, SIGNAL(mousePosChanged(double,double)), this, SLOT(updateMousePosLabel(double,double)));

    connect(allDone,SIGNAL(triggered()),this,SLOT(onQuit()));
   
    mousePosLabel = new QLabel("x: 0, y: 0",this);
   
    QVBoxLayout *viewLayout = new QVBoxLayout(this);
    //QWidget *barWidget = new QWidget(this);
    QHBoxLayout *barLayout = new QHBoxLayout();
    viewLayout->addLayout(barLayout);  
    //viewLayout->addWidget(barWidget);  
    viewLayout->addWidget(graphW);

    mousePosLabel->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Fixed);
    toolbar->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Minimum);    
    //mousePosLabel->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Minimum);    
    //barWidget->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Maximum);
    
    graphW->setSizePolicy(QSizePolicy::Ignored,QSizePolicy::Ignored);
    
    barLayout->addWidget(toolbar);
    barLayout->addSpacerItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Minimum));
    barLayout->addWidget(mousePosLabel);

    graphW->autoResize=true;         
    graphW->setDraw(PD_draw::Draw,graphW->p);    
  }

  ~PDV_interface()
  {

  }
private:
  bool regularQuit;

private slots:

  void closeEvent(QCloseEvent *event)
  {
    if (!regularQuit) printf("Selected threshold: none\n");
  }
  
  void onQuit()
  {
    std::vector<double> level;
    double nsig=graphW->getLevels(level);

    if (level[0]<=0)
      {
	printf("Selected threshold: none\n");
      }
    else
      {
	printf("Selected threshold:");
	if (nsig>=0) 
	  {
	    printf(" -nsig %g\n",nsig);
	    
	    //printf("-cutR [%g",level[0]);
	    //for (unsigned int i=1;i<level.size();i++)
	    //  if (level[i]>=0) printf(" %g",level[i]);
	    //printf("]\n");
	    
	  }
	else printf(" -cut %g\n",level[0]);
      }
    regularQuit=true;
    close();
  }

  void updateMousePosLabel(double x, double y)
  {
    static char txt[255];
    sprintf(txt,"x=%5.5g, y=%5.5g",x,y);
    mousePosLabel->setText(QString(txt));
    
  }

  void uncheckPlotRatio()
  {
    plotRatio->setChecked(false);
  }

  void helpPopup()
  {
    helpDialog->show();
  }
  
};

#endif
