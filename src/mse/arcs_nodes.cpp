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
#include <float.h> 
#include "arcs_nodes.hxx"

const double NDcomplex_arc::infinity=DBL_MAX;
const double NDcomplex_node::infinity=DBL_MAX;


void NDcomplex_arc::write(std::ofstream &str,std::map<NDcomplex_node *, long> &node2id, std::map<NDcomplex_arc *, long> &arc2id, const arc_it &null_arc)
{
  long i[6];
  
  i[0]=node2id[&(*n[0])];
  i[1]=node2id[&(*n[1])];
  if (next[0]!=null_arc) i[2]=arc2id[&(*next[0])]; else i[2]=-1;
  if (next[1]!=null_arc) i[3]=arc2id[&(*next[1])]; else i[3]=-1;
  if (prev[0]!=null_arc) i[4]=arc2id[&(*prev[0])]; else i[4]=-1;
  if (prev[1]!=null_arc) i[5]=arc2id[&(*prev[1])]; else i[5]=-1;

  str.write((const char *) i, 6*sizeof(long));

}

void NDcomplex_arc::read(std::ifstream &str, std::vector<node_it> &id2node, std::vector<arc_it> &id2arc, const arc_it &null_arc)
{
  long i[6];
  
  str.read((char *) i, 6*sizeof(long));
  
  n[0]=id2node[i[0]];
  n[1]=id2node[i[1]];
  if (i[2]!=-1) next[0] = id2arc[i[2]]; else next[0]=null_arc;
  if (i[3]!=-1) next[1] = id2arc[i[3]]; else next[1]=null_arc;
  if (i[4]!=-1) prev[0] = id2arc[i[4]]; else prev[0]=null_arc;
  if (i[5]!=-1) prev[1] = id2arc[i[5]]; else prev[1]=null_arc;
  
}

void NDcomplex_node::write(std::ofstream &str,std::map<NDcomplex_node *, long> &node2id, std::map<NDcomplex_arc *, long> &arc2id, const arc_it &null_arc, const node_it &null_node)
{
  long i[2];
  long j[2];

  if (p_pair != null_node) i[0]=node2id[&(*p_pair)]; else i[0]=-1;
  if (arc_list != null_arc) i[1]=arc2id[&(*arc_list)]; else i[1]=-1;
  if (ref_up != null_node) j[0]=node2id[&(*ref_up)]; else j[0]=-1;
  if (ref_down != null_node) j[1]=node2id[&(*ref_down)]; else j[1]=-1;
  
  str.write((const char *) &flags, sizeof(char));
  cellID.write(str);
  str.write((const char *) &type, sizeof(char));
  str.write((const char *) &val, sizeof(double));
  str.write((const char *) &val2, sizeof(double));
  str.write((const char *) &numArcs, sizeof(int));
  str.write((const char *) &numUpArcs, sizeof(int));
  str.write((const char *) i, 2*sizeof(long));  
  str.write((const char *) j, 2*sizeof(long));  
}


void NDcomplex_node::read(std::ifstream &str, std::vector<node_it> &id2node, std::vector<arc_it> &id2arc, const arc_it &null_arc, const node_it &null_node)
{
  long i[2];
  long j[2];

  str.read((char *) &flags, sizeof(char));
  cellID.read(str);
  str.read((char *) &type, sizeof(char));
  str.read((char *) &val, sizeof(double));
  str.read((char *) &val2, sizeof(double));
  str.read((char *) &numArcs, sizeof(int));
  str.read((char *) &numUpArcs, sizeof(int));
  str.read((char *) i,2*sizeof(long));
  str.read((char *) j,2*sizeof(long));

  if (i[0]!=-1) p_pair=id2node[i[0]]; else p_pair=null_node;
  if (i[1]!=-1) arc_list=id2arc[i[1]]; else arc_list=null_arc;
  if (j[0]!=-1) ref_up=id2node[j[0]]; else ref_up=null_node;
  if (j[1]!=-1) ref_down=id2node[j[1]]; else ref_down=null_node;
  
}
