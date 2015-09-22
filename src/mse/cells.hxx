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
#ifndef CELL_HXX_
#define CELL_HXX_

#include <iostream>
#include <utility>
#include <cmath>
#include <limits>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std::rel_ops;

template <class T = unsigned long>
class cellIdentity {
public:
  
  static const T TYPE_NBITS = 4;
  static const T TYPE_DEC = (8*sizeof(T)-TYPE_NBITS);
  static const T TYPE_MASK = ( (((T)1)<<TYPE_NBITS) -1 ) << TYPE_DEC;
  static const T ID_MASK = ~TYPE_MASK;
  static const T UNSET =(((T)1)<<TYPE_DEC)-1;
  static const T MAXIMUM = UNSET-1;

private:
    
  T cell_ID;
  T make_cell(T type, T id) 
  {
    return ((type) << TYPE_DEC) + id;
  }

public :
  
  typedef T storageType;
  typedef char typeT;
  typedef long idT; 
  
  typeT type() const {return (unsigned char) ((cell_ID&TYPE_MASK) >> TYPE_DEC);}
  idT id() const {return cell_ID&ID_MASK;}
  bool notSet() const {return (id()==UNSET);}
 
  void set(char type, long id)
  {
    cell_ID = make_cell(type,id);
  }

  void setType(char type)
  {
    cell_ID = make_cell(type,id());
  }

  void setId(char id)
  {
    cell_ID = make_cell(type(),id);
  }

  void set(T val)
  {
    cell_ID = val;
  }

  T get()
  {
    return cell_ID;
  }

  std::pair<typeT,idT> const getAsPair()
  {
    return std::make_pair(type(),id());
  }

  double const getAsDouble()
  {
    if (id()==UNSET) return -std::numeric_limits<double>::max();
    else if (id()==MAXIMUM) return std::numeric_limits<double>::max();
    else return static_cast<double>(type())/10.+static_cast<double>(id());
  }

  cellIdentity(T val)
  {
    cell_ID = val;
  }

  cellIdentity()
  {
    //cell_ID = UNSET;
  }

  cellIdentity(char type, long id)
  {
    //printf("%ld %ld %ld %ld\n",TYPE_NBITS,TYPE_DEC,TYPE_MASK,ID_MASK);
    //printf("%ld %ld ",(unsigned long) type, (unsigned long) id);
    this->set(type,id);
    //printf("-> %ld %ld\n",(unsigned long) this->type(), (unsigned long) this->id());
    
  }

  ~cellIdentity() {}

  bool operator< (const cellIdentity<T> &other) const {return cell_ID < other.cell_ID;}
  bool operator== (const cellIdentity<T> &other) const {return cell_ID == other.cell_ID;}
  T operator() (void) const {return cell_ID;}
  cellIdentity<T> &operator=(T const &val) {this->cell_ID=val; return *this;}

  void write(std::ofstream &str)
  {
    str.write((const char *) &cell_ID, sizeof(T));
  }

  void read(std::ifstream &str)
  {
    str.read((char *) &cell_ID, sizeof(T));
  }
};

#endif
