/*------------------------------------------------------------------------------
 *                                                                              
 *        Alinea.Nema : Nema Model
 *                                                                              
 *        Copyright 2010 INRA - CIRAD
 *                                                                              
 *        File author(s): Elmer Ccopa Rivera <elmer.ccopa-rivera@avignon.inra.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>         
 *                                                                              
 *        Distributed under the CeCILL-C License.                               
 *        See accompanying file LICENSE.txt 
 *                                                                              
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr                    
 *       
 *        $Id: $
 *                                                                       
 *-----------------------------------------------------------------------------*/

#ifndef __CLASS_WNEMA_UTIL__
#define __CLASS_WNEMA_UTIL__

// Boost.Python Wrapper utilitites

#include <boost/python.hpp>
#include <vector>

template< typename T >
boost::python::list 
convert_vector_to_list( std::vector< T > v )
{
  boost::python::list res;
  int i = 0;
  for( i= 0; i < v.size(); i++ )
    res.append(v[i]); 
  return res;
}

template< typename T >
std::vector< T > 
convert_list_to_vector( boost::python::list l )
{
  int i = 0;
  std::vector<T> res;
  for( i= 0; i < l.attr("__len__")(); i++ )
      res.push_back(boost::python::extract<T>(l[i]));
  return res;
}

template < typename T >
std::vector< std::vector < T > > listlist_vectorvector( boost::python::list ll)
{
  std::vector< std::vector <T> > res;
  int i, j;
  for( i= 0; i <  ll.attr("__len__")(); i++ )
    {
      boost::python::list l = boost::python::extract< boost::python::list > (ll[i]);
      res.push_back(convert_list_to_vector<T>(l));
    }
  return res;
}

template< typename T >
boost::python::list 
convert_vectorofvector_to_listoflist( std::vector< std::vector <T> > v )
{
  boost::python::list res,l;
  std::vector< std::vector <T> > sv;
  int i,j;
  for( i= 0; i < v.size(); i++ ) {
    /*     sv=s[i]; */
    /*     l = boost::python::list(); */
    /*     for( j= 0; j < sv.size(); j++ ) */
    /*       l.append(sv[j]); */
    /*     res.append(l); */
    res.append(convert_vector_to_list <T> (v[i])); 
  }
  return res;
}

#define DEFINE_GETSET(_CLASS,PROP) \
  void Set##_CLASS##PROP( _CLASS & e, list vals ) {\
    e.Set##PROP( convert_list_to_vector <double> ( vals ) );} \
  list Get##_CLASS##PROP( _CLASS & e ) {\
    return convert_vector_to_list <double> (e.Get##PROP() ) ;} \

#define DEFINE_GETVECTOR(_CLASS,PROP) \
  list Get##_CLASS##PROP( _CLASS & e ) {\
    return convert_vector_to_list <double> (e.Get##PROP() ) ;} \

#define DEFINE_GETVECTOROFVECTOR(_CLASS,PROP) \
  list Get##_CLASS##PROP( _CLASS & e ) {\
    return  convert_vectorofvector_to_listoflist <double> (e.Get##PROP() ) ;} \

#define DEFINE_GETSET_PUBLIC(_CLASS,PROP) \
  void Set##_CLASS##PROP( _CLASS & e, list vals ) {\
    e.PROP = convert_list_to_vector <double> ( vals ) ;} \
  list Get##_CLASS##PROP( _CLASS & e ) {\
    return convert_vector_to_list <double> (e.PROP ) ;} \



#define ADD_VECTOR_PROPERTY(_CLASS,PROP,DOC) \
  .add_property(#PROP, Get##_CLASS##PROP, Set##_CLASS##PROP, DOC ) \

#define ADD_VECTOR_OUTPUT(_CLASS,PROP,DOC) \
  .add_property(#PROP, Get##_CLASS##PROP, DOC ) \

#endif
