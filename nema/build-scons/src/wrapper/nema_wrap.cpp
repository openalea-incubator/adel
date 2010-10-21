/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *                        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: stat_tool_wrap.cpp 9099 2010-06-08 09:03:00Z pradal $
 *
 *-----------------------------------------------------------------------------*/



/* WRAPPER Boost.python for stat_tool class */
#include "export_plant.h"

#include <boost/python.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 103400
#include <boost/python/docstring_options.hpp>
#endif

using namespace boost::python;

// Define python module "_stat_tool"
BOOST_PYTHON_MODULE(_nema)
{
  //show_user_defined : true
  //show_signatures : false
#if BOOST_VERSION >= 103400
  docstring_options doc_options(true, false);
#endif

    class_grain();
    class_root();
    class_environment();
    class_entity();

    class_plant();
}

