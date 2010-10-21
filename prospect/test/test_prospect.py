# -*- python -*-
#
#       Alinea.Prospect: Prospect package
#
#       Copyright 2010 INRA 
#
#       File author(s): Alexis Comar <alexis.comar@avignon.inra.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#

__doc__="""

"""

__license__= "Cecill-C"
__revision__=" $Id $"



from alinea.prospect.prospect import prospect

class TestSceneObject:
    """Test the SceneObject class"""

    def setup_method(self, method):
        self.sceneobj = None

    def teardown_method(self, method):
        self.sceneobj = None

    def test_name(self):
        assert True

