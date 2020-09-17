# **************************************************************************
# *
# * Authors: Jorge Garcia Condado (jorgeschool@gmail.com)
# *
# * BCU, Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import pwem

from chimera import Plugin as chimera_plugin

_logo = "icon.jpg"
_references = ['CROLL2018']

class Plugin(pwem.Plugin):

    @classmethod
    def defineBinaries(cls, env):
        """ Install ISOLDE with Chimerax toolshed command """
        from scipion.install.funcs import \
            VOID_TGZ  # Local import to avoid having scipion-app installed when building the package.

        pathToChimera = chimera_plugin.getProgram()
        installPluginsCommand = [("%s --nogui --exit " 
                                  "--cmd 'toolshed install isolde; exit'" % pathToChimera, [])]
        env.addPackage('isolde', version='1.0',
                       tar=VOID_TGZ,
                       default=True,
                       commands=installPluginsCommand)
