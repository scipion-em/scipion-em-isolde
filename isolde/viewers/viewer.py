# **************************************************************************
# *
# * Authors: Jorge Garcia Condado (jgcondado@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os

from ..protocols.protocol_isolde import ProtIsolde

from pwem.convert import Ccp4Header
from pwem.emlib.image import ImageHandler
from pwem.objects import Volume
from pwem.objects import Transform

from pwem.viewers.viewer_chimera import (Chimera,
                                         sessionFile)
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer


class IsoldeViewer(Viewer):
    """ Visualize the output of protocols protocol_fit and protocol_operate """
    _environments = [DESKTOP_TKINTER]
    _targets = [ProtIsolde]

    def _visualize(self, obj, **args):
        """ Visualize any saved pdb and map if none were saved
        show the input files.
        """
        _inputVolFlag = False
        _inputPDBFlag = False
        directory = self.protocol._getExtraPath()

        fnCmd = self.protocol._getTmpPath("chimera_output.cxc")
        f = open(fnCmd, 'w')
        f.write('cd %s\n' % os.getcwd())

        counter = 0
        # Find all saved maps and pdbs from protocl. If none
        # are found show the input files to the protocol
        for filename in sorted(os.listdir(directory)):
            if filename.endswith(".mrc"):
                _inputVolFlag = True
                counter += 1
                volFileName = os.path.join(directory, filename)
                vol = Volume()
                vol.setFileName(volFileName)

                # fix mrc header
                ccp4header = Ccp4Header(volFileName, readHeader=True)
                sampling = ccp4header.computeSampling()
                origin = Transform()
                shifts = ccp4header.getOrigin()
                origin.setShiftsTuple(shifts)
                f.write("open %s\n" % volFileName)
                f.write("volume #%d style surface voxelSize %f\n"
                        "volume #%d origin %0.2f,%0.2f,%0.2f\n"
                        % (counter, sampling, counter, shifts[0], shifts[1], shifts[2]))
                f.write("volume #%d level %0.3f\n"
                        % (counter, 0.001))
                # Set volume to translucent
                f.write("volume #%d transparency 0.5\n" % counter)

        for filename in os.listdir(directory):
            if filename.endswith(".pdb") or filename.endswith(".cif"):
                _inputPDBFlag = True
                path = os.path.join(directory, filename)
                f.write("open %s\n" % path)

   #     if !_inputVolFlag:
   #         pass

        f.close()

        # run in the background
        Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")
        return []