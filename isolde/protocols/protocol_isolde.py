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


"""
Describe your python module here:
This module will provide the traditional Hello world example
"""

from pyworkflow.protocol.params import (PointerParam,
                                        BooleanParam)
from pwem.protocols import EMProtocol
from pwem.viewers.viewer_chimera import Chimera
from chimera import Plugin as chimera

class ProtIsolde(EMProtocol):
    """ Protocol to run ISOLDE within Chimera """
    _label = 'isolde operate'

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form, doHelp=True):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=False,
                      important=True,
                      help="Volume to process")
        form.addParam('pdbFileToBeRefined', PointerParam,
                      pointerClass="AtomStruct", allowsNull=False,
                      important=True,
                      label='Atomic structure',
                      help="PDBx/mmCIF file that you can save after operating "
                           "with it.")
        form.addParam('addH', BooleanParam,
                      label='Add hydrogens to PDB', allowsNull=True,
                      default = True,
                      help="Automatically add hydrogens to PDB")
        form.addParam('hideHC', BooleanParam,
                      label='Hide non-polar hydrogens', allowsNull=True,
                      default=True,
                      help="Automatically hide non-polar hydrogens of PDB")
        form.addParam('restrainLigands', BooleanParam,
                      label='Restrain ligands', allowsNull=True,
                      default=True,
                      help="Automatically restrain ligands in simulation to"
                            " avoid them being sent away flying")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runChimeraStep')
        self._insertFunctionStep('createOutputStep')
    # --------------------------- STEPS functions -----------------------------
    def runChimeraStep(self):
        self.writeChimeraScript()
        args = '--script %s' % (self._getExtraPath('scriptChimera.py'))

        Chimera.runProgram(chimera.getProgram(), args)

    def createOutputStep(self):
        pass

    # --------------------------- UTILS functions ----------------------------
    def writeChimeraScript(self):
        f = open(self._getExtraPath('scriptChimera.py'), "w")
        f.write("from chimerax.core.commands import run\n")
        f.write("run(session, 'open %s')\n" % self.pdbFileToBeRefined.get().getFileName())
        f.write("run(session, 'open %s')\n" % self.inputVolume.get().getFileName())
        f.write("run(session, 'clipper assoc #2 to #1')\n")
        f.write("run(session, 'isolde start')\n")
        if self.addH:
            f.write("run(session, 'addh')\n")
        if self.hideHC:
            f.write("run(session, 'hide HC')\n")
        f.close()

