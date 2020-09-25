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
This module implements a protocol to run ISOLDE from Scipion.
"""

import os

from chimera.constants import CHIMERA_CONFIG_FILE

try:
    from pwem.objects import AtomStruct
except ImportError:
    from pwem.objects import PdbFile as AtomStruct

from pyworkflow.protocol.params import (PointerParam,
                                        BooleanParam)
from pyworkflow.utils.properties import Message

from pwem.protocols import EMProtocol
from pwem.objects import Volume
from pwem.convert.headers import Ccp4Header
from pwem.objects import Transform
from pwem.viewers.viewer_chimera import (Chimera,
                                         sessionFile,
                                         chimeraMapTemplateFileName,
                                         chimeraScriptFileName,
                                         chimeraPdbTemplateFileName)
from chimera import Plugin as chimera

import configparser

class ProtIsolde(EMProtocol):
    """ Protocol to run ISOLDE within Chimera """
    _label = 'isolde operate'

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form, doHelp=True):
        """ Define parts of form consisting mainly of two parts
        1. Input volume and input atomic structure
        2. Flags for simulation
        """
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
        if doHelp:
            form.addSection(label='Help')
            form.addLine('''To save: scipionwrite [model #n] [prefix stringAddedToFilename]
            Reference the pdb and map objects directly not the DataManager
            Type 'help scipionwrite' in chimera command line for details''')

        return form

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runChimeraStep')
        self._insertFunctionStep('createOutputStep')
    # --------------------------- STEPS functions -----------------------------
    def runChimeraStep(self):
        """ Run Chimera script to start simulation and enable scipionwrite
        """
        self.writeChimeraScript()
        args = '--script %s' % (os.path.abspath(self._getExtraPath('scriptChimera.py')))

        # Go to extra dir and save there the output of
        # scipionwrite
        # f.write('cd %s' % os.path.abspath(
        #    self._getExtraPath()))
        # save config file with information
        # this is information is pased from scipion to chimerax
        config = configparser.ConfigParser()
        _chimeraPdbTemplateFileName = \
            os.path.abspath(self._getExtraPath(
                chimeraPdbTemplateFileName))
        _chimeraMapTemplateFileName = \
            os.path.abspath(self._getExtraPath(
                chimeraMapTemplateFileName))
        _sessionFile = os.path.abspath(
            self._getExtraPath(sessionFile))
        protId = self.getObjId()
        config['chimerax'] = {'chimerapdbtemplatefilename':
                                  _chimeraPdbTemplateFileName % protId,
                              'chimeramaptemplatefilename':
                                  _chimeraMapTemplateFileName % protId,
                              'sessionfile': _sessionFile,
                              'enablebundle': True,
                              'protid': self.getObjId()}
                              # set to True when
                              # protocol finished
                              # viewers will check this configuration file
        with open(self._getExtraPath(CHIMERA_CONFIG_FILE),
                  'w') as configfile:
            config.write(configfile)

        cwd = os.path.abspath(self._getExtraPath())
        Chimera.runProgram(chimera.getProgram(), args, cwd=cwd)

    def createOutputStep(self):
        """ Copy the PDB structure and register the output object.
        """

        # Check vol and pdb files
        directory = self._getExtraPath()
        for filename in sorted(os.listdir(directory)):
            if filename.endswith(".mrc"):
                volFileName = os.path.join(directory, filename)
                vol = Volume()
                vol.setFileName(volFileName)

                # fix mrc header
                ccp4header = Ccp4Header(volFileName, readHeader=True)
                sampling = ccp4header.computeSampling()
                origin = Transform()
                shifts = ccp4header.getOrigin()
                origin.setShiftsTuple(shifts)
                vol.setOrigin(origin)
                vol.setSamplingRate(sampling)
                keyword = filename.split(".mrc")[0]
                kwargs = {keyword: vol}
                self._defineOutputs(**kwargs)

            if filename.endswith(".pdb") or filename.endswith(".cif"):
                path = os.path.join(directory, filename)
                pdb = AtomStruct()
                pdb.setFileName(path)
                if filename.endswith(".cif"):
                    keyword = filename.split(".cif")[0].replace(".","_")
                else:
                    keyword = filename.split(".pdb")[0].replace(".", "_")
                kwargs = {keyword: pdb}
                self._defineOutputs(**kwargs)

        # upodate config file flag enablebundle
        # so scipionwrite is disabled
        config = configparser.ConfigParser()
        config.read(self._getExtraPath(CHIMERA_CONFIG_FILE))
        config.set('chimerax', 'enablebundle', 'False')
        with open(self._getExtraPath(CHIMERA_CONFIG_FILE),
                  'w') as configfile:
            config.write(configfile)

    # --------------------------- UTILS functions ----------------------------
    def writeChimeraScript(self):
        """ Builds script file to open the volume and pdb and associate them as
        well as set some further parameters for the simulation.
        """
        f = open(self._getExtraPath('scriptChimera.py'), "w")
        f.write("from chimerax.core.commands import run\n")
        f.write("run(session, 'open %s')\n" % os.path.abspath(self.pdbFileToBeRefined.get().getFileName()))
        f.write("run(session, 'open %s')\n" % os.path.abspath(self.inputVolume.get().getFileName()))
        f.write("run(session, 'clipper assoc #2 to #1')\n")
        f.write("run(session, 'isolde start')\n")
        if self.addH:
            f.write("run(session, 'addh')\n")
        if self.hideHC:
            f.write("run(session, 'hide HC')\n")
        if self.restrainLigands:
            f.write("run(session, 'isolde restrain ligands #1')\n")
        f.close()

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        summary.append("Input Volume provided: %s"
                       % self.inputVolume.get().getFileName())
        summary.append("Input PDB provided: %s"
                       % self.pdbFileToBeRefined.get().getFileName())
        if self.getOutputsSize() > 0:
            directory = self._getExtraPath()
            counter = 1
            summary.append("Produced files:")
            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".pdb") or filename.endswith(".cif"):
                    summary.append("PDB: %s" % filename)
            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".mrc"):
                    summary.append("Map: %s" % filename)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_FILES)
        return summary

    def _citations(self):
        return ['CROLL2018']