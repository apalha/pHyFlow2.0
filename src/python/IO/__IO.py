__doc__ = """

Input/Output modules for pHyFlow.


Description
-----------
Implements all necessary functions for saving pHyFlow data into the format
of plotting software.

Implemented Formats
-------------------
    VTK files, including PVD and VTU files.


:First added:   2014-02-18
:Last updated:  2014-02-19
:Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
:License:       GNU GPL version 3 or any later version

"""

"""
Reviews:
-------
    (apalha, 2014-02-19) Changed the File class to accept with the operator '<<' objects of the type
                         pHyFlow.vortex.VortexBlobs. Additionally the option to save only the circulation
                         of the vortex blobs only or also the induced velocity at the blobs was added. The
                         time stamp procedure was updated. If the dataObject used in file << dataObject has
                         an attribute dataObject.t then this is used as the time stamp of the PVD files. If
                         not, the internal tStep is used.


"""


__all__ = ['File']

# External packages
import numpy as _numpy
import dolfin as _dolfin
import os as _os
import lxml.etree as _lET
from __evtk.hl import pointsToVTK as _pointsToVTK
from __evtk.hl import linesToVTK as _linesToVTK

# pHyFlow packages
from pHyFlow.blobs import Blobs as _BlobsClass
from pHyFlow.panels import Panels as _PanelsClass
from pHyFlow.eulerian import EulerianSolver as _EulerianSolverClass

class File:
    r"""
    Class used to save all the objects of pHyFlow.

    If the filetype is 'pvd' then the object is saved into VTK file to be read by Paraview, MayaVI,
    VTK or any other software that reads VTK files.

    If the filetype is 'xml' then the object is saved into xml file containing all data that enables the
    reconstruction of the object at the moment it was saved.

    Files saved are not a single file. For example 'pvd' files save the pvd file itself and all the associated
    'vtu' files containing the data at each saves time instant. 'xml' files do the same, the main xml file is
    saved and additional files containing data are also saved.

    Usage
    -----
    .. code-block:: python

        file = File(filename)

    Parameters
    ----------
    filename : string
               Contains the full path to the file where data is to be saved.
               .pvd files save the data to plot using VTK format
               .xml files save the data to be retrieved at another moment

    Attributes
    ----------
    __extension : string
                  The extension of the file to save. Can be '.pvd' or '.xml'
    __ filename : string
                  The filename of the file to save, excluding the extension.
    __filename__path : string
                       The full path of the file to save including the filename. The extension is not included
    __filename_path_full : string
                           The full path of the file to save including the filename and extension.
    __ fileType : string
                  The type of file to save. Can be 'pvd' or 'xml'
    __collection : lxml.etree.SubElement object
                   The data files used in both the '.pvd' and '.xml' file types are contained inside a collection
                   subelement. This is just a xml container for the time stamps of data.
    __root : lxml.etree.Element object
             Both '.pvd. and '.xml' file types are xml files. root is the root of the xml tree structure of the files.
    __tStep : int
              The time step stamp of the saved files. Each time data is saved the main file is updated with a new set of data
              corresponding to the current time step.

    :First Added:   2014-02-18
    :Last Modified: 2014-02-19
    :Copyright:     Copyright (C) 2013 Artur Palha, **pHyFlow**
    :License:       GNU GPL version 3 or any later version

    """
    def __init__(self,filename):
        # save the filename
        self.__filename_path_full = filename

        # extract the full path and base name and the extension
        self.__filename_path,self.__extension = _os.path.splitext(filename)

        # extract the basename of the file
        self.__filename = _os.path.basename(self.__filename_path)

        # check the extension of the file
        if self.__extension == '.pvd':
            self.__fileType = 'pvd'
        else:
            raise ValueError('File must be of type .pvd, it is of type: %s.' % self.__extension)

        # initialize the root and collection variables
        self.__root = None # for now both are None because no object has been associated to the file
        self.__collection =  None
        
        self.__dolfinObjects = None

        # initialize the time step
        self.__tStep = 0


    def __updatePVD(self,dataObject):
        r"""
        Update the '.pvd' file by adding a reference to the datafile containing the corresponding current state of
        the dataObject.'

        Usage
        -----
        .. code-block :: python

            self.__updatePVD(dataObject)

        Parameters
        ----------
        dataObject : pHyFlow.vortex.VortexBlobs
                     The dataObject to be saved to the file.

        Returns
        -------

        Attributes
        ----------

        :First Added:   2014-02-18
        :Last Modified: 2014-03-06
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        # check if the time instant stamp is part of dataObject
        if hasattr(dataObject,'t'): # if it has use the time stamp of the object
            timestep = ('%.6f' % dataObject.t) # get the string representing the current time instant
        else: # if not use the internal time step of the file
            timestep = ('%d' % self.__tStep)

        # determine the type of dataObject and save it accordingly

        if isinstance(dataObject,_BlobsClass): # it is a vortex blob, therefore save vorticity and velocity if option is set
            if self.__root == None:
                # generate the pvd file xml tree
                # the pvd file is an xml file which start with a VTK element and then
                # a collection element
                self.__root = _lET.Element('VTKFile',type='Collection',version='0.1')
                self.__collection = _lET.SubElement(self.__root,'Collection')

            # add the filename of the current time step to the PVD file
            _lET.SubElement(self.__collection,'DataSet',timestep=timestep,group='',part='0',file=self.__filename + ('%09d' % self.__tStep) + '.vtu')

            # generate the xml tree
            pvdTree = _lET.ElementTree(self.__root)
            pvdTree.write(self.__filename_path_full,pretty_print=True,xml_declaration=False,encoding=None)

        elif isinstance(dataObject,_PanelsClass): # it is a panel body, therefore save circulation
            if self.__root == None:
                # generate the pvd file xml tree
                # the pvd file is an xml file which start with a VTK element and then
                # a collection element

                # initialize the root and collection arrays in order to store the root and collection associated to each
                # panel body
                self.__root = []
                self.__collection = []

                # add the vtk file to each of the bodies' pvd files
                for body in range(0,dataObject.nBodies):
                    self.__root.append(_lET.Element('VTKFile',type='Collection',version='0.1'))
                    self.__collection.append(_lET.SubElement(self.__root[body],'Collection'))

            # add the filename of the current time step to the PVD file
            for body,bodyName in enumerate(dataObject.geometryKeys): # loop over all the panel bodies
                _lET.SubElement(self.__collection[body],'DataSet',timestep=timestep,group='',part='0',file=self.__filename + ('_%s_%09d' % (bodyName,self.__tStep)) + '.vtu')

                # generate the xml tree
                pvdTree = _lET.ElementTree(self.__root[body])
                pvdTree.write(self.__filename_path + ('_%s' % bodyName) + self.__extension,pretty_print=True,xml_declaration=False,encoding=None)

        else:
            raise TypeError('Only dataObjects of type pHyFlow.blobs.Blobs or pHyFlow.panels.Panels can be plotted.')


        # update the current time step and time instant
        self.__tStep += 1


    def __saveVTU(self,dataObject):
        r"""
        Save a new '.vtu' file corresponding to the current state of the dataObject.'

        Usage
        -----
        .. code-block :: python

            self.__saveVTU(dataObject)

        Parameters
        ----------
        dataObject : pHyFlow.vortex.VortexBlobs
                     The dataObject to be saved to the file.

        Returns
        -------

        Attributes
        ----------

        :First Added:   2014-02-18
        :Last Modified: 2014-02-19
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """

        # determine the type of dataObject and save it accordingly
        if isinstance(dataObject,_BlobsClass): # it is a vortex blob, therefore save vorticity and velocity if option is set
            if dataObject.plotVelocity == True:
                # since the velocity is to be plotted, compute it first
                vx,vy = dataObject.evaluateVelocity(None,None) # corresponds to evaluating velocity at the blobs
                _pointsToVTK(self.__filename_path+('%09d' % self.__tStep), dataObject.x, dataObject.y, _numpy.zeros(dataObject.x.shape), scalars={"g": dataObject.g}, vectors={"v": (vx,vy,_numpy.zeros(dataObject.g.shape))})
            else: # do not plot the velocity
                _pointsToVTK(self.__filename_path+('%09d' % self.__tStep), dataObject.x, dataObject.y, _numpy.zeros(dataObject.x.shape), scalars={"g": dataObject.g}, vectors=None)

        elif isinstance(dataObject,_PanelsClass): # it is a panel body, therefore save circulation
            for body,bodyName in enumerate(dataObject.geometryKeys): # loop over all the panel bodies
                _linesToVTK(self.__filename_path + ('_%s_%09d' % (bodyName,self.__tStep)), dataObject.xyPanelGlobal[0][body], dataObject.xyPanelGlobal[1][body], _numpy.zeros(dataObject.nPanels[body]+1), scalars={"g": dataObject.sPanel[body]})
        else:
            raise TypeError('Only dataObjects of type pHyFlow.vortex.VortexBlobs can be plotted.')


    def __saveDolfinObject(self, dataObject):
        r"""
        Function to save the dolfin data object
        """
        
        # Make dolfin data file object
        if self.__dolfinObjects is None:
            self.__dolfinObjects = {}
            
            # Make velocity export function
            if dataObject.plotVelocity == True:
                self.__dolfinObjects['velocity'] = _dolfin.File(self.__filename_path + '_velocity.pvd', 'compressed')
            # Make pressure export function                
            if dataObject.plotPressure == True:
                self.__dolfinObjects['pressure'] = _dolfin.File(self.__filename_path + '_pressure.pvd', 'compressed')
            # Make vorticity export function                                
            if dataObject.plotVorticity == True:
                self.__dolfinObjects['vorticity'] = _dolfin.File(self.__filename_path + '_vorticity.pvd', 'compressed')
                
        # Export velocity
        if dataObject.plotVelocity == True:
            self.__dolfinObjects['velocity'] << (dataObject._EulerianSolver__solver.u1, dataObject.t)
        # Export pressure            
        if dataObject.plotPressure == True:
            self.__dolfinObjects['pressure'] << (dataObject._EulerianSolver__solver.p1, dataObject.t)
        # Export vorticity            
        if dataObject.plotVorticity == True:
            self.__dolfinObjects['vorticity'] << (dataObject._EulerianSolver__solver.vorticity(), dataObject.t)
        

    def __lshift__(self,dataObject):
        r"""
        Save the object to the file. If data is already saved, add the new data with the current time stamp.

        Usage
        -----
        .. code-block :: python

            file << dataObject

        Parameters
        ----------
        dataObject : pHyFlow.vortex.VortexBlobs
                     The dataObject to be saved to the file.

        Returns
        -------

        Attributes
        ----------

        :First Added:   2014-02-18
        :Last Modified: 2014-02-19
        :Copyright:     Copyright (C) 2014 Artur Palha, **pHyFlow**
        :License:       GNU GPL version 3 or any later version

        """
        if self.__fileType == 'pvd':
            if isinstance(dataObject, _EulerianSolverClass):
                self.__saveDolfinObject(dataObject)
            else:
                self.__saveVTU(dataObject) # first save the vtu file with the new data
                self.__updatePVD(dataObject) # second update the pvd file containing the information of the new time step