# ***********************************************************************************
# * Copyright 2010 - 2012 Paulo A. Herrera. All rights reserved.                    * 
# *                                                                                 *
# * Redistribution and use in source and binary forms, with or without              *
# * modification, are permitted provided that the following conditions are met:     *
# *                                                                                 *
# *  1. Redistributions of source code must retain the above copyright notice,      *
# *  this list of conditions and the following disclaimer.                          *
# *                                                                                 *
# *  2. Redistributions in binary form must reproduce the above copyright notice,   *
# *  this list of conditions and the following disclaimer in the documentation      *
# *  and/or other materials provided with the distribution.                         *
# *                                                                                 *
# * THIS SOFTWARE IS PROVIDED BY PAULO A. HERRERA ``AS IS'' AND ANY EXPRESS OR      *
# * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF    *
# * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO      *
# * EVENT SHALL <COPYRIGHT HOLDER> OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,        *
# * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,  *
# * BUT NOT LIMITED TO, PROCUREMEN OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,    *
# * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY           *
# * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING  *
# * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS              *
# * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    *
# ***********************************************************************************

# **************************************
# *  High level Python library to      *
# *  export data to binary VTK file.   *
# **************************************

from vtk import * # VtkFile, VtkUnstructuredGrid, etc.
import numpy as np

# =================================
#       Helper functions      
# =================================
def _addDataToFile(vtkFile, cellData=None, pointData=None):
    # Point data
    if pointData <> None:
        # the scalars and vectors have to be separated because the
        # data has to be open for both scalars and vectors. Otherwise,
        # at the second time vtkFile.openData is invoked the previous
        # data is deleted
        if pointData['scalars'] <> None:
            keysPointScalars = pointData['scalars'].keys()
        else:
                keysPointScalars = [None]
        
        if pointData['vectors'] <> None:
            keysPointVectors = pointData['vectors'].keys()
        else:
            keysPointVectors = [None]
		
        vtkFile.openData("Point", scalars = keysPointScalars[0], vectors = keysPointVectors[0])
        
        if pointData['scalars'] <> None:
            for key in keysPointScalars:
                data = pointData['scalars'][key]
                vtkFile.addData(key, data)
            
	    
        if pointData['vectors'] <> None:
            for key in keysPointVectors:
                data = pointData['vectors'][key]
                vtkFile.addData(key, data)

        vtkFile.closeData("Point")

    # Cell data
    if cellData <> None:
        # the scalars and vectors have to be separated because the
        # data has to be open for both scalars and vectors. Otherwise,
        # at the second time vtkFile.openData is invoked the previous
        # data is deleted
        if cellData['scalars'] <> None:
            keysCellsScalars = cellData['scalars'].keys()
        else:
            keysCellsScalars = [None]

        if cellData['vectors'] <> None:
            keysCellsVectors = cellData['vectors'].keys()
        else:
            keysCellsVectors = [None]

        vtkFile.openData("Cell", scalars = keysCellsScalars[0], vectors = keysCellsVectors[0])

        if cellData['scalars'] <> None:
            for key in keysCellsScalars:
                data = cellData['scalars'][key]
                vtkFile.addData(key, data)


        if cellData['vectors'] <> None:
            for key in keysCellsVectors:
                data = cellData['vectors'][key]
                vtkFile.addData(key, data)

        vtkFile.closeData("Cell")

def _appendDataToFile(vtkFile, cellData=None, pointData=None):
    # Append data to binary section
    if pointData <> None:
        if pointData['scalars'] <> None:
            keys = pointData['scalars'].keys()
            for key in keys:
                data = pointData['scalars'][key]
                vtkFile.appendData(data)
        
        if pointData['vectors'] <> None:
            keys = pointData['vectors'].keys()
            for key in keys:
                data = pointData['vectors'][key]
                vtkFile.appendData(data)
        
    if cellData <> None:
        if cellData['scalars'] <> None:
            keys = cellData['scalars'].keys()
            for key in keys:
                data = cellData['scalars'][key]
                vtkFile.appendData(data)

        if cellData['vectors'] <> None:
            keys = cellData['vectors'].keys()
            for key in keys:
                data = cellData['vectors'][key]
                vtkFile.appendData(data)

def _requiresLargeVTKFileSize(vtkFileType, numPoints, numCells, pointData=None, cellData=None):
    """
        If the block size will be close to or larger than what a uint32 can hold we
        need to use VTK's large file ability. note that this is only
        available in VTK 6.0 and greater.
        Contributed by Andy Bauer.

        PARAMETERS:
            vtkFileType: name of the type of grid. If it is
                         "VtkUnstructuredGrid" or VtkStructuredGrid"
                         then it will add extra space needed
                         for storing points. If it is "VtkUnstructuredGrid" it wil
                         also add extra space for cell connectivities.
            numPoints: the number of points in the data set. It is only used
                       if vtkFileType is "VtkUnstructuredGrid" or VtkStructuredGrid".
            numCells: the number of cells in the data set. It is only used
                      if vtkFileType is "VtkUnstructuredGrid".
            cellData: dictionary containing arrays with cell centered data.
                      Keys should be the names of the data arrays.
                      Arrays must have the same dimensions in all directions and must contain
                      only scalar data.
            pointData: dictionary containing arrays with node centered data.
                       Keys should be the names of the data arrays.
                       Arrays must have same dimension in each direction and
                       they should be equal to the dimensions of the cell data plus one and
                       must contain only scalar or vector data. Scalar data is given by the
                       key 'scalars' and vector data is given by the key 'vectors'.
	    
        RETURNS:
            True if the large file format should be used.

    """
    sum = 0
    if pointData <> None:
        if pointData['scalars'] <> None:
            # do the counting for scalar data
            keys = pointData['scalars'].keys()
            for key in keys:
                data = pointData['scalars'][key]
                sum = sum + data.size*data.dtype.itemsize

        if pointData['vectors'] <> None:
            # do the counting for vector data
            keys = pointData['vectors'].keys()
            for key in keys:
                data = pointData['vectors'][key]
                for vectorComponent in range(len(data)): # loop over the components of the vectors
                    sum = sum + data[vectorComponent].size*data[vectorComponent].dtype.itemsize

    if cellData <> None:
        # keys = cellData.keys()
        # for key in keys:
        #     data = cellData[key]
        #     sum = sum + data.size*data.dtype.itemsize


        if cellData['scalars'] <> None:
            # do the counting for scalar data
            keys = cellData['scalars'].keys()
            for key in keys:
                data = cellData['scalars'][key]
                sum = sum + data.size*data.dtype.itemsize

        if cellData['vectors'] <> None:
            # do the counting for vector data
            keys = cellData['vectors'].keys()
            for key in keys:
                data = cellData['vectors'][key]
                for vectorComponent in range(len(data)): # loop over the components of the vectors
                    sum = sum + data[vectorComponent].size*data[vectorComponent].dtype.itemsize
   
    if vtkFileType == "VtkUnstructuredGrid":
        sum = sum + numPoints*3*8 + numCells*8*8
    elif vtkFileType == "VtkStructuredGrid":
        sum = sum + numPoints*3*8

    # leave a little wiggle room instead of using 4294967295
    return sum > 4000000000

# =================================
#       High level functions      
# =================================
def imageToVTK(path, origin = (0.0,0.0,0.0), spacing = (1.0,1.0,1.0), cellData = None, pointData = None ):
    """ Exports data values as a rectangular image.
        
        PARAMETERS:
            path: name of the file without extension where data should be saved.
            origin: grid origin (default = (0,0,0))
            spacing: grid spacing (default = (1,1,1))
            cellData: dictionary containing arrays with cell centered data.
                      Keys should be the names of the data arrays.
                      Arrays must have the same dimensions in all directions and must contain 
                      only scalar data.
            nodeData: dictionary containing arrays with node centered data.
                      Keys should be the names of the data arrays.
                      Arrays must have same dimension in each direction and 
                      they should be equal to the dimensions of the cell data plus one and
                      must contain only scalar data.
         
         RETURNS:
            Full path to saved file.

        NOTE: At least, cellData or pointData must be present to infer the dimensions of the image.
    """
    assert (cellData <> None or pointData <> None)
    
    # Extract dimensions
    start = (0,0,0)
    end = None
    if cellData <> None:
        keys = cellData.keys()
        data = cellData[keys[0]]
        end = data.shape
    elif pointData <> None:
        keys = pointData.keys()
        data = pointData[keys[0]]
        end = data.shape
        end = (end[0] - 1, end[1] - 1, end[2] - 1)

    largeFileFlag = _requiresLargeVTKFileSize("VtkImageData", numPoints = 0, numCells = 0, pointData = pointData, cellData = cellData)

    # Write data to file
    w = VtkFile(path, VtkImageData, largeFileFlag)
    w.openGrid(start = start, end = end, origin = origin, spacing = spacing)
    w.openPiece(start = start, end = end)
    _addDataToFile(w, cellData, pointData)
    w.closePiece()
    w.closeGrid()
    _appendDataToFile(w, cellData, pointData)
    w.save()
    return w.getFileName()

# ==============================================================================
def gridToVTK(path, x, y, z, cellData = None, pointData = None):
    """
        Writes data values as a rectilinear or rectangular grid.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: coordinates of the nodes of the grid. They can be 1D or 3D depending if
                     the grid should be saved as a rectilinear or logically structured grid, respectively.
                     Arrays should contain coordinates of the nodes of the grid.
                     If arrays are 1D, then the grid should be Cartesian, i.e. faces in all cells are orthogonal.
                     If arrays are 3D, then the grid should be logically structured with hexahedral cells.
                     In both cases the arrays dimensions should be equal to the number of nodes of the grid.
            cellData: dictionary containing arrays with cell centered data.
                      Keys should be the names of the data arrays.
                      Arrays must have the same dimensions in all directions and must contain 
                      only scalar data.
            pointData: dictionary containing arrays with node centered data.
                       Keys should be the names of the data arrays.
                       Arrays must have same dimension in each direction and 
                       they should be equal to the dimensions of the cell data plus one and
                       must contain only scalar data.

        RETURNS:
            Full path to saved file.

    """
    # Extract dimensions
    start = (0,0,0)
    nx = ny = nz = 0

    if (x.ndim == 1 and y.ndim == 1 and z.ndim == 1):
        nx, ny, nz = x.size - 1, y.size - 1, z.size - 1
        isRect = True
        ftype = VtkRectilinearGrid
        largeFileFlag = _requiresLargeVTKFileSize("VtkRectilinearGrid", 0, 0, pointData, cellData)
    elif (x.ndim == 3 and y.ndim == 3 and z.ndim == 3):
        s = x.shape
        nx, ny, nz = s[0] - 1, s[1] - 1, s[2] - 1
        isRect = False
        ftype = VtkStructuredGrid
        largeFileFlag = _requiresLargeVTKFileSize("VtkStructuredGrid", numPoints = x.size, numCells = 0, pointData = pointData, cellData = cellData)
    else:
        assert(False)
    end = (nx, ny, nz)


    w =  VtkFile(path, ftype, largeFileFlag)
    w.openGrid(start = start, end = end)
    w.openPiece(start = start, end = end)

    if isRect:
        w.openElement("Coordinates")
        w.addData("x_coordinates", x)
        w.addData("y_coordinates", y)
        w.addData("z_coordinates", z)
        w.closeElement("Coordinates")
    else:
        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")

    _addDataToFile(w, cellData, pointData)
    w.closePiece()
    w.closeGrid()
    # Write coordinates
    if isRect:
        w.appendData(x).appendData(y).appendData(z)
    else:
        w.appendData( (x,y,z) )
    # Write data
    _appendDataToFile(w, cellData, pointData)
    w.save()
    return w.getFileName()


# ==============================================================================
def pointsToVTK(path, x, y, z, scalars=None,vectors=None):
    """
        Export points and associated data as an unstructured grid.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: 1D arrays with coordinates of the points.
            scalars: dictionary with scalar variables associated to each point.
                     Keys should be the names of the variable stored in each array.
                     All arrays must have the same number of elements and equal to the number
                     of points.
	        vectors: dictionary with vector variables associated to each point.
                     Keys should be the names of the variable stored in each array.
                     All arrays must have the same number of elements and equal to the number
                     of points. The vector variables are tuples of numpy.array. Can be 2D or 3D.

        RETURNS:
            Full path to saved file.

    """
    assert (x.size == y.size == z.size)
    npoints = x.size
    
    # create some temporary arrays to write grid topology
    offsets = np.arange(start = 1, stop = npoints + 1, dtype = 'int32')   # index of last node in each cell
    connectivity = np.arange(npoints, dtype = 'int32')                    # each point is only connected to itself
    cell_types = np.empty(npoints, dtype = 'uint8') 
   
    cell_types[:] = VtkVertex.tid
	
    largeFileFlag = _requiresLargeVTKFileSize("VtkUnstructuredGrid", numPoints = npoints, numCells = npoints, pointData = {'scalars':scalars,'vectors':vectors}, cellData = None)

    w = VtkFile(path, VtkUnstructuredGrid, largeFileFlag)
    w.openGrid()
    w.openPiece(ncells = npoints, npoints = npoints)
    
    w.openElement("Points")
    w.addData("points", (x,y,z))
    w.closeElement("Points")
    w.openElement("Cells")
    w.addData("connectivity", connectivity)
    w.addData("offsets", offsets)
    w.addData("types", cell_types)
    w.closeElement("Cells")
    
    _addDataToFile(w, cellData = None, pointData = {'scalars':scalars,'vectors':vectors})

    w.closePiece()
    w.closeGrid()
    w.appendData( (x,y,z) )
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)

    _appendDataToFile(w, cellData = None, pointData = {'scalars':scalars,'vectors':vectors})

    w.save()
    return w.getFileName()


def linesToVTK(path, x, y, z, scalars=None):
    """
        Export line and associated data as an unstructured grid.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: 1D arrays with coordinates of the points that make up the line.
                     The points given go along the line.
            scalars: dictionary with scalar variables associated to each segment
                     of the line.
                     Keys should be the names of the variable stored in each array.
                     All arrays must have the same number of elements and equal to the number
                     of panels (number of points - 1).

        RETURNS:
            Full path to saved file.

    """

    assert(x.size == y.size == z.size)
    nPoints = x.size
    nLines = nPoints-1

    # create some temporary arrays to write the topology of the grid

    # generate the connectivity of the grid
    # this array specifies the points in each line segment (each cell)
    # this is done in the following way considering the line:
    #  P1      P2      P3      P4
    #   X-------X-------X-------X
    #
    #      L1      L2      L3
    #
    # connectivity will then be:
    #   [1 2 2 3 3 4], since cell L1 has points 1 and 2, cell L2 chas points 2 and 3, etc.
    connectivity = np.zeros(2*nLines,dtype = 'int32')
    connectivity[0::2] = np.arange(0,nPoints-1)
    connectivity[1::2] = np.arange(1,nPoints)

    # generate the offsets of the grid
    # this array specifies where the start of each cell is in the connectivity array
    # taking the previous example we have that cell 1 starts at 0 (this is not
    # specified because cell 1 always starts at 0), cell 2 starts at 2, cell 3
    # starts at 4 and then the last number is the size of the connectivity grid + 1
    offsets = np.arange(2,2*nLines+1,2,dtype = 'int32')

    # generate the array of cell types
    cell_types = np.empty(nLines, dtype = 'uint8')
    cell_types[:] = VtkLine.tid

    # determine if a large file should be generated or not
    largeFileFlag = _requiresLargeVTKFileSize("VtkUnstructuredGrid", numPoints = nPoints, numCells = nLines, pointData = None, cellData = {'scalars':scalars,'vectors':None})

    # open the vtk file to start adding the data
    # the line is a vtk unstructure grid with two points
    w = VtkFile(path, VtkUnstructuredGrid, largeFileFlag)
    # open the grid to add the mesh
    w.openGrid()
    w.openPiece(ncells = nLines, npoints = nPoints)

    # add the meta data of the points
    w.openElement('Points')
    w.addData('points',(x,y,z))
    w.closeElement("Points")

    # add the meta data of the cells
    w.openElement("Cells")
    w.addData("connectivity", connectivity)
    w.addData("offsets", offsets)
    w.addData("types", cell_types)
    w.closeElement("Cells")

    # add the metada of the cell data
    _addDataToFile(w, cellData = {'scalars': scalars, 'vectors': None}, pointData = None)

    # finish metadata input
    w.closePiece()
    w.closeGrid()

    # add the points data to the end of the xml file
    w.appendData( (x,y,z) )
    # add the connectivity, offsets and cell type data to the end of the xml file
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)
    # add the actual data to the end of the xml file
    _appendDataToFile(w,cellData={'scalars': scalars, 'vectors': None}, pointData = None)

    # save the vtk xml file
    w.save()

