#!/usr/bin/env python

import colorcet as cc
import numpy as np
import os
import pyvista as pv
import pyvistaqt as pvqt
import sys
import vtk

os.environ["QT_API"] = "pyqt5"

from enum import Enum
from qtpy import QtWidgets

FAR, TRIAL, VALID = 0, 1, 2

class PlotMode(Enum):
    NumericalSolution = 0
    TrueSolution = 1
    PointwiseError = 2

class Eikonal(Enum):
    DirectArrival = 0
    OFaceReflection = 1
    NFaceReflection = 2

class Field(Enum):
    T = 0
    Tx = 1
    Ty = 2
    Tz = 3
    Txx = 4
    Txy = 5
    Txz = 6
    Tyx = 7
    Tyy = 8
    Tyz = 9
    Tzx = 10
    Tzy = 11
    Tzz = 12
    Origin = 13

class MainWindow(pvqt.MainWindow):
    def __init__(self, parent=None, show=True):
        QtWidgets.QMainWindow.__init__(self, parent)

        self.pointIndex = [None, None, None, None]

        self.loadData()

        # start setting up main frame and layout

        self.frame = QtWidgets.QFrame()
        hlayout = QtWidgets.QHBoxLayout()

        # set up first panel

        frame1 = QtWidgets.QFrame()
        vlayout1 = QtWidgets.QVBoxLayout()

        self.reloadDataPushButton = QtWidgets.QPushButton('Reload data')
        self.makeParameterFrame()

        self.reloadDataPushButton.clicked.connect(self.loadData)

        vlayout1.addWidget(self.reloadDataPushButton)
        vlayout1.addWidget(self.parameterFrame)
        vlayout1.addStretch(1)

        frame1.setLayout(vlayout1)
        hlayout.addWidget(frame1, stretch=1)

        # set up second panel

        frame2 = QtWidgets.QFrame()
        vlayout2 = QtWidgets.QVBoxLayout()

        self.makePlotModeGroupBox()
        self.makeEikonalSelectGroupBox()
        self.makeFieldSelectGroupBox()

        vlayout2.addWidget(self.plotModeGroupBox)
        vlayout2.addWidget(self.eikonalSelectGroupBox)
        vlayout2.addWidget(self.fieldSelectGroupBox)
        vlayout2.addStretch(1)

        frame2.setLayout(vlayout2)
        hlayout.addWidget(frame2, stretch=1)

        # set up third panel

        frame3 = QtWidgets.QFrame()
        vlayout3 = QtWidgets.QVBoxLayout()

        self.makePointSelectionFrame()
        self.togglePointsFrame()

        vlayout3.addWidget(self.pointSelectionFrame)
        vlayout3.addWidget(self.togglePointsFrame)
        vlayout3.addStretch(1)

        frame3.setLayout(vlayout3)
        hlayout.addWidget(frame3, stretch=1)

        # set up plotter

        self.plotter = pvqt.QtInteractor(self.frame)

        self.signal_close.connect(self.plotter.close)

        hlayout.addWidget(self.plotter.interactor, stretch=10)

        # finish setting up main layout

        self.frame.setLayout(hlayout)
        self.setCentralWidget(self.frame)

        # start plotting

        self.first_plot = True
        self.updatePlot()
        if show:
            self.show()

    def makeParameterFrame(self):
        self.parameterFrame = QtWidgets.QFrame()

        nLineEdit = QtWidgets.QLineEdit()
        maxvolLineEdit = QtWidgets.QLineEdit()
        rfacLineEdit = QtWidgets.QLineEdit()
        phipLineEdit = QtWidgets.QLineEdit()
        spLineEdit = QtWidgets.QLineEdit()
        widthLineEdit = QtWidgets.QLineEdit()
        heightLineEdit = QtWidgets.QLineEdit()

        formLayout = QtWidgets.QFormLayout()

        formLayout.addRow("-n", nLineEdit)
        formLayout.addRow("--maxvol", maxvolLineEdit)
        formLayout.addRow("--rfac", rfacLineEdit)
        formLayout.addRow("--phip", phipLineEdit)
        formLayout.addRow("--sp", spLineEdit)
        formLayout.addRow("--width", widthLineEdit)
        formLayout.addRow("--height", heightLineEdit)

        self.parameterFrame.setLayout(formLayout)

    def makePlotModeGroupBox(self):
        self.plotModeGroupBox = QtWidgets.QGroupBox("Plot mode:")

        plotNumerical = QtWidgets.QRadioButton("Numerical solution")
        plotGroundtruth = QtWidgets.QRadioButton("True solution")
        plotError = QtWidgets.QRadioButton("Pointwise error")
        plotNumerical.setChecked(True)

        plotNumerical.toggled.connect(self.updatePlotOnSelect)
        plotGroundtruth.toggled.connect(self.updatePlotOnSelect)
        plotError.toggled.connect(self.updatePlotOnSelect)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(plotNumerical)
        layout.addWidget(plotGroundtruth)
        layout.addWidget(plotError)
        layout.addStretch(1)

        self.plotModeGroupBox.setLayout(layout)

    def makeEikonalSelectGroupBox(self):
        self.eikonalSelectGroupBox = QtWidgets.QGroupBox("Eikonal to plot:")

        plotDirect = QtWidgets.QRadioButton("Direct arrival")
        plotOReflection = QtWidgets.QRadioButton("Reflection (o-face)")
        plotNReflection = QtWidgets.QRadioButton("Reflection (n-face)")
        plotDirect.setChecked(True)

        plotDirect.toggled.connect(self.updatePlotOnSelect)
        plotOReflection.toggled.connect(self.updatePlotOnSelect)
        plotNReflection.toggled.connect(self.updatePlotOnSelect)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(plotDirect)
        layout.addWidget(plotOReflection)
        layout.addWidget(plotNReflection)
        layout.addStretch(1)

        self.eikonalSelectGroupBox.setLayout(layout)

    def makeFieldSelectGroupBox(self):
        self.fieldSelectGroupBox = QtWidgets.QGroupBox("Field to plot:")

        plotT = QtWidgets.QRadioButton("T")
        plotTx = QtWidgets.QRadioButton("∂T/∂x")
        plotTy = QtWidgets.QRadioButton("∂T/∂y")
        plotTz = QtWidgets.QRadioButton("∂T/∂z")
        plotTxx = QtWidgets.QRadioButton("∂²T/∂x²")
        plotTxy = QtWidgets.QRadioButton("∂²T/∂y∂x")
        plotTxz = QtWidgets.QRadioButton("∂²T/∂z∂x")
        plotTyx = QtWidgets.QRadioButton("∂²T/∂x∂y")
        plotTyy = QtWidgets.QRadioButton("∂²T/∂y²")
        plotTyz = QtWidgets.QRadioButton("∂²T/∂z∂y")
        plotTzx = QtWidgets.QRadioButton("∂²T/∂x∂z")
        plotTzy = QtWidgets.QRadioButton("∂²T/∂y∂z")
        plotTzz = QtWidgets.QRadioButton("∂²T/∂z²")
        plotOrigin = QtWidgets.QRadioButton("Origin")
        plotT.setChecked(True)

        plotT.toggled.connect(self.updatePlotOnSelect)
        plotTx.toggled.connect(self.updatePlotOnSelect)
        plotTy.toggled.connect(self.updatePlotOnSelect)
        plotTz.toggled.connect(self.updatePlotOnSelect)
        plotTxx.toggled.connect(self.updatePlotOnSelect)
        plotTxy.toggled.connect(self.updatePlotOnSelect)
        plotTxz.toggled.connect(self.updatePlotOnSelect)
        plotTyx.toggled.connect(self.updatePlotOnSelect)
        plotTyy.toggled.connect(self.updatePlotOnSelect)
        plotTyz.toggled.connect(self.updatePlotOnSelect)
        plotTzx.toggled.connect(self.updatePlotOnSelect)
        plotTzy.toggled.connect(self.updatePlotOnSelect)
        plotTzz.toggled.connect(self.updatePlotOnSelect)
        plotOrigin.toggled.connect(self.updatePlotOnSelect)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(plotT)
        layout.addWidget(plotTx)
        layout.addWidget(plotTy)
        layout.addWidget(plotTz)
        layout.addWidget(plotTxx)
        layout.addWidget(plotTxy)
        layout.addWidget(plotTxz)
        layout.addWidget(plotTyx)
        layout.addWidget(plotTyy)
        layout.addWidget(plotTyz)
        layout.addWidget(plotTzx)
        layout.addWidget(plotTzy)
        layout.addWidget(plotTzz)
        layout.addWidget(plotOrigin)
        layout.addStretch(1)

        self.fieldSelectGroupBox.setLayout(layout)

    def makePointSelectionFrame(self):
        self.pointSelectionFrame = QtWidgets.QGroupBox("Points to highlight:")

        index1LineEdit = QtWidgets.QLineEdit()
        index2LineEdit = QtWidgets.QLineEdit()
        index3LineEdit = QtWidgets.QLineEdit()
        index4LineEdit = QtWidgets.QLineEdit()

        def pointIndexUpdated(i, lineEdit):
            s = lineEdit.text()
            try:
                self.pointIndex[i] = int(s)
                if 0 <= self.pointIndex[i] < self.verts.shape[0]:
                    self.updatePlot()
                else:
                    print(f'index {self.pointIndex[i]} out of range [0, {self.verts.shape[0]})')
                    self.pointIndex[i] = None
            except:
                print(f'string "{s}" is not a valid index')
                self.pointIndex[i] = None

        index1LineEdit.editingFinished.connect(
            lambda: pointIndexUpdated(0, index1LineEdit))

        index2LineEdit.editingFinished.connect(
            lambda: pointIndexUpdated(1, index2LineEdit))

        index3LineEdit.editingFinished.connect(
            lambda: pointIndexUpdated(2, index3LineEdit))

        index4LineEdit.editingFinished.connect(
            lambda: pointIndexUpdated(3, index4LineEdit))

        formLayout = QtWidgets.QFormLayout()
        formLayout.addRow(f'Point #1 index:', index1LineEdit)
        formLayout.addRow(f'Point #2 index:', index2LineEdit)
        formLayout.addRow(f'Point #3 index:', index3LineEdit)
        formLayout.addRow(f'Point #4 index:', index4LineEdit)

        self.pointSelectionFrame.setLayout(formLayout)

    def togglePointsFrame(self):
        self.togglePointsFrame = QtWidgets.QFrame()

        showDirectPointsCheckBox = QtWidgets.QCheckBox()
        showDiffPointsCheckBox = QtWidgets.QCheckBox()

        showDirectPointsCheckBox.setCheckState(2)
        showDiffPointsCheckBox.setCheckState(2)

        formLayout = QtWidgets.QFormLayout()
        formLayout.addRow('Show direct zone:', showDirectPointsCheckBox)
        formLayout.addRow('Show diff. zone:', showDiffPointsCheckBox)

        showDirectPointsCheckBox.stateChanged.connect(self.updatePlot)
        showDiffPointsCheckBox.stateChanged.connect(self.updatePlot)

        self.showDirectZone = lambda: showDirectPointsCheckBox.checkState() == 2
        self.showDiffZone = lambda: showDiffPointsCheckBox.checkState() == 2

        self.togglePointsFrame.setLayout(formLayout)

    def reloadData(self):
        self.loadData()
        self.updatePlot()
        self.show()

    def loadData(self):
        # load vertices and cells of the tetrahedron mesh discretizing
        # the domain
        self.verts = np.fromfile('verts.bin', dtype=np.float64).reshape(-1, 3)
        self.cells = np.fromfile('cells.bin', dtype=np.uintp).reshape(-1, 4)

        # make a version of it in PyVista
        self.mesh = pv.UnstructuredGrid({vtk.VTK_TETRA: self.cells}, self.verts)

        # load the problem specification
        with open('spec.txt', 'r') as f:
            self.spec = {k: v.strip() for k, v in map(lambda s: s.split(':'), f)}
        for k, Type in {
                'verbose': bool,
                'visualize': bool,
                'maxvol': float,
                'n': float,
                'w': float,
                'h': float,
                'R': float}.items():
            self.spec[k] = Type(self.spec[k])

        # print it out
        #
        # TODO: validate it and make sure the values agree with the
        # requested parameters
        print('problem specification (3d wedge):')
        for k, v in self.spec.items():
            print(f'- {k}: {v}')

        ## LOAD DATA FOR DIRECT EIKONAL

        # numerical solution
        self.direct_jet = np.fromfile('direct_jet.bin', dtype=np.float64).reshape(-1, 4)
        self.direct_T = self.direct_jet[:, 0]
        self.direct_grad_T = self.direct_jet[:, 1:]
        self.direct_hess_T = np.fromfile('direct_hess.bin', dtype=np.float64).reshape(-1, 3, 3)

        # final states
        self.direct_state = np.fromfile('direct_state.bin', dtype=np.intc)

        self.direct_par_l = np.fromfile('direct_par_l.bin', dtype=np.uintp).reshape(-1, 3)

        # approximate origins
        self.direct_origin = np.fromfile('direct_origin.bin', dtype=np.float64)

        # true solution ("gt")
        self.direct_jet_gt = np.fromfile('direct_jet_gt.bin', dtype=np.float64).reshape(-1, 13)
        self.direct_T_gt = self.direct_jet_gt[:, 0]
        self.direct_grad_T_gt = self.direct_jet_gt[:, 1:4]
        self.direct_hess_T_gt = self.direct_jet_gt[:, 4:].reshape(-1, 3, 3)

        # errors
        self.direct_error_jet = self.direct_jet - self.direct_jet_gt[:, :4]
        self.direct_error_T = self.direct_error_jet[:, 0]
        self.direct_error_grad_T = self.direct_error_jet[:, 1:]
        self.direct_error_hess_T = self.direct_hess_T - self.direct_hess_T_gt

        ## LOAD DATA FOR o-FACE REFLECTION

        # numerical solution
        self.o_refl_jet = np.fromfile('o_refl_jet.bin', dtype=np.float64).reshape(-1, 4)
        self.o_refl_T = self.o_refl_jet[:, 0]
        self.o_refl_grad_T = self.o_refl_jet[:, 1:]
        self.o_refl_hess_T = np.fromfile('o_refl_hess.bin', dtype=np.float64).reshape(-1, 3, 3)

        # final states
        self.o_refl_state = np.fromfile('o_refl_state.bin', dtype=np.intc)

        self.o_refl_par_l = np.fromfile('o_refl_par_l.bin', dtype=np.uintp).reshape(-1, 3)

        # approximate origins
        self.o_refl_origin = np.fromfile('o_refl_origin.bin', dtype=np.float64)

        # true solution
        self.o_refl_jet_gt = np.fromfile('o_refl_jet_gt.bin', dtype=np.float64).reshape(-1, 13)
        self.o_refl_T_gt = self.o_refl_jet_gt[:, 0]
        self.o_refl_grad_T_gt = self.o_refl_jet_gt[:, 1:4]
        self.o_refl_hess_T_gt = self.o_refl_jet_gt[:, 4:].reshape(-1, 3, 3)

        # errors
        self.o_refl_error_jet = self.o_refl_jet - self.o_refl_jet_gt[:, :4]
        self.o_refl_error_T = self.o_refl_error_jet[:, 0]
        self.o_refl_error_grad_T = self.o_refl_error_jet[:, 1:]
        self.o_refl_error_hess_T = self.o_refl_hess_T - self.o_refl_hess_T_gt

        ## LOAD DATA FOR n-FACE REFLECTION

        # numerical solution
        self.n_refl_jet = np.fromfile('n_refl_jet.bin', dtype=np.float64).reshape(-1, 4)
        self.n_refl_T = self.n_refl_jet[:, 0]
        self.n_refl_grad_T = self.n_refl_jet[:, 1:]
        self.n_refl_hess_T = np.fromfile('n_refl_hess.bin', dtype=np.float64).reshape(-1, 3, 3)

        # final states
        self.n_refl_state = np.fromfile('n_refl_state.bin', dtype=np.intc)

        self.n_refl_par_l = np.fromfile('n_refl_par_l.bin', dtype=np.uintp).reshape(-1, 3)

        # approximate origins
        self.n_refl_origin = np.fromfile('n_refl_origin.bin', dtype=np.float64)

        # true solution
        self.n_refl_jet_gt = np.fromfile('n_refl_jet_gt.bin', dtype=np.float64).reshape(-1, 13)
        self.n_refl_T_gt = self.n_refl_jet_gt[:, 0]
        self.n_refl_grad_T_gt = self.n_refl_jet_gt[:, 1:4]
        self.n_refl_hess_T_gt = self.n_refl_jet_gt[:, 4:].reshape(-1, 3, 3)

        # errors
        self.n_refl_error_jet = self.n_refl_jet - self.n_refl_jet_gt[:, :4]
        self.n_refl_error_T = self.n_refl_error_jet[:, 0]
        self.n_refl_error_grad_T = self.n_refl_error_jet[:, 1:4]
        self.n_refl_error_hess_T = self.n_refl_hess_T - self.n_refl_hess_T_gt

    def getPlotMode(self):
        children = self.plotModeGroupBox.findChildren(QtWidgets.QRadioButton)
        index = next(i for i, _ in enumerate(children) if _.isChecked())
        return PlotMode(index)

    def getSelectedEikonal(self):
        children = self.eikonalSelectGroupBox.findChildren(QtWidgets.QRadioButton)
        index = next(i for i, _ in enumerate(children) if _.isChecked())
        return Eikonal(index)

    def getSelectedField(self):
        children = self.fieldSelectGroupBox.findChildren(QtWidgets.QRadioButton)
        index = next(i for i, _ in enumerate(children) if _.isChecked())
        return Field(index)

    def updatePlotOnSelect(self, checked):
        if checked:
            self.updatePlot()

    def getCurrentValues(self):
        plotMode = self.getPlotMode()
        eikonalSel = self.getSelectedEikonal()
        fieldSel = self.getSelectedField()
        if plotMode is PlotMode.NumericalSolution:
            if eikonalSel is Eikonal.DirectArrival:
                if fieldSel is Field.T:        return self.direct_T
                elif fieldSel is Field.Tx:     return self.direct_grad_T[:, 0]
                elif fieldSel is Field.Ty:     return self.direct_grad_T[:, 1]
                elif fieldSel is Field.Tz:     return self.direct_grad_T[:, 2]
                elif fieldSel is Field.Txx:    return self.direct_hess_T[:, 0, 0]
                elif fieldSel is Field.Txy:    return self.direct_hess_T[:, 0, 1]
                elif fieldSel is Field.Txz:    return self.direct_hess_T[:, 0, 2]
                elif fieldSel is Field.Tyx:    return self.direct_hess_T[:, 1, 0]
                elif fieldSel is Field.Tyy:    return self.direct_hess_T[:, 1, 1]
                elif fieldSel is Field.Tyz:    return self.direct_hess_T[:, 1, 2]
                elif fieldSel is Field.Tzx:    return self.direct_hess_T[:, 2, 0]
                elif fieldSel is Field.Tzy:    return self.direct_hess_T[:, 2, 1]
                elif fieldSel is Field.Tzz:    return self.direct_hess_T[:, 2, 2]
                elif fieldSel is Field.Origin: return self.direct_origin
            elif eikonalSel is Eikonal.OFaceReflection:
                if fieldSel is Field.T:        return self.o_refl_T
                elif fieldSel is Field.Tx:     return self.o_refl_grad_T[:, 0]
                elif fieldSel is Field.Ty:     return self.o_refl_grad_T[:, 1]
                elif fieldSel is Field.Tz:     return self.o_refl_grad_T[:, 2]
                elif fieldSel is Field.Txx:    return self.o_refl_hess_T[:, 0, 0]
                elif fieldSel is Field.Txy:    return self.o_refl_hess_T[:, 0, 1]
                elif fieldSel is Field.Txz:    return self.o_refl_hess_T[:, 0, 2]
                elif fieldSel is Field.Tyx:    return self.o_refl_hess_T[:, 1, 0]
                elif fieldSel is Field.Tyy:    return self.o_refl_hess_T[:, 1, 1]
                elif fieldSel is Field.Tyz:    return self.o_refl_hess_T[:, 1, 2]
                elif fieldSel is Field.Tzx:    return self.o_refl_hess_T[:, 2, 0]
                elif fieldSel is Field.Tzy:    return self.o_refl_hess_T[:, 2, 1]
                elif fieldSel is Field.Tzz:    return self.o_refl_hess_T[:, 2, 2]
                elif fieldSel is Field.Origin: return self.o_refl_origin
            elif eikonalSel is Eikonal.NFaceReflection:
                if fieldSel is Field.T:        return self.n_refl_T
                elif fieldSel is Field.Tx:     return self.n_refl_grad_T[:, 0]
                elif fieldSel is Field.Ty:     return self.n_refl_grad_T[:, 1]
                elif fieldSel is Field.Tz:     return self.n_refl_grad_T[:, 2]
                elif fieldSel is Field.Txx:    return self.n_refl_hess_T[:, 0, 0]
                elif fieldSel is Field.Txy:    return self.n_refl_hess_T[:, 0, 1]
                elif fieldSel is Field.Txz:    return self.n_refl_hess_T[:, 0, 2]
                elif fieldSel is Field.Tyx:    return self.n_refl_hess_T[:, 1, 0]
                elif fieldSel is Field.Tyy:    return self.n_refl_hess_T[:, 1, 1]
                elif fieldSel is Field.Tyz:    return self.n_refl_hess_T[:, 1, 2]
                elif fieldSel is Field.Tzx:    return self.n_refl_hess_T[:, 2, 0]
                elif fieldSel is Field.Tzy:    return self.n_refl_hess_T[:, 2, 1]
                elif fieldSel is Field.Tzz:    return self.n_refl_hess_T[:, 2, 2]
                elif fieldSel is Field.Origin: return self.n_refl_origin
        elif plotMode is PlotMode.TrueSolution:
            if eikonalSel is Eikonal.DirectArrival:
                if fieldSel is Field.T:        return self.direct_T_gt
                elif fieldSel is Field.Tx:     return self.direct_grad_T_gt[:, 0]
                elif fieldSel is Field.Ty:     return self.direct_grad_T_gt[:, 1]
                elif fieldSel is Field.Tz:     return self.direct_grad_T_gt[:, 2]
                elif fieldSel is Field.Txx:    return self.direct_hess_T_gt[:, 0, 0]
                elif fieldSel is Field.Txy:    return self.direct_hess_T_gt[:, 0, 1]
                elif fieldSel is Field.Txz:    return self.direct_hess_T_gt[:, 0, 2]
                elif fieldSel is Field.Tyx:    return self.direct_hess_T_gt[:, 1, 0]
                elif fieldSel is Field.Tyy:    return self.direct_hess_T_gt[:, 1, 1]
                elif fieldSel is Field.Tyz:    return self.direct_hess_T_gt[:, 1, 2]
                elif fieldSel is Field.Tzx:    return self.direct_hess_T_gt[:, 2, 0]
                elif fieldSel is Field.Tzy:    return self.direct_hess_T_gt[:, 2, 1]
                elif fieldSel is Field.Tzz:    return self.direct_hess_T_gt[:, 2, 2]
                elif fieldSel is Field.Origin: return self.direct_origin
            elif eikonalSel is Eikonal.OFaceReflection:
                if fieldSel is Field.T:        return self.o_refl_T_gt
                elif fieldSel is Field.Tx:     return self.o_refl_grad_T_gt[:, 0]
                elif fieldSel is Field.Ty:     return self.o_refl_grad_T_gt[:, 1]
                elif fieldSel is Field.Tz:     return self.o_refl_grad_T_gt[:, 2]
                elif fieldSel is Field.Txx:    return self.o_refl_hess_T_gt[:, 0, 0]
                elif fieldSel is Field.Txy:    return self.o_refl_hess_T_gt[:, 0, 1]
                elif fieldSel is Field.Txz:    return self.o_refl_hess_T_gt[:, 0, 2]
                elif fieldSel is Field.Tyx:    return self.o_refl_hess_T_gt[:, 1, 0]
                elif fieldSel is Field.Tyy:    return self.o_refl_hess_T_gt[:, 1, 1]
                elif fieldSel is Field.Tyz:    return self.o_refl_hess_T_gt[:, 1, 2]
                elif fieldSel is Field.Tzx:    return self.o_refl_hess_T_gt[:, 2, 0]
                elif fieldSel is Field.Tzy:    return self.o_refl_hess_T_gt[:, 2, 1]
                elif fieldSel is Field.Tzz:    return self.o_refl_hess_T_gt[:, 2, 2]
                elif fieldSel is Field.Origin: return self.o_refl_origin
            elif eikonalSel is Eikonal.NFaceReflection:
                if fieldSel is Field.T:        return self.n_refl_T_gt
                elif fieldSel is Field.Tx:     return self.n_refl_grad_T_gt[:, 0]
                elif fieldSel is Field.Ty:     return self.n_refl_grad_T_gt[:, 1]
                elif fieldSel is Field.Tz:     return self.n_refl_grad_T_gt[:, 2]
                elif fieldSel is Field.Txx:    return self.n_refl_hess_T_gt[:, 0, 0]
                elif fieldSel is Field.Txy:    return self.n_refl_hess_T_gt[:, 0, 1]
                elif fieldSel is Field.Txz:    return self.n_refl_hess_T_gt[:, 0, 2]
                elif fieldSel is Field.Tyx:    return self.n_refl_hess_T_gt[:, 1, 0]
                elif fieldSel is Field.Tyy:    return self.n_refl_hess_T_gt[:, 1, 1]
                elif fieldSel is Field.Tyz:    return self.n_refl_hess_T_gt[:, 1, 2]
                elif fieldSel is Field.Tzx:    return self.n_refl_hess_T_gt[:, 2, 0]
                elif fieldSel is Field.Tzy:    return self.n_refl_hess_T_gt[:, 2, 1]
                elif fieldSel is Field.Tzz:    return self.n_refl_hess_T_gt[:, 2, 2]
                elif fieldSel is Field.Origin: return self.n_refl_origin
        elif plotMode is PlotMode.PointwiseError:
            if eikonalSel is Eikonal.DirectArrival:
                if fieldSel is Field.T:        return self.direct_error_T
                elif fieldSel is Field.Tx:     return self.direct_error_grad_T[:, 0]
                elif fieldSel is Field.Ty:     return self.direct_error_grad_T[:, 1]
                elif fieldSel is Field.Tz:     return self.direct_error_grad_T[:, 2]
                elif fieldSel is Field.Txx:    return self.direct_error_hess_T[:, 0, 0]
                elif fieldSel is Field.Txy:    return self.direct_error_hess_T[:, 0, 1]
                elif fieldSel is Field.Txz:    return self.direct_error_hess_T[:, 0, 2]
                elif fieldSel is Field.Tyx:    return self.direct_error_hess_T[:, 1, 0]
                elif fieldSel is Field.Tyy:    return self.direct_error_hess_T[:, 1, 1]
                elif fieldSel is Field.Tyz:    return self.direct_error_hess_T[:, 1, 2]
                elif fieldSel is Field.Tzx:    return self.direct_error_hess_T[:, 2, 0]
                elif fieldSel is Field.Tzy:    return self.direct_error_hess_T[:, 2, 1]
                elif fieldSel is Field.Tzz:    return self.direct_error_hess_T[:, 2, 2]
                elif fieldSel is Field.Origin: return self.direct_origin
            elif eikonalSel is Eikonal.OFaceReflection:
                if fieldSel is Field.T:        return self.o_refl_error_T
                elif fieldSel is Field.Tx:     return self.o_refl_error_grad_T[:, 0]
                elif fieldSel is Field.Ty:     return self.o_refl_error_grad_T[:, 1]
                elif fieldSel is Field.Tz:     return self.o_refl_error_grad_T[:, 2]
                elif fieldSel is Field.Txx:    return self.o_refl_error_hess_T[:, 0, 0]
                elif fieldSel is Field.Txy:    return self.o_refl_error_hess_T[:, 0, 1]
                elif fieldSel is Field.Txz:    return self.o_refl_error_hess_T[:, 0, 2]
                elif fieldSel is Field.Tyx:    return self.o_refl_error_hess_T[:, 1, 0]
                elif fieldSel is Field.Tyy:    return self.o_refl_error_hess_T[:, 1, 1]
                elif fieldSel is Field.Tyz:    return self.o_refl_error_hess_T[:, 1, 2]
                elif fieldSel is Field.Tzx:    return self.o_refl_error_hess_T[:, 2, 0]
                elif fieldSel is Field.Tzy:    return self.o_refl_error_hess_T[:, 2, 1]
                elif fieldSel is Field.Tzz:    return self.o_refl_error_hess_T[:, 2, 2]
                elif fieldSel is Field.Origin: return self.o_refl_origin
            elif eikonalSel is Eikonal.NFaceReflection:
                if fieldSel is Field.T:        return self.n_refl_error_T
                elif fieldSel is Field.Tx:     return self.n_refl_error_grad_T[:, 0]
                elif fieldSel is Field.Ty:     return self.n_refl_error_grad_T[:, 1]
                elif fieldSel is Field.Tz:     return self.n_refl_error_grad_T[:, 2]
                elif fieldSel is Field.Txx:    return self.n_refl_error_hess_T[:, 0, 0]
                elif fieldSel is Field.Txy:    return self.n_refl_error_hess_T[:, 0, 1]
                elif fieldSel is Field.Txz:    return self.n_refl_error_hess_T[:, 0, 2]
                elif fieldSel is Field.Tyx:    return self.n_refl_error_hess_T[:, 1, 0]
                elif fieldSel is Field.Tyy:    return self.n_refl_error_hess_T[:, 1, 1]
                elif fieldSel is Field.Tyz:    return self.n_refl_error_hess_T[:, 1, 2]
                elif fieldSel is Field.Tzx:    return self.n_refl_error_hess_T[:, 2, 0]
                elif fieldSel is Field.Tzy:    return self.n_refl_error_hess_T[:, 2, 1]
                elif fieldSel is Field.Tzz:    return self.n_refl_error_hess_T[:, 2, 2]
                elif fieldSel is Field.Origin: return self.n_refl_origin
        raise RuntimeError(f"couldn't get values for combination: {plotMode}, {eikonalSel}, {fieldSel}")

    def getMask(self):
        showDirect = self.showDirectZone()
        showDiff = self.showDiffZone()
        eikonalSel = self.getSelectedEikonal()
        if eikonalSel is Eikonal.DirectArrival:
            origin = self.direct_origin
        elif eikonalSel is Eikonal.OFaceReflection:
            origin = self.o_refl_origin
        elif eikonalSel is Eikonal.NFaceReflection:
            origin = self.n_refl_origin
        else:
            assert False
        if showDirect and showDiff:
            mask = np.ones(origin.shape, dtype=np.bool_)
        elif showDirect:
            mask = origin >= 0.5
        elif showDiff:
            mask = origin <= 0.5
        else:
            mask = np.zeros(origin.shape, dtype=np.bool_)
        assert mask.size == self.verts.shape[0]
        return mask

    def getCMap(self):
        if self.getSelectedField() is Field.Origin:
            return cc.cm.fire
        elif self.getSelectedField() is Field.T and \
           self.getPlotMode() is not PlotMode.PointwiseError:
            return cc.cm.rainbow
        else:
            return cc.cm.coolwarm

    def getCLim(self, values):
        plotMode = self.getPlotMode()
        if plotMode == PlotMode.PointwiseError:
            print((values == 0).sum())
            vmax = np.nanmax(abs(values))
            return (-vmax, vmax)
        else:
            return (np.nanmin(values), np.nanmax(values))

    def updatePlot(self):
        values = self.getCurrentValues()
        values[~np.isfinite(values)] = np.nan

        mask = self.getMask()

        if mask.size > 0:
            self.poly_data = pv.PolyData(self.verts[mask])
            self.poly_data['values'] = values[mask]

        if not self.first_plot:
            old_camera = self.plotter.camera.copy()

        self.plotter.clear()

        self.plotter.add_mesh(
            self.mesh,
            show_edges=True,
            opacity=0.25,
        )

        if mask.size > 0:
            self.plotter.add_mesh(
                self.poly_data,
                scalars='values',
                cmap=self.getCMap(),
                clim=self.getCLim(values),
                # nan_color='magenta',
                nan_opacity=0.0,
            )

        for i in range(4):
            l = self.pointIndex[i]
            if l is not None:
                x = self.verts[l]
                r = 0.025
                c = cc.cm.rainbow(float(i)/float(3))
                self.plotter.add_mesh(pv.Sphere(r, x), color=c)

        if self.first_plot:
            self.first_plot = False
        else:
            self.plotter.camera = old_camera

        def set_picked_point_index(grid):
            x = self.plotter.picked_point
            dists = np.sqrt(np.sum((x - self.verts)**2, axis=1))
            i = np.argmin(dists)
            self.picked_point_index = i
            print(f'selected vertex {self.picked_point_index}')
            print(f'- parents:', end='')
            if self.getSelectedEikonal() == Eikonal.DirectArrival:
                print(self.direct_par_l[i])
            elif self.getSelected() == Eikonal.OFaceReflection:
                print(self.o_refl_par_l[i])
            else:
                print(self.n_refl_par_l[i])

        self.plotter.enable_point_picking(callback=set_picked_point_index, font_size=12)

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())
