#!/usr/bin/env python3
#coding utf-8

import sys 
import math

import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

import lat
import lat.read_atoms

SCALE = 20
RADIUS = 5
WIDTH = 5
BIG_INT = 1000000

class Atom():
    def __init__(self, atomNumber,
                       atomGroup,
                       atomType,
                       atomCharge,
                       x, y, z):
        self.values = {}
        self.values['atomNumber'] = atomNumber
        self.values['atomGroup'] = atomGroup
        self.values['atomType'] = atomType
        self.values['atomCharge'] = atomCharge
        self.values['X'] = x
        self.values['Y'] = y
        self.values['Z'] = z

class Bounds():
    def __init__(self, xlo, xhi,
                       ylo, yhi,
                       zlo, zhi):
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi
        self.lx = xhi - xlo
        self.ly = yhi - ylo
        self.lz = zhi - zlo

class Bond():
    def __init__(self, bondNumber,
                       bondType,
                       bondAtomOneNumber,
                       bondAtomTwoNumber):
        self.values = {}
        self.values['bondNumber'] = bondNumber
        self.values['bondType'] = bondType
        self.values['bondAtomOneNumber'] = bondAtomOneNumber
        self.values['bondAtomTwoNumber'] = bondAtomTwoNumber
        self.values['bondBegin'] = {}
        self.values['bondEnd'] = {}

    def calculateEnds(self, atomsList, molecule):
        l = len(molecule)
        atomsInMolecule = [molecule[i].values['atomNumber'] for i in range(l)]
        for atom in molecule:
            if atom.values['atomNumber'] == self.values['bondAtomOneNumber']:
                self.values['bondBegin']['X'] = atom.values['X']
                self.values['bondBegin']['Y'] = atom.values['Y']
                self.values['bondBegin']['Z'] = atom.values['Z']
            if atom.values['atomNumber'] == self.values['bondAtomTwoNumber']:
                self.values['bondEnd']['X'] = atom.values['X']
                self.values['bondEnd']['Y'] = atom.values['Y']
                self.values['bondEnd']['Z'] = atom.values['Z']

class Molecule():
    def __init__(self, moleculeNumber, systemName):
        if systemName == '10x20':
            self.hydrogenTypes = [12, 13] # Интересен только скелет
            chainsPerCell = 10
            self.polymerLength = 382
            delta = 1560
            systemSize = 5380
            cellNumber = moleculeNumber // chainsPerCell
            inCellNumber = moleculeNumber % chainsPerCell
            startAtomNumber = (cellNumber * systemSize + 
                               delta + 
                               inCellNumber * self.polymerLength) + 1
            self.startAtomNumber = startAtomNumber

    def makeSkeleton(self, atomsList, bounds):
        """Формирует скелет цепочки с учётом пересечения границ"""
        molecule = []

        lx = bounds.lx
        ly = bounds.ly
        lz = bounds.lz

        startX = atomsList[self.startAtomNumber - 1].values['X']
        startY = atomsList[self.startAtomNumber - 1].values['Y']
        startZ = atomsList[self.startAtomNumber - 1].values['Z']
        molecule.append(Atom(self.startAtomNumber,
                             None,
                             None,
                             None,
                             startX,
                             startY,
                             startZ))

        startAtomNumber = self.startAtomNumber
        for i in range(self.polymerLength - 1):
            nearestAtomNumber = startAtomNumber + 1
            nearestAtomX = atomsList[nearestAtomNumber - 1].values['X']
            nearestAtomY = atomsList[nearestAtomNumber - 1].values['Y']
            nearestAtomZ = atomsList[nearestAtomNumber - 1].values['Z']
            r2old = BIG_INT
            coords = [0, 0, 0]
            for x in [-1, 0, 1]:
                nx = nearestAtomX + x * lx
                dx = abs(nx - startX)
                for y in [-1, 0, 1]:
                    ny = nearestAtomY + y * ly
                    dy = abs(ny - startY)
                    for z in [-1, 0, 1]:
                        nz = nearestAtomZ + z * lz
                        dz = abs(nz - startZ)
                        r2new = dx**2 + dy**2 + dz**2
                        if r2new < r2old:
                            r2old = r2new 
                            coords = [x, y, z]
            if not (atomsList[nearestAtomNumber - 1].values['atomType'] in
                    self.hydrogenTypes):
                molecule.append(Atom(nearestAtomNumber,
                                     None,
                                     None,
                                     None,
                                     nearestAtomX + coords[0] * lx,
                                     nearestAtomY + coords[1] * ly,
                                     nearestAtomZ + coords[2] * lz))
            startAtomNumber = nearestAtomNumber
            startX = nearestAtomX + coords[0] * lx                    
            startY = nearestAtomY + coords[1] * ly        
            startZ = nearestAtomZ + coords[2] * lz
        self.molecule = molecule

    def chooseBonds(self, atomsList, bondsList):
        """Отбирает те связи, что входят в данную молекулу
           Также проставляет каждой связи координаты концов"""
        atomsInMolecule = [atom.values['atomNumber'] for atom in self.molecule]
        bondsInMolecule = []
        for bond in bondsList:
            if (bond.values['bondAtomOneNumber'] in atomsInMolecule and
                bond.values['bondAtomTwoNumber'] in atomsInMolecule):
                bondsInMolecule.append(bond)
                bond.calculateEnds(None, self.molecule)
        self.bondsInMolecule = bondsInMolecule

    def makeBonds(self, atomsList, bondsList):
        for bond in self.bondsInMolecule:
            bond.calculateEnds(atomsList)

    def computeRanges(self):
        maxX = -BIG_INT
        minX = BIG_INT
        maxY = -BIG_INT
        minY = BIG_INT
        maxZ = -BIG_INT
        minZ = BIG_INT
        for atom in self.molecule:
            if atom.values['X'] < minX:
                minX = atom.values['X']
            if atom.values['X'] > maxX:
                maxX = atom.values['X']
            if atom.values['Y'] < minY:
                minY = atom.values['Y']
            if atom.values['Y'] > maxY:
                maxY = atom.values['Y']
            if atom.values['Z'] < minZ:
                minZ = atom.values['Z']
            if atom.values['Z'] > maxZ:
                maxZ = atom.values['Z']
        self.ranges = { 'minX': minX,
                        'maxX': maxX,
                        'minY': minY,
                        'maxY': maxY,
                        'minZ': minZ,
                        'maxZ': maxZ }

class MainWidget(QWidget):
    """Рисует молекулу (атомы и связи)"""
    def __init__(self, molecule,
                       coord1='X', coord2='Y'):
        super().__init__()
        self.molecule = molecule
        self.radius = RADIUS
        self.scale = SCALE
        self.width = WIDTH
        self.coord1 = coord1
        self.coord2 = coord2
        x = (self.scale * (self.molecule.ranges['max' + self.coord1] - 
                           self.molecule.ranges['min' + self.coord1]) +
             2 * self.radius + 
             self.width)
        y = (self.scale * (self.molecule.ranges['max' + self.coord2] -
                           self.molecule.ranges['min' + self.coord2]) +
             2 * self.radius + 
             self.width)
        self.resize(x, y)
        self.takeScreenshot()

    def paintEvent(self, e):
        qp = QPainter()
        qp.begin(self)
        brush = QBrush(QColor(0, 0, 0), Qt.SolidPattern)
        pen = QPen(brush, 5)
        qp.setBrush(brush)
        qp.setPen(pen)
        indices = { '1': 'X',
                    '2': 'Y',
                    '3': 'Z' }
        for atom in self.molecule.molecule:
            qp.drawEllipse(QPoint(self.scale * (atom.values[self.coord1] -
                                                self.molecule.ranges['min' +
                                                                     self.coord1]) +
                                  self.radius + 
                                  self.width / 2,
                                  self.scale * (atom.values[self.coord2] -
                                                self.molecule.ranges['min' +
                                                                     self.coord2]) +
                                  self.radius + 
                                  self.width / 2),
                           self.radius, self.radius)
        for bond in self.molecule.bondsInMolecule:
            qp.drawLine(self.scale * (bond.values['bondBegin'][self.coord1] -
                                      self.molecule.ranges['min' + self.coord1]) +
                        self.radius + 
                                  self.width / 2,
                        self.scale * (bond.values['bondBegin'][self.coord2] -
                                      self.molecule.ranges['min' + self.coord2]) +
                        self.radius + 
                                  self.width / 2,
                        self.scale * (bond.values['bondEnd'][self.coord1] -
                                      self.molecule.ranges['min' + self.coord1]) +
                        self.radius + 
                                  self.width / 2,
                        self.scale * (bond.values['bondEnd'][self.coord2] -
                                      self.molecule.ranges['min' + self.coord2]) +
                        self.radius + 
                                  self.width / 2)

        qp.end()

    def takeScreenshot(self, fname='XY.png', format='png'):
        p = self.grab();
        self.moleculeNumber = 1
        fname2 = str(self.moleculeNumber) + self.coord1 + self.coord2 + '.png'
        p.save(fname2, format, -1)

def main():
    app = QApplication(sys.argv)

    fname = '/home/anton/DataExamples/10x20.data'
    if fname == '/home/anton/DataExamples/10x20.data':
        systemName = '10x20'
        moleculeNumber = 90

    i = 0
    atomsList = []
    bondsList = []
    print('Parsing ', fname);
    [atoms, bounds, bonds, angles] = lat.read_atoms.read_atoms(fname)
    for atom in atoms:
        i += 1
        print('Preparing atoms: ', i, ' / ', len(bonds))
        at = Atom(atom[0], atom[1], atom[2],
                  atom[3], atom[4], atom[5], atom[6])
        atomsList.append(at)
    bounds = Bounds(bounds[0], bounds[1], bounds[2],
                    bounds[3], bounds[4], bounds[5])
    i = 0
    for bond in bonds:
        i += 1
        print('Preparing bonds: ', i, ' / ', len(bonds))
        bo = Bond(bond[0], bond[1], bond[2], bond[3])
        #bo.calculateEnds(atomsList)
        bondsList.append(bo)
    i = 0
    for i in range(1):
        print('Preparing molecules: ', i + 1, ' / ', moleculeNumber)
        mol = Molecule(i, systemName)
        mol.makeSkeleton(atomsList, bounds)
        mol.chooseBonds(atomsList, bondsList)
        mol.computeRanges()

        w = MainWidget(mol, 'X', 'Y')
        w = MainWidget(mol, 'X', 'Z')
        w = MainWidget(mol, 'Y', 'Z')

    app.quit()

main()
