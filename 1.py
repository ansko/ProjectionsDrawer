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
        self.atomNumber = atomNumber
        self.atomGroup = atomGroup
        self.atomType = atomType
        self.atomCharge = atomCharge
        self.x = x
        self.y = y
        self.z = z

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
        self.bondNumber = bondNumber
        self.bondType = bondType
        self.bondAtomOneNumber = bondAtomOneNumber
        self.bondAtomTwoNumber = bondAtomTwoNumber

    def calculateEnds(self, atomsList, molecule):
        atomsInMolecule = [molecule[i][0] for i in range(len(molecule))]
        for atom in molecule:
            if atom[0] == self.bondAtomOneNumber:
                self.bondBegin = [atom[1], atom[2], atom[3]]
            if atom[0] == self.bondAtomTwoNumber:
                self.bondEnd = [atom[1], atom[2], atom[3]]

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

        startX = atomsList[self.startAtomNumber - 1].x
        startY = atomsList[self.startAtomNumber - 1].y
        startZ = atomsList[self.startAtomNumber - 1].z
        molecule.append([self.startAtomNumber, startX, startY, startZ])

        startAtomNumber = self.startAtomNumber
        for i in range(self.polymerLength - 1):
            nearestAtomNumber = startAtomNumber + 1
            nearestAtomX = atomsList[nearestAtomNumber - 1].x
            nearestAtomY = atomsList[nearestAtomNumber - 1].y
            nearestAtomZ = atomsList[nearestAtomNumber - 1].z
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
            if not atomsList[nearestAtomNumber - 1].atomType in self.hydrogenTypes:
                molecule.append([nearestAtomNumber,
                                 nearestAtomX + coords[0] * lx,
                                 nearestAtomY + coords[1] * ly,
                                 nearestAtomZ + coords[2] * lz])
            startAtomNumber = nearestAtomNumber
            startX = nearestAtomX + coords[0] * lx                    
            startY = nearestAtomY + coords[1] * ly        
            startZ = nearestAtomZ + coords[2] * lz
        self.molecule = molecule

    def chooseBonds(self, atomsList, bondsList):
        """Отбирает те связи, что входят в данную молекулу
           Также проставляет каждой связи координаты концов"""
        atomsInMolecule = [atom[0] for atom in self.molecule]
        bondsInMolecule = []
        for bond in bondsList:
            if (bond.bondAtomOneNumber in atomsInMolecule and
                bond.bondAtomTwoNumber in atomsInMolecule):
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
            if atom[1] < minX:
                minX = atom[1]
            if atom[1] > maxX:
                maxX = atom[1]
            if atom[2] < minY:
                minY = atom[2]
            if atom[2] > maxY:
                maxY = atom[2]
            if atom[3] < minZ:
                minZ = atom[3]
            if atom[3] > maxZ:
                maxZ = atom[3]
        self.ranges = (minX, maxX, minY, maxY, minZ, maxZ)

class MainWidget(QWidget):
    """Рисует молекулу (атомы и связи)"""
    def __init__(self, molecule,
                       coord1=1, coord2=2):
        super().__init__()
        self.molecule = molecule
        self.radius = RADIUS
        self.scale = SCALE
        self.width = WIDTH
        self.coord1 = coord1
        self.coord2 = coord2
        x = (self.scale * (self.molecule.ranges[1] - 
                           self.molecule.ranges[0]) +
             2 * self.radius + 
             self.width)
        y = (self.scale * (self.molecule.ranges[3] -
                           self.molecule.ranges[2]) +
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
        for atom in self.molecule.molecule:
            qp.drawEllipse(QPoint(self.scale * (atom[1] -
                                                self.molecule.ranges[0]) +
                                  self.radius + 
                                  self.width / 2,
                                  self.scale * (atom[2] -
                                                self.molecule.ranges[2]) +
                                  self.radius + 
                                  self.width / 2),
                           self.radius, self.radius)
        for bond in self.molecule.bondsInMolecule:
            qp.drawLine(self.scale * (bond.bondBegin[0] -
                                      self.molecule.ranges[0]) +
                        self.radius + 
                                  self.width / 2,
                        self.scale * (bond.bondBegin[1] -
                                      self.molecule.ranges[2]) +
                        self.radius + 
                                  self.width / 2,
                        self.scale * (bond.bondEnd[0] -
                                      self.molecule.ranges[0]) +
                        self.radius + 
                                  self.width / 2,
                        self.scale * (bond.bondEnd[1] -
                                      self.molecule.ranges[2]) +
                        self.radius + 
                                  self.width / 2)

        qp.end()

    def takeScreenshot(self, fname='./1.png', format='png'):
        p = self.grab();
        p.save(fname, format, -1)

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

        w = MainWidget(mol, 1, 2)
        w = MainWidget(mol, 1, 3)
        w = MainWidget(mol, 2, 3)

    app.quit()

main()
