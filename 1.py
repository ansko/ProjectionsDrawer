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
SYSTEM_NAME = '10x20'
PHASE = 'polymer'

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
    def __init__(self, moleculeNumber):
        self.moleculeNumber = moleculeNumber
        modifiersPerCell = 12
        modifierLength = 70
        if SYSTEM_NAME == '10x20':
            self.hydrogenTypes = [8, 12, 13] # Интересен только скелет
            chainsPerCell = 10
            self.polymerLength = 382
            delta = 1560
            systemSize = 5380

            # for polymer
            if PHASE == 'polymer':
                cellNumber = moleculeNumber // chainsPerCell
                inCellNumber = moleculeNumber % chainsPerCell
                startAtomNumber = (cellNumber * systemSize + 
                                   delta + 
                                   inCellNumber * self.polymerLength) + 1
            # for modifier
            elif PHASE == 'modifier':
                cellNumber = moleculeNumber // modifiersPerCell
                inCellNumber = moleculeNumber % modifiersPerCell
                delta = 720
                startAtomNumber = (cellNumber * systemSize + 
                                   delta + 
                                   inCellNumber * modifierLength) + 1
        self.startAtomNumber = startAtomNumber

    def makeSkeleton(self, atomsList, bounds):
        """Формирует скелет цепочки с учётом пересечения границ"""
        molecule = []

        lx = bounds.lx
        ly = bounds.ly
        lz = bounds.lz

        num = self.startAtomNumber - 1
        startX = atomsList[num].values['X']
        startY = atomsList[num].values['Y']
        startZ = atomsList[num].values['Z']
        if not atomsList[num].values['atomType'] in self.hydrogenTypes:
            molecule.append(Atom(self.startAtomNumber,
                                 atomsList[num].values['atomGroup'],
                                 atomsList[num].values['atomType'],
                                 atomsList[num].values['atomCharge'],
                                 startX,
                                 startY,
                                 startZ))

        startAtomNumber = self.startAtomNumber
        if PHASE == 'polymer': # for polymer
            l = self.polymerLength - 1
        elif PHASE == 'modifier': # for modifier
            l = 69
        for i in range(l):
            num = startAtomNumber
            nearestAtomX = atomsList[num].values['X']
            nearestAtomY = atomsList[num].values['Y']
            nearestAtomZ = atomsList[num].values['Z']
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
            if not atomsList[num].values['atomType'] in self.hydrogenTypes:
                molecule.append(Atom(num + 1,
                                     atomsList[num].values['atomGroup'],
                                     atomsList[num].values['atomType'],
                                     atomsList[num].values['atomCharge'],
                                     nearestAtomX + coords[0] * lx,
                                     nearestAtomY + coords[1] * ly,
                                     nearestAtomZ + coords[2] * lz))
            startAtomNumber = num + 1
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
        indices = { '1': 'X',
                    '2': 'Y',
                    '3': 'Z' }
        brush = QBrush(QColor(0, 0, 0), Qt.SolidPattern)
        pen = QPen(brush, 5)
        qp.setPen(pen)
        dxyz = self.molecule.ranges
        for bond in self.molecule.bondsInMolecule:
            qp.drawLine(self.scale * (bond.values['bondBegin'][self.coord1] -
                                      dxyz['min' + self.coord1]) +
                        self.radius + 
                                  self.width / 2,
                        self.scale * (bond.values['bondBegin'][self.coord2] -
                                      dxyz['min' + self.coord2]) +
                        self.radius + 
                                  self.width / 2,
                        self.scale * (bond.values['bondEnd'][self.coord1] -
                                      dxyz['min' + self.coord1]) +
                        self.radius + 
                                  self.width / 2,
                        self.scale * (bond.values['bondEnd'][self.coord2] -
                                      dxyz['min' + self.coord2]) +
                        self.radius + 
                                  self.width / 2)
        for atom in self.molecule.molecule:
            if SYSTEM_NAME == '10x20':
                if atom.values['atomType'] in [4, 5, 6, 7, 15]:
                    brush = QBrush(QColor(255, 0, 0), Qt.SolidPattern)
                elif atom.values['atomType'] in [9, 16, 17]:
                    brush = QBrush(QColor(0, 255, 0), Qt.SolidPattern)
                elif atom.values['atomType'] in [10, 11, 14]:
                    brush = QBrush(QColor(105, 105, 105), Qt.SolidPattern)
                else: # неописанные атомы заметного розового цвета
                    print(atom.values['atomType'])
                    brush = QBrush(QColor(255, 0, 255), Qt.SolidPattern)
            qp.setBrush(brush)
            pen = QPen(brush, 5)
            qp.setPen(pen)
            qp.drawEllipse(QPoint(self.scale * (atom.values[self.coord1] -
                                                dxyz['min' + self.coord1]) +
                                  self.radius + 
                                  self.width / 2,
                                  self.scale * (atom.values[self.coord2] -
                                                dxyz['min' + self.coord2]) +
                                  self.radius + 
                                  self.width / 2),
                           self.radius, self.radius)
        qp.end()

    def takeScreenshot(self, fname='XY.png', format='png'):
        p = self.grab();
        fname2 = ('pics/' + 
                  str(self.molecule.moleculeNumber) +
                  self.coord1 +
                  self.coord2 +
                  '.png')
        p.save(fname2, format, -1)

def main():
    app = QApplication(sys.argv)

    fname = '/home/anton/DataExamples/10x20.data'
    if SYSTEM_NAME == '10x20':
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
        bondsList.append(bo)
    i = 0
    if PHASE == 'polymer': # for polymer
        l = 90
    elif PHASE == 'modifier': # for modifier
        l = 108
    for i in range(l): # for modifier
        print('Preparing molecules: ', i + 1, ' / ', moleculeNumber)
        mol = Molecule(i)
        mol.makeSkeleton(atomsList, bounds)
        mol.chooseBonds(atomsList, bondsList)
        mol.computeRanges()

        w = MainWidget(mol, 'X', 'Y')
        w = MainWidget(mol, 'X', 'Z')
        w = MainWidget(mol, 'Y', 'Z')

    app.quit()

main()
