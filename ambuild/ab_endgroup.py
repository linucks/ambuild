'''
Created on Jan 15, 2013

@author: abbietrewin
'''
# import collections
# import copy
# import csv
import logging
# import os
# 
import numpy as np

ENDGROUPBONDED = '*'
logger = logging.getLogger()

class EndGroup(object):

    def __init__(self):

        self.blocked = False # Used for indicating that this endGroup is unbonded but not free
        self.bonded = False
        self._endGroupType = None
        self.fragment = None
        self.fragmentEndGroupIdx = None
        self.blockEndGroupIdx = None
        self.fragmentCapIdx = None
        self.blockCapIdx = None
        self.capBondLength = None
        self.fragmentDihedralIdx = -1
        self.blockDihedralIdx = -1
        self.fragmentUwIdx = -1
        self.blockUwIdx = -1
        return

    def block(self):
        f = True
        b = True
        if self.fragment.block is None:
            b = False
        if self.fragment is None:
            f = False
        if not b and f:
            raise RuntimeError("None Block {0} Fragment {1}\n{2}".format(b, f, self))
        return self.fragment.block

    def capIdx(self):
        """Return the index of the endGroup atom in external block indices"""
        return self.blockCapIdx

    def bondedCatalyst(self):
        """Return True if this endGroup belongs to a catalyst that is bonded to another catalyst"""
        return self.fragment.catalyst and self._endGroupType.endswith(ENDGROUPBONDED)

    def coord(self,endGroup=True):
        """Need to think about an API for accessing coordinates for endGroups
        This just hacks in returning the endGroup.
        """
        return self.fragment._coords[self.fragmentEndGroupIdx]

    def dihedralIdx(self):
        """Return the index of the dihedral atom in external block indices"""
        return self.blockDihedralIdx

    def endGroupIdx(self):
        """Return the index of the endGroup atom in external block indices"""
        return self.blockEndGroupIdx

    def fragmentType(self):
        return self.fragment.fragmentType

    def free(self):
        return not self.bonded and not self.blocked

    def setBonded(self, bond):
        self.bonded = True
        self.fragment.addBond(self, bond)
        return

    def unBond(self, bondEndGroup):
        self.bonded = False
        # HACK WE REMOVE ALL SUFFIXES
        for eg in self.fragment.endGroups():
            if eg._endGroupType.endswith(ENDGROUPBONDED):
                logger.debug("unBond ENDGROUPBONDED")
                eg._endGroupType = eg._endGroupType.rstrip(ENDGROUPBONDED)
        self.fragment.delBond(self.type())
        # Unmask cap and set the coordinate to the coordinate of the last block atom
        # NOTE - NEED TO SCALE BY CORRECT LENGTH
        if hasattr(bondEndGroup, 'coord'):
            # Reposition the cap atom based on the bond vector
            #self.fragment._coords[self.fragmentCapIdx] = bondEndGroup.coord()
            # Get vector from this endGroup to the other endGroup
            egPos = self.fragment._coords[self.fragmentEndGroupIdx]
            v1 = bondEndGroup.coord() - egPos
            # Now get unit vector
            uv = v1 / np.linalg.norm(v1)
            # calculate noew position
            self.fragment._coords[self.fragmentCapIdx] = egPos + (uv * self.capBondLength)

        # Unhide the cap atom
        self.fragment.masked[self.fragmentCapIdx] = False
        # Mark capAtom as unBonded so that it won't be included in the optimisation
        self.fragment.unBonded[self.fragmentCapIdx] = True

        if self.fragmentUwIdx != -1:
            raise RuntimeError("Cannot unbond masked endGroups yet!")
            self.fragment.masked[ self.fragmentUwIdx ] = True
        self.fragment.update()
        return

    def type(self):
        """The type of the endGroup"""
        return "{0}:{1}".format(self.fragment.fragmentType, self._endGroupType)

    def updateAncillaryIndices(self, cap2endGroup):
        """The block has been updated so we need to update our block indices based on where the
        fragment starts in the block"""
        if self.fragment.masked[self.fragmentCapIdx]:
            assert self.bonded
            # If this endGroup is involved in a bond, we want to get the block index of the
            # opposite endGroup as this has now become our cap atom
            self.blockCapIdx = cap2endGroup[(self.fragment, self.fragmentCapIdx)]
        else:
            self.blockCapIdx = self.fragment._int2ext[self.fragmentCapIdx] + self.fragment.blockIdx

        # -1 means no dihedral or uw atom set
        if self.fragmentDihedralIdx != -1:
            if self.fragment.masked[self.fragmentDihedralIdx]:
                # Now work out which atom this is bonded to
                self.blockDihedralIdx = cap2endGroup[(self.fragment, self.fragmentDihedralIdx)]
            else:
                self.blockDihedralIdx = self.fragment._int2ext[self.fragmentDihedralIdx] + self.fragment.blockIdx

        # -1 means no uw atom and if it is masked then there will not be an external index
        if self.fragmentUwIdx != -1 and not self.fragment.masked[self.fragmentUwIdx]:
            # assert not self.fragment.isMasked( self.fragmentUwIdx )
            self.blockUwIdx = self.fragment._int2ext[self.fragmentUwIdx] + self.fragment.blockIdx
        return

    def updateEndGroupIndex(self):
        """The block has been updated so we need to update our block indices based on where the
        fragment starts in the block"""
        self.blockEndGroupIdx = self.fragment._int2ext[self.fragmentEndGroupIdx] + self.fragment.blockIdx
        return

    def __str__(self):
        """List the data attributes of this object"""
#         me = {}
#         for slot in dir(self):
#             attr = getattr(self, slot)
#             if not slot.startswith("__") and not ( isinstance(attr, types.MethodType) or
#               isinstance(attr, types.FunctionType) ):
#                 me[slot] = attr
#
#         t = []
#         for k in sorted(me):
#             t.append( str( ( k, me[k] ) ) )
#
#         return "{0} : {1}".format(self.__repr__(), ",".join( t ) )
        #EndGroup:
        s = "{0} {1}:{2} {3}({4})->{5}({6}) {7}->{8}".format(self.type(), id(self.fragment), id(self),
                                                    self.fragmentEndGroupIdx, self.fragment._symbols[self.fragmentEndGroupIdx],
                                                    self.fragmentCapIdx,self.fragment._symbols[self.fragmentCapIdx],
                                                    self.blockEndGroupIdx, self.blockCapIdx)

        return s
