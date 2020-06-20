class Bond(object):
    """An object to hold all info on a bond
    """

    def __init__(self, endGroup1, endGroup2):
        self.endGroup1 = endGroup1
        self.endGroup2 = endGroup2
        return

    def engage(self):
        assert self.endGroup1 != self.endGroup2
        assert self.isPossible()
        assert not self.endGroup1.fragment == self.endGroup2.fragment
        self.endGroup1.setBonded(self)
        self.endGroup2.setBonded(self)
        self.endGroup1.block().bondBlock(self)

    def isInternalBond(self):
        return self.endGroup1.block() == self.endGroup2.block()

    def isPossible(self):
        return self.endGroup1.free() and self.endGroup2.free()

    @property
    def rootFragment(self):
        return self.endGroup1.fragment

    @property
    def rootId(self):
        return self.endGroup1.block().id

    def separate(self):
        self.endGroup1.unBond(self.endGroup2)
        self.endGroup2.unBond(self.endGroup1)

    @property
    def targetId(self):
        return self.endGroup2.block().id

    @property
    def targetFragment(self):
        return self.endGroup2.fragment

    def __str__(self):
        """List the data attributes of this object"""
        return "Bond {0}: {1}:{2} -> {3}:{4}".format(
            id(self), self.rootId, self.endGroup1, self.targetId, self.endGroup2
        )
