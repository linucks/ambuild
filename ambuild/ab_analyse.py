import csv
import time


class Analyse():
    def __init__(self, cell, logfile="ambuild.csv"):

        self.fieldnames = ['step',
                           'type',
                           'tot_time',
                           'time',
                           'num_frags',
                           'num_particles',
                           'num_blocks',
                           'density',
                           'num_free_endGroups',
                           'potential_energy',
                           'num_tries',
                           'fragment_types',
                           'file_count',
                            ]
        self.cell = cell

        self.step = 0
        self._startTime = time.time()
        self._stepTime = time.time()

        # Need to create an initial entry as we query the previous one for any data we don't have
        d = {}
        for f in self.fieldnames:
            if f == 'time':
                d[ f ] = self._startTime
            # elif f == 'file_count':
            elif f == 'type':
                d[ f ] = 'init'
            else:
                d[ f ] = 0

        self.last = d

        self.logfile = logfile
        self._logWriter = csv.DictWriter(open(self.logfile, 'w'), self.fieldnames)

        self._logWriter.writeheader()

        return

    def start(self):
        """Called whenever we start a step"""
        assert self._stepTime == None
        assert self.last
        self.step += 1
        self._stepTime = time.time()
        return

    def stop(self, stype, d={}):
        """Called at the end of a step with the data to write to the csv file"""

        new = {}

        for f in self.fieldnames:
            if f == 'type':
                new[ f ] = stype
            elif f == 'step':
                new [ f ] = self.step
            elif f == 'time':
                new[ f ] = time.time() - self._stepTime
            elif f == 'tot_time':
                new[ f ] = time.time() - self._startTime
                self._stepTime
            elif f == 'num_blocks':
                new[ f ] = self.cell.numBlocks()
            elif f == 'num_frags':
                new[ f ] = self.cell.numFragments()
            elif f == 'num_particles':
                new[ f ] = self.cell.numAtoms()
            elif f == 'num_free_endGroups':
                new[ f ] = self.cell.numFreeEndGroups()
            elif f == 'density':
                new[ f ] = self.cell.density()
            elif f == 'fragment_types':
                new[ f ] = str(self.cell.fragmentTypes())
            elif f == 'file_count':
                new[ f ] = self.cell._fileCount
            elif f in d:
                new[ f ] = d[ f ]
            else:
                new[ f ] = self.last[ f ]

        self._logWriter.writerow(new)

        self.last = new
        self._stepTime = None
        self.start()
        return
