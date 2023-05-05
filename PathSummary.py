__author__ = 'marti'
from sys import maxint
from math import ceil
import settings

class PathSummary:
    def __init__(self):
        self.pathsComputeTime = 0
        self.pathsNumber = 0
        self.noOfInvariants = 0
        self.maxPathLong = 0
        self.minPathLong = maxint
        self.avrgPathLong = 0
        self.percentilPath20 = 0
        self.percentilPath50 = 0
        self.percentilPath70 = 0
        self.percentilPath90 = 0
        self.pathsLong = []
        self.paths_limit_raised = False

    def to_list(self):
        ret = []
        if not self.paths_limit_raised:
            ret.append('%.5f s' % self.pathsComputeTime)
            ret.append('%i' % self.pathsNumber)
        else:
            ret.append('MAXPATHSRAISED')
            ret.append('>' + str(settings.MAX_PATHS))
        ret.append('%i' % self.noOfInvariants)
        ret.append('%i' % self.maxPathLong)
        if self.minPathLong == maxint:
            ret.append("NA")
        else:
            ret.append('%i' % self.minPathLong)
        ret.append('%.3f' % self.avrgPathLong)
        ret.append('%i' % self.percentilPath20)
        ret.append('%i' % self.percentilPath50)
        ret.append('%i' % self.percentilPath70)
        ret.append('%i' % self.percentilPath90)
        return ret


def percentilePreSorted(data, percentile):
    size = len(data)
    return data[int(ceil((size * percentile) / 100)) - 1]
