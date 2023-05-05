__author__ = 'marti'
import urllib2
import csv


def getCompoundInfo(compound, compFile):
    try:
        print 'Downloading compound info for ' + compound + '...'
        url = urllib2.urlopen(str("http://www.genome.jp/dbget-bin/www_bget?-f+k+compound+" + compound))
        content = url.read()
        if len(content) < 2:
            return False
        with open(compFile, "a") as f:
            f.write(content)
        return True
    except:
        print "Error getting " + compound + " info. Please note that you need an active internet connection for this" \
                                            "process. Plesase check that and try again. Aborting process..."
        return False

def loadFromCSV(csvFile, delimiterChar, columnNumber, rowsToSkip, max):
    ret = []
    with open(csvFile, 'rb') as csvFile:
        reader = csv.reader(csvFile, delimiter=delimiterChar, quotechar='"', quoting=csv.QUOTE_ALL)
        for x in xrange(0, rowsToSkip):
            next(reader)
        for row in reader:
            ret.append(str(row[columnNumber]))
            max -= 1
            if max == 0:
                break
    return ret

