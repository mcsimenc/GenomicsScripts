#!/usr/bin/env python3

import sys


def help():
    print("""
    Usage:
    ------------
    gff2circos-heatmap.py -gff <path> -window <int> 
        [<path> -scafLens | -fasta <path>] [-scafList <path>]

    Description:
    ------------
    Takes two files as input, a GFF3 file and either (1) a two-column 
    tab-delimited file where the first column contains scaffold names 
    and the second column contains scaffold lengths, or (2) a FASTA
    file from which sequence lengths will be obtained. The -window 
    parameter sets the number of bases in the sliding window used to 
    determine the density of features and corresponding colors in
    Circos. Different GFF3 features types are not distinguished from
    each other; therefore the input GFF3 file should contain only one
    type of feature (e.g. genes of interest). Overlapping features are
    merged and feature strand is not taken into account. Features are
    all treated as if they are on a single strand.

    Required parameters:
    ------------
    -gff          <path>    Path to input maker gff3 file (mandatory)
    -scafLens     <path>    Scaffold lengths. 2-column, (1) scaf name 
                            and (2) scaffold lengths
    -fasta        <path>    FASTA file corresponding to the seq names 
                            in the -gff Alternative to using -scafLens
    -window       <int>     The number of bases in the sliding window 
    -scafList     <path>    Restrict output to these scaffolds. A file
                            with a list of scaffold names
    Output:
    ------------
    scaf    start    stop    density
""", file=sys.stderr)


def parseLenFile(filepath):
    """Parses two-column input file returns dict with first column item 
    as key and second column item as value. Keys are expected to be
    unique.
    """
    lenDct = {}
    with open(filepath) as fl:
        for line in fl:
            contents = line.strip().split('\t')
            lenDct[contents[0]] = int(contents[1])

    return lenDct

def fasta2LenDct(filepath, scafList=None):
    """
    Calculates the lengths of seqs in fasta format file provided with 
    -fasta and returns a dictionary with seq headers as keys and 
    lengths as values. Doesn't do any checks that the FASTA file is
    the correct format. Whitespace in sequence headers will not work
    in Circos.

    Arguments:
    ----------
    filepath: path to a FASTA file
    scafList: optional list of scaffolds to restrict output to
    """
    # initialize variables
    lenDct = {}
    header = None
    removedLastScaf = False
    if scafList:
        scafList = set(scafList)
    # open the fasta file and read it line by line, parsing sequence
    # headers and then counting the number of bases in the sequence
    with open(filepath) as fl:
        for line in fl:
            # parse sequence headers
            if line.startswith('>'):
                header = line.strip()[1:]
                removedLastScaf = False
            # a list of scaffolds has been provided. remove the current
            # sequence header from the set of sequence headers of which
            # to restrict output. stop parsing the fasta file if the
            # set of sequence headers is empty
            elif header and scafList != None and not removedLastScaf:
                if header not in scafList:
                    if scafList == set():
                        break
                    continue
                else:
                    scafList.remove(header)
                    removedLastScaf = True
            # add the number of characters in each line of sequence to
            # the running total
            else:
                if header in lenDct:
                    lenDct[header] += len(line.strip())
                else:
                    lenDct[header] = len(line.strip())
    # return a dictionary: {header:length}
    return lenDct

            

def gff2circosHeatmap(filepath, 
                      scafLens, 
                      windowLen, 
                      scafList=None):
    """1. Reads scaf, start, and start for features in GFF3 file into 
       memory
    2. Creates all windows for a given scaf
    3. Loops through the features and adds them to whatever scaf they
       are part of

    Arguments:
    ----------
    filepath: Path to a GFF3 file
    scafLens: Dictionary with lengths of sequences represented in the
              GFF3 file
    """

    def whichWindow(coord, windowLen):
        """Takes a coordinate and the length of the window as arguments
        and returns the window number.
        """
        return coord//windowLen

    def genWindowCoords(n, windowLen):
        """Takes the window number n and the length of the window as
        arguments and returns a 2-element tuple of the start and end
        coordinates of that window.
        """
        return (n * windowLen, (n + 1) * windowLen - 1)

    def mergeCoords(A, B):
        """Takes two tuples of coordinates A and B as arguments; A must
        have a start coordinate before or equal to the start coordinate
        of B. If they do not overlap then A and B are returned as input
        and if they overlap the minimum and maximum values are returned.

        let A = (a1, a2), B = (b1, b2) | a1<=b1, a1<=a2, b1<=b2

        case 1: a2<=b1 ---> output A and B

        case 2: b1<a2 && b2>a2 ---> output (a1, b2)

        case 3: b2<=a2 ---> output A
        """

        assert min(A) <= min(B), ("tuples given to mergeCoords in wrong order: "
                                    "A={0}, B={1}").format(A,B)
        if min(B) >= max(A):
            return ((A, B), 0)
        elif min(B) < max(A) and max(B) > max(A):
            output = (min(A),max(B))
            return ((output, output), 1)
        elif max(B) <= max(A):
            return ((A, A), 2)
        else:
            raise Exception(("Unexpected result from mergeCoords(A,B) using "
                               " A={0}, B={1}").format(A,B))

    # open gff file and read it line by line, adding start and end
    # coordinates for each feature to the scaffold on which it resides.
    with open(filepath) as fl:
        gffFeats = {}
        for line in fl:
            # skip commented lines
            if not line.startswith('#'):
                contents = line.strip().split('\t')
                scaf = contents[0]
                start = int(contents[3]) - 1
                end = int(contents[4]) - 1
                if scaf in gffFeats:
                    gffFeats[scaf].append((start, end))
                else:
                    gffFeats[scaf] = [(start, end)]
        # for each scaffold, sort features by length then merge
        # all overlapping adjacent features
        for scaf in gffFeats:
            gffFeats[scaf] = sorted(gffFeats[scaf], key=lambda x:x[0])
            newScafFeats = []
            currentFeat = gffFeats[scaf][0]
            i = 0
            while i < len(gffFeats[scaf]) - 1:
                mergeResult = mergeCoords(currentFeat, gffFeats[scaf][i + 1])
                # adjacent features do not overlap: save the first as
                # a feature in the processed merged set and set the
                # second as the next feature to check for overlaps
                if mergeResult[1] == 0:
                    newScafFeats.append(mergeResult[0][0])
                    currentFeat = mergeResult[0][1]
                # features overlap: merge them and continue to check
                # if they overlap with the next adjacent feature
                elif mergeResult[1] == 1:
                    currentFeat = mergeResult[0][0]
                elif mergeResult[1] == 2:
                    currentFeat = mergeResult[0][0]
                else:
                    raise Exception("Unexpected result from mergeResult block")
                i += 1
            newScafFeats.append(currentFeat)
            gffFeats[scaf] = newScafFeats
        # if a scaffold list to restrict output to hasn't been provided
        # then set the scaffolds to output as the scaffolds with
        # features
        if scafList == None:
            scafList = sorted(list(gffFeats.keys()))
        # for each scaffold create a set of windows, the last of which will be
        # the length from the end of the second to last window to the end of
        # the scaffold
        for scaf in scafList:
            windows = {}
            # for cumulative length of all features in window
            windowFeatLens = {}
            lastWindowEnd = 0
            currentWindow = 0
            # generate windows for each scaffold
            while lastWindowEnd < scafLens[scaf] - windowLen + 1:
                windows[currentWindow] = genWindowCoords(currentWindow, windowLen)
                windowFeatLens[currentWindow] = 0
                lastWindowEnd = windows[currentWindow][1]
                currentWindow += 1
            assert lastWindowEnd <= scafLens[scaf], ("During window generation "
                                  "window extends beyond length of scaffold.")
            # generate the window for the last scaffold (may be shorter than
            # the window size)
            if lastWindowEnd < scafLens[scaf]:
                windows[currentWindow] = (lastWindowEnd + 1, scafLens[scaf])
                windowFeatLens[currentWindow] = 0
            # calculate proportion of window occupied by any feature by summing
            # the lengths of the each feature. This is output as the density
            # 4th column in the Circos heatmap file
            if scaf in gffFeats:
                for feature in gffFeats[scaf]:
                    startWindow = whichWindow(feature[0], windowLen)
                    endWindow = whichWindow(feature[1], windowLen)
                    windowSpread = endWindow - startWindow
                    assert windowSpread >= 0, ("windowSpread < 0, check "
                                                "whichWhindow() implementation")
                    # feature is fully contained in one window: add the length
                    # of this feature to the running total
                    if windowSpread == 0:
                        windowFeatLens[startWindow] += (feature[1] - feature[0])
                        assert windowFeatLens[startWindow] <= windowLen, ("1: "
                                   "window coverage greater than window length")
                    # feature is spread across two or more winodws: add the
                    # length of the fraction of this feature within each
                    # window to the running totals of each window
                    elif windowSpread >= 1:
                        windowFeatLens[startWindow] += (windows[startWindow][1] 
                                                        - feature[0])
                        windowFeatLens[endWindow] += (feature[1] 
                                                      -  windows[endWindow][0])
                        assert windowFeatLens[startWindow] <= windowLen, ("2.start: "
                                  "window coverage greater than window length")
                        assert windowFeatLens[endWindow] <= windowLen, ("2.end:"
                                 " window coverage greater than window length")
                        # feature is spread across three or more windows:
                        # for each internal window mark the density as 100%
                        if windowSpread > 1:
                            for i in range(startWindow+1, endWindow):
                                windowFeatLens[i] += windowLen
                                assert windowFeatLens[i] <= windowLen, ("3. "
                                  "window coverage greater than window length")
                # this scaffold contains features. output Circos heatmap line
                for i in range(len(windowFeatLens)):
                    propFeatInWindow = windowFeatLens[i]/windowLen
                    windowStart = windows[i][0]
                    windowEnd = windows[i][1]
                    print('{0}\t{1}\t{2}\t{3:.10f}'.format(scaf, windowStart, 
                                                  windowEnd, propFeatInWindow))
            # this scaffold contains no features. output lines with 0 for every
            # value in the density column
            else:
                for i in windows:
                    print('{0}\t{1}\t{2}\t{3:.10f}'.format(scaf, windows[i][0], 
                                                             windows[i][1], 0))
            for i in windowFeatLens:
                assert windowFeatLens[i] <= windowLen, ("4. window coverage"
                                                 " greater than window length")


if __name__ == '__main__':
    args = sys.argv
    # print help if command line arguments are missing
    if ('-h' in args 
     or len(args) < 7 
     or ('-gff' not in args 
        or ('-scafLens' not in args and '-fasta' not in args) 
        or '-window' not in args)):
        help()
        sys.exit()
    # read in gff file
    gff_filepath = args[args.index('-gff') +1]
    # read in list of scaffolds to restrict output to
    if '-scafList' in args:
        scafListFl = args[args.index('-scafList') +1]
        scafList = open(scafListFl).read().strip().split('\n')
    else:
        scafList = None
    # calculate scaffold lengths from fasta file
    if '-fasta' in args:
        fasta_filepath = args[args.index('-fasta') +1]
        scafLens = fasta2LenDct(fasta_filepath, scafList)
    # parse file containing scaffold lengths
    elif '-scafLens' in args:
        scafLens_filepath = args[args.index('-scafLens') +1]
        scafLens = parseLenFile(scafLens_filepath)
    # read the window length
    windowLen = int(args[args.index('-window') +1])
    # calculate and output the Circos heatmap
    gff2circosHeatmap(gff_filepath, scafLens, windowLen, scafList)
