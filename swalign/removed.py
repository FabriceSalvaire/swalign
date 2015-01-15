####################################################################################################

import StringIO

####################################################################################################

class ScoringMatrix(object):

    '''
    Read scoring matrix from a file or string

    Matrix should be space-delimited in a format like:

      A C G T
    A 1 0 0 0
    C 0 1 0 0
    G 0 0 1 0
    T 0 0 0 1

    Rows and Columns must be in the same order

    '''

    ##############################################

    def __init__(self, filename=None, text=None, wildcard_score=0):

        assert filename or text

        if filename:
            fs = open(filename)
        else:
            fs = StringIO.StringIO(text)

        self.scores = []
        self.bases = None
        self.wildcard_score = wildcard_score

        for line in fs:
            if line[0] == '#':
                continue

            if not self.bases:
                self.bases = line.split()
                self.base_count = len(self.bases)
            else:
                cols = line.split()
                self.scores.extend([int(x) for x in cols[1:]])

        fs.close()

    ##############################################

    def score(self, one, two, wildcard=None):

        if self.wildcard_score and wildcard and (one in wildcard or two in wildcard):
            return self.wildcard_score

        one_idx = 0
        two_idx = 0
        for i, b in enumerate(self.bases):
            if b == one:
                one_idx = i
            if b == two:
                two_idx = i

        return self.scores[(one_idx * self.base_count) + two_idx]

####################################################################################################

def fasta_gen(fname):

    def gen():
        seq = ''
        name = ''
        comments = ''

        if fname == '-':
            f = sys.stdin
            name = 'stdin'
        else:
            f = open(fname)

        for line in f:
            if line[0] == '>':
                if name and seq:
                    yield (name, seq, comments)

                spl = line[1:].strip().split(' ', 1)
                name = spl[0]
                if len(spl) > 1:
                    comments = spl[1]
                else:
                    comments = ''

                seq = ''
            else:
                seq += line.strip()

        if name and seq:
            yield (name, seq, comments)

        if fname != '-':
            f.close()

    return gen

####################################################################################################

def seq_gen(name, seq):

    def gen():
        yield (name, seq, '')

    return gen

####################################################################################################

def extract_region(comments):

    ref = None
    start = None
    # start_offset = 0
    # end_offset = 0

    try:
        attrs = comments.split(' ')
        for attr in attrs:
            if '=' in attr:
                k, v = attr.split('=')
                if k == 'range':
                    spl = v.split(':')
                    ref = spl[0]
                    start, end = [int(x) for x in spl[1].split('-')]
                # elif k == "5'pad":
                #     start_offset = int(v)
                # elif k == "3'pad":
                #     end_offset = int(v)
    except:
        pass

    if ref and start:
        return (ref, start - 1, '%s:%s-%s' % (ref, start, end))

    return None

####################################################################################################

__revcomp = {}
for a, b in zip('atcgATCGNn', 'tagcTAGCNn'):
    __revcomp[a] = b
__cache = {}

####################################################################################################

def revcomp(seq):

    if seq in __cache:
        return __cache[seq]

    ret = []
    for s in seq.upper()[::-1]:
        ret.append(__revcomp[s])

    __cache[seq] = ''.join(ret)

    return __cache[seq]

####################################################################################################
# 
# End
# 
####################################################################################################
