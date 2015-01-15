####################################################################################################

''' Simple Smith-Waterman aligner '''

####################################################################################################

import sys

####################################################################################################

class IdentityScoringMatrix(object):

    ##############################################

    def __init__(self, match=1, mismatch=-1):

        self._match = match
        self._mismatch = mismatch

    ##############################################

    def score(self, one, two, wildcard=None):

        if wildcard and (one in wildcard or two in wildcard):
            return self._match
        else:
            if one == two:
                return self._match
            else:
                return self._mismatch

####################################################################################################

NucleotideScoringMatrix = IdentityScoringMatrix

####################################################################################################

class MatrixCell(object):

    ##############################################

    def __init__(self, score=0, op=' ', run_length=0):

        self.score = score
        self.op = op
        self.run_length = run_length

    ##############################################

    def set(self, other):

        self.score = other.score
        self.op = other.op
        self.run_length = other.run_length

####################################################################################################

class Matrix(object):

    ##############################################

    def __init__(self, number_of_rows, number_of_cols):

        self.number_of_rows = number_of_rows
        self.number_of_cols = number_of_cols
        self._values = [[MatrixCell() for c in xrange(number_of_cols)]
                        for r in xrange(number_of_rows)]

    ##############################################

    def get(self, row, col):
        return self._values[row][col]

    ##############################################

    def set(self, row, col, value):
        self._values[row][col].set(value)

####################################################################################################

class LocalAlignment(object):

    ##############################################

    def __init__(self, scoring_matrix,
                 gap_penalty=-1, gap_extension_penalty=-1, gap_extension_decay=0.0,
                 prefer_gap_runs=True,
                 verbose=False, wildcard=None):

        self.scoring_matrix = scoring_matrix
        self.gap_penalty = gap_penalty
        self.gap_extension_penalty = gap_extension_penalty
        self.gap_extension_decay = gap_extension_decay
        self.verbose = verbose
        self.prefer_gap_runs = prefer_gap_runs
        self.wildcard = wildcard

    ##############################################

    def align(self, ref, query, ref_name='ref', query_name='query', rc=False):

        orig_ref = ref
        orig_query = query

        ref = ref.upper()
        query = query.upper()

        matrix = Matrix(len(query) +1, len(ref) +1)
        for row in xrange(1, matrix.number_of_rows):
            matrix.set(row, 0, MatrixCell(0, 'i', 0))
        for col in xrange(1, matrix.number_of_cols):
            matrix.set(0, col, MatrixCell(0, 'd', 0))

        max_value = 0
        max_row = 0
        max_col = 0

        # calculate matrix
        for row in xrange(1, matrix.number_of_rows):
            for col in xrange(1, matrix.number_of_cols):

                left = matrix.get(row, col -1) # Deletion
                up = matrix.get(row -1, col) # Insertion
                previous = matrix.get(row -1, col -1) # Match/Mismatch

                # Match/Mismatch
                mm_value = previous.score + self.scoring_matrix.score(query[row -1], ref[col -1], self.wildcard)

                ins_run = 0
                del_run = 0

                if up.op == 'i':
                    ins_run = up.run_length
                    if up.score == 0:
                        # no penalty to start the alignment
                        ins_value = 0
                    else:
                        # if not self.gap_extension_decay:
                        ins_value = up.score + self.gap_extension_penalty
                        # else:
                        #     ins_value = (up.score +
                        #                  min(0, self.gap_extension_penalty + ins_run * self.gap_extension_decay))
                else:
                    ins_value = up.score + self.gap_penalty

                if left.op == 'd':
                    del_run = left.run_length
                    if left.score == 0:
                        # no penalty to start the alignment
                        del_value = 0
                    else:
                        # if not self.gap_extension_decay:
                        del_value = left.score + self.gap_extension_penalty
                        # else:
                        #     del_value = (left.score +
                        #                  min(0, self.gap_extension_penalty + del_run * self.gap_extension_decay))

                else:
                    del_value = left.score + self.gap_penalty
                    
                cell_value = max(mm_value, del_value, ins_value, 0)

                if not self.prefer_gap_runs:
                    # clear run length
                    ins_run = 0
                    del_run = 0

                if del_run and cell_value == del_value:
                    cell = MatrixCell(cell_value, 'd', del_run +1)
                elif ins_run and cell_value == ins_value:
                    cell = MatrixCell(cell_value, 'i', ins_run +1)
                elif cell_value == mm_value:
                    cell = MatrixCell(cell_value, 'm', 0)
                # prefer_gap_runs
                elif cell_value == del_value:
                    cell = MatrixCell(cell_value, 'd', 1)
                elif cell_value == ins_value:
                    cell = MatrixCell(cell_value, 'i', 1)
                else:
                    # ???
                    cell = MatrixCell(0, 'x', 0)

                if cell.score >= max_value:
                    max_value = cell.score
                    max_row = row
                    max_col = col

                matrix.set(row, col, cell)

        # backtrack
        row = max_row
        col = max_col
        value = max_value

        op = ''
        aln = []

        path = []
        while True:
            item = matrix.get(row, col)

            if item.score <= 0:
                break

            path.append((row, col))
            aln.append(item.op)

            if item.op == 'm':
                row -= 1
                col -= 1
            elif item.op == 'i':
                row -= 1
            elif item.op == 'd':
                col -= 1
            else:
                break

        aln.reverse()

        if self.verbose:
            print '-'*80
            self.dump_matrix(ref, query, matrix, path)
            print
            print aln
            print
            print 'max:', (max_row, max_col), max_value
            print '-'*80

        cigar = _reduce_cigar(aln)

        return Alignment(orig_query, orig_ref, row, col, cigar, max_value,
                         ref_name, query_name, rc, self.wildcard)

    ##############################################

    def dump_matrix(self, ref, query, matrix, path, show_row=-1, show_col=-1):

        output = sys.stdout

        output.write('      -      ')
        output.write('       '.join(ref))
        output.write('\n')
        for row in xrange(matrix.number_of_rows):
            if row == 0:
                output.write('-')
            else:
                output.write(query[row -1])

            for col in xrange(matrix.number_of_cols):
                if show_row == row and show_col == col:
                    output.write('       *')
                else:
                    output.write(' %5s%s%s' % (matrix.get(row, col).score, matrix.get(row, col).op,
                                               '$' if (row, col) in path else ' '))
            output.write('\n')

####################################################################################################

def _reduce_cigar(operations):

    count = 1
    last = None
    ret = []
    for op in operations:
        if last and op == last:
            count += 1
        elif last:
            ret.append((count, last.upper()))
            count = 1
        last = op

    if last:
        ret.append((count, last.upper()))
    return ret

####################################################################################################

def _cigar_str(cigar):

    out = ''
    for num, op in cigar:
        out += '%s%s' % (num, op)
    return out

####################################################################################################

class Alignment(object):

    ##############################################

    def __init__(self, query, ref, q_pos, r_pos, cigar, score, ref_name='', query_name='',
                 rc=False, wildcard=None):

        self.query = query
        self.ref = ref
        self.q_pos = q_pos
        self.r_pos = r_pos
        self.cigar = cigar
        self.score = score
        self.r_name = ref_name
        self.q_name = query_name
        self.rc = rc
        self.wildcard = wildcard

        self.r_offset = 0
        self.r_region = None

        self.orig_query = query
        self.query = query.upper()

        self.orig_ref = ref
        self.ref = ref.upper()

        q_len = 0
        r_len = 0

        self.matches = 0
        self.mismatches = 0

        i = self.r_pos
        j = self.q_pos

        for count, op in self.cigar:
            if op == 'M':
                q_len += count
                r_len += count
                for k in xrange(count):
                    if self.query[j] == self.ref[i]:
                        self.matches += 1
                    else:
                        self.mismatches += 1
                    i += 1
                    j += 1

            elif op == 'I':
                q_len += count
                j += count
                self.mismatches += count
            elif op == 'D':
                r_len += count
                i += count
                self.mismatches += count

        self.q_end = q_pos + q_len
        self.r_end = r_pos + r_len
        if self.mismatches + self.matches > 0:
            self.identity = float(self.matches) / (self.mismatches + self.matches)
        else:
            self.identity = 0

    ##############################################

    def set_ref_offset(self, ref, offset, region):

        self.r_name = ref
        self.r_offset = offset
        self.r_region = region

    ##############################################

    @property
    def extended_cigar_str(self):

        qpos = 0
        rpos = 0
        ext_cigar_str = ''
        working = []
        for count, op in self.cigar:
            if op == 'M':
                for k in xrange(count):
                    if self.query[self.q_pos + qpos + k] == self.ref[self.r_pos + rpos + k]:
                        ext_cigar_str += 'M'
                    else:
                        ext_cigar_str += 'X'
                qpos += count
                rpos += count

            elif op == 'I':
                qpos += count
                ext_cigar_str += 'I' * count
            elif op == 'D':
                rpos += count
                ext_cigar_str += 'D' * count

            working = _reduce_cigar(ext_cigar_str)

        out = ''
        for num, op in working:
            out += '%s%s' % (num, op)
        return out

    ##############################################

    @property
    def cigar_str(self):
        return _cigar_str(self.cigar)

    ##############################################

    def dump(self, wrap=None, out=sys.stdout):

        i = self.r_pos
        j = self.q_pos

        q = ''
        m = ''
        r = ''
        qlen = 0
        rlen = 0

        for count, op in self.cigar:
            if op == 'M':
                qlen += count
                rlen += count
                for k in xrange(count):
                    q += self.orig_query[j]
                    r += self.orig_ref[i]
                    if self.query[j] == self.ref[i] or (self.wildcard and
                                                        (self.query[j] in self.wildcard or
                                                         self.ref[i] in self.wildcard)):
                        m += '|'
                    else:
                        m += '.'

                    i += 1
                    j += 1
            elif op == 'D':
                rlen += count
                for k in xrange(count):
                    q += '-'
                    r += self.orig_ref[i]
                    m += ' '
                    i += 1
            elif op == 'I':
                qlen += count
                for k in xrange(count):
                    q += self.orig_query[j]
                    r += '-'
                    m += ' '
                    j += 1

            elif op == 'N':
                q += '-//-'
                r += '-//-'
                m += '    '

        if self.q_name:
            out.write('Query: %s%s (%s nt)\n' % (self.q_name, ' (reverse-compliment)'
                                                 if self.rc else '', len(self.query)))
        if self.r_name:
            if self.r_region:
                out.write('Ref  : %s (%s)\n\n' % (self.r_name, self.r_region))
            else:
                out.write('Ref  : %s (%s nt)\n\n' % (self.r_name, len(self.ref)))

        poslens = [self.q_pos +1, self.q_end +1,
                   self.r_pos + self.r_offset +1,
                   self.r_end + self.r_offset +1]
        maxlen = max([len(str(x)) for x in poslens])

        q_pre = 'Query: %%%ss ' % maxlen
        r_pre = 'Ref  : %%%ss ' % maxlen
        m_pre = ' ' * (8 + maxlen)

        rpos = self.r_pos
        if not self.rc:
            qpos = self.q_pos
        else:
            qpos = self.q_end

        while q and r and m:
            if not self.rc:
                out.write(q_pre % (qpos +1))  # pos is displayed as 1-based
            else:
                out.write(q_pre % (qpos))  # revcomp is 1-based on the 3' end

            if wrap:
                qfragment = q[:wrap]
                mfragment = m[:wrap]
                rfragment = r[:wrap]

                q = q[wrap:]
                m = m[wrap:]
                r = r[wrap:]
            else:
                qfragment = q
                mfragment = m
                rfragment = r

                q = ''
                m = ''
                r = ''

            out.write(qfragment)
            if not self.rc:
                for base in qfragment:
                    if base != '-':
                        qpos += 1
            else:
                for base in qfragment:
                    if base != '-':
                        qpos -= 1

            if not self.rc:
                out.write(' %s\n' % qpos)
            else:
                out.write(' %s\n' % (qpos +1))

            out.write(m_pre)
            out.write(mfragment)
            out.write('\n')
            out.write(r_pre % (rpos + self.r_offset +1))
            out.write(rfragment)
            for base in rfragment:
                if base != '-':
                    rpos += 1
            out.write(' %s\n\n' % (rpos + self.r_offset))

        out.write("Score: %s\n" % self.score)
        out.write("Matches: %s (%.1f%%)\n" % (self.matches, self.identity * 100))
        out.write("Mismatches: %s\n" % (self.mismatches,))
        out.write("CIGAR: %s\n" % self.cigar_str)

####################################################################################################
# 
# End
# 
####################################################################################################
