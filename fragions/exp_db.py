import argparse
import os
import fnmatch
import re
import shelve

SEPARATOR = '--gc0p4Jq0M2Yt08jU534c0p'

SECTION = 'Content-Type: application/x-Mascot; name="{0}"'

PEP_REGEX = re.compile(r'^q(\d+)'      # nativeid
                        '_p(\d+)'      # rank
                        '=\d+?'        # missed cleavages
                        ',([\.\d]+?)'  # peptide mass
                        ',([-\.\d]+?)' # delta mass
                        ',\d+?'        # n ions matches
                        ',(\w+?)'      # peptide string
                        ',\d+?'        # peaks used for Ions1
                        ',(\d+?)'      # modstring
                        ',([\.\d]+?)'  # score
                        ',\d+?'        # ion series found
                        ',\d+?'        # peaks used for Ions2
                        ',\d+?'        # peaks used for Ions3
                        ';(.+)$'       # protein accessions string
                        )
SCAN_REGEX = re.compile(r'FinneganScanNumber%3a%20(\d+)')
RAW_FN = re.compile(r'(RawFile|period)%3a%20(.+raw)')

class DatParser(object):
    def __init__(self, dat_fh, db,
                 sections=('summary', 'peptides', 'query1')):
        self.dat_fh = dat_fh
        self.db = db
        self.offsets = self._get_section_offsets(sections)

    def _get_section_offsets(self, sections):
        offsets = dict()
        # Iterating over file object can't get offsets right

        section_lines = [SECTION.format(s) for s in sections]
        while True:
            line = self.dat_fh.readline()
            if not line:
                return offsets
            elif line.startswith("Content-Type"):
                line = line.strip()
                for sect_line in section_lines:
                    if sect_line == line:
                        index = section_lines.index(sect_line)
                        section = sections[index]
                        offsets[section] = self.dat_fh.tell()
                        break

    def parse_queries(self):
        self.query_scan = [None]
        self.dat_fh.seek(self.offsets['query1'])
        for line in self.dat_fh:
            if line.startswith("title="):
                line = line.strip()
                scan = re.search(SCAN_REGEX, line).group(1)
                raw_fn = re.search(RAW_FN, line).group(1)
                self.query_scan.append(scan)
            elif line.startswith("Ions1="):
                peaks = list()
                mascot_ions = line.strip().split("=")[1].split(",")
                for pair in mascot_ions:
                    mz, inten = pair.split(':')
                    peak = float(mz), float(inten)
                    peaks.append(peak)
                db_id = scan
                # shelve without writeback can't mutate values directly
                if self.db.has_key(db_id):
                    temp = self.db[db_id]
                    temp.update({'peaks': peaks, 'raw_fn': raw_fn})
                    self.db[db_id] = temp
                else:
                    self.db[db_id] = {'peaks': peaks}
            elif line.strip() == SECTION.format("index"):
                return

    def parse_summary(self):
        self.dat_fh.seek(self.offsets['summary'])
        query = 1
        for line in self.dat_fh:
            if line.startswith('qexp'):
                scan = self.query_scan[query]
                mz, charge = line.strip().split('=')[1].split(',')
                temp = self.db[scan]
                temp.update({'mz': mz, 'charge': charge})
                self.db[scan] = temp
                query += 1
            elif line.strip() == SEPARATOR:
                return

    def parse_peptides(self):
        self.dat_fh.seek(self.offsets['peptides'])
        for line in self.dat_fh:
            match = re.match(PEP_REGEX, line)
            if match:
                rank = int(match.group(2))
                if rank == 1:
                    scan = self.query_scan[int(match.group(1))]
                    target = is_target(match.group(8))
                    score = match.group(7)
                    temp = self.db[scan]
                    temp.update({'is_target': target, 'score': score})
                    self.db[scan] = temp

    def __call__(self):
        self.parse_queries()
        self.parse_summary()
        self.parse_peptides()

def is_target(accs):
    return any(acc.startswith('"IPI:IPI')
                for acc in accs.split(','))

def parse_dats(dats_path):
    for root, dirs, files in os.walk(dats_path):
        for fn in files:
            if fnmatch.fnmatch(fn, u"*.dat"):
                yield os.path.join(root, fn)
def main():
    parser = argparse.ArgumentParser(description='Generates a shelve db '
            'with all spectra from dat files.')
    parser.add_argument('dat')
    parser.add_argument('-db', '--dbname', default='exp.db')

    args = parser.parse_args()
    if args.dbname in os.listdir('.'):
        raise ValueError('{0} exists!'.format(args.dbname))
    db = shelve.open(args.dbname)
    with open(args.dat) as dat_fh:
        DatParser(dat_fh, db)()
