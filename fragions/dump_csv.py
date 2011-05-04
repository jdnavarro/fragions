import sys
import argparse
import shelve
import csv

import numpy

def dump_csv(db, fh):
    fields = ('is_target', 'scan', 'raw_fn', 'mz', 'charge', 'score',
              'sum_inten', 'std', 'mean', 'peaks')
    header = dict()
    for field in fields:
        header[field] = field
    writer = csv.DictWriter(fh, fields)
    writer.writerow(header)
    for scan, psm in db.iteritems():
        inten = numpy.array([i[1] for i in psm['peaks']])
        psm['sum_inten'] = inten.sum()
        psm['std'] = inten.std()
        psm['mean'] = inten.mean()
        psm['scan'] = scan
        writer.writerow(psm)

def main():
    parser = argparse.ArgumentParser(description='Dumps the results in a '
                                                 'CSV file')
    parser.add_argument('-db', '--dbname', default='exp.db')
    parser.add_argument('outfile', nargs='?',
                         type=argparse.FileType('w'),
                         default=sys.stdout)
    args = parser.parse_args()
    db = shelve.open(args.dbname)
    dump_csv(db, args.outfile)
