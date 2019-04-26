'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
import sys

from sbc_ngs import ice, pathway


def main(args):
    '''main method.'''
    if args[0] == 'pathway':
        pathway.main(args[1:])
    elif args[0] == 'ice':
        ice.main(args[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
