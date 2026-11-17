#!/usr/bin/env python
"""Command line of RFs_Ophiolite_Codes.

python rfo.py --help              the commands in the order of the workflow
python rfo.py <command> --help    the options of one command
"""

from rf_ophiolite.cli import main

if __name__ == '__main__':
    main()
