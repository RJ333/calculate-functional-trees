#!/usr/bin/env python

"""
This script removes the all characters from phylip that leads to problem with RAxML.
"""

import re
import sys

with open(sys.argv[1]) as f:
    for line in f:
        print(re.sub(r"[\[\]();,:']", "", line).rstrip())

