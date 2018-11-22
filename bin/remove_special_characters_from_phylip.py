#!/usr/bin/env python

import re
import sys

with open(sys.argv[1]) as f:
    for line in f:
        print(re.sub(r"[\[\]();,:']", "", line).rstrip())

