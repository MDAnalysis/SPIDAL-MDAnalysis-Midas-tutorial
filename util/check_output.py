#!/uar/bin/env python
# check if output files are missing

import os.path
import json

if __name__ =="__main__":
   m = json.load(open("manifest.json"))
   for line in m:
        print("{0:20s}     {1}".format(line[0], "OK" if os.path.exists(line[0]) else "MISSING"))
