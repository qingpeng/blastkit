import sys
dir = '/Users/t/dev/blastkit/lib'

if dir not in sys.path:
    sys.path.insert(0, dir)

import blastkit_config
