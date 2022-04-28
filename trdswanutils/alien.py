import json
import collections
from pathlib import Path
from pprint import pprint

import alienpy.alien

jalien = None

def alien_init():
    global jalien

    alienpy.alien.setup_logging()

    try:
        # info = %system alien-token-info
        # print("AliEn token info:", info)

        jalien = alienpy.alien.AliEn()
        jalien.ProcessMsg('whoami -v')

    except SystemExit:
        print("ERROR initializing AliEn session")
        print("Make sure you have a valid token by opening a SWAN terminal and running `alien-token-init`")

def alien_ls(filepattern):
    global jalien

    AlienFileInfo = collections.namedtuple('AlienFileInfo', ['path', 'size'])
    ret = jalien.run(f'ls -l {filepattern}').ansdict
    if 'results' in ret:
        return [AlienFileInfo(Path(x['path']), int(x['size'])) for x in ret['results']]
    else:
        return list()


alien_init()
