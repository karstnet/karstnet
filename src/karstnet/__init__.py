# Note: importing the package with "import karstnet as kn"
# allows to access:

# - kn.__version__ : version (string) defined in karstnet._version.py
from karstnet._version import __version__

# - kn.<...> : all stuff defined in karstnet/base.py (e.g. the main class: kn.KGraph)
from karstnet.base import *

# - kn.utils.<...> : all stuff defined in karstnet/utils.py
from karstnet import tools
from karstnet import utils

# - kn.io.<...> : all stuff defined in karstnet/_io/*.py
from karstnet import clean
from karstnet import geom
from karstnet import io_func
from karstnet import misc
from karstnet import view

