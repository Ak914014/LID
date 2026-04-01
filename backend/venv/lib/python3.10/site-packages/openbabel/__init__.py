import sys
import warnings
import os
if os.name == 'nt':
    base_dir = os.path.abspath(os.path.dirname(__file__))
    os.add_dll_directory(os.path.join(base_dir, "bin"))


from . import openbabel

__version__ = "3.1.0"

if "ON" == "ON":
    import os
    base_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["BABEL_LIBDIR"] = os.path.join(base_dir, "lib", "openbabel", __version__)
    os.environ["BABEL_DATADIR"] = os.path.join(base_dir, "share", "openbabel", __version__)

_thismodule = sys.modules[__name__]

class OBProxy(object):
    def __getattr__(self, name):
        if hasattr(_thismodule, name):
            return getattr(_thismodule, name)
        elif hasattr(openbabel, name):
            warnings.warn('"import openbabel" is deprecated, instead use "from openbabel import openbabel"')
            return getattr(openbabel, name)
        else:
            raise AttributeError

sys.modules[__name__] = OBProxy()
