from _typeshed import Incomplete
from collections.abc import Generator
from contextlib import contextmanager
from ctypes import Structure, Union
import types

import numpy as np

INTEROPLIBDIR: Incomplete
INTEROPLIB_FILENAME: str
INTEROPLIB: str
ENVIRONPATH: str
REMOTE_MODULE_ON: bool

def initLibraryEnv(remoteArgs) -> None: ...
@contextmanager
def environ(env) -> Generator[None]: ...

class Session(Structure): ...

class LumApiSession:
    iapi: Incomplete
    handle: Incomplete
    __doc__: str
    def __init__(self, iapiArg, handleArg) -> None: ...

class LumString(Structure): ...
class LumMat(Structure): ...
class LumNameValuePair(Structure): ...
class LumStruct(Structure): ...
class LumList(Structure): ...
class ValUnion(Union): ...
class Any(Structure): ...

def lumWarning(message) -> None: ...
def initLib(remoteArgs): ...

class LumApiError(Exception):
    value: Incomplete
    def __init__(self, value) -> None: ...

def verifyConnection(handle): ...

biopen = open

def extractsHostnameAndPort(remoteArgs): ...
def open(
    product,
    key: Incomplete | None = None,
    hide: bool = False,
    serverArgs={},
    remoteArgs={},
): ...
def close(handle) -> None: ...
def evalScript(handle, code, verifyConn: bool = False) -> None: ...
def getVar(handle, varname, verifyConn: bool = False): ...
def putString(handle, varname, value, verifyConn: bool = False) -> None: ...
def putMatrix(handle, varname, value, verifyConn: bool = False) -> None: ...
def putDouble(handle, varname, value, verifyConn: bool = False) -> None: ...
def putStruct(handle, varname, values, verifyConn: bool = False) -> None: ...
def putList(handle, varname, values, verifyConn: bool = False) -> None: ...
def packMatrix(handle, value): ...
def unpackMatrix(handle, value): ...
def isIntType(value): ...

class MatrixDatasetTranslator:
    @staticmethod
    def applyConventionToStruct(d) -> None: ...
    @staticmethod
    def createStructMemberPreTranslators(d): ...

class PointDatasetTranslator:
    @staticmethod
    def applyConventionToStruct(
        d, geometryShape, paramShape, removeScalarDim
    ) -> None: ...
    @staticmethod
    def createStructMemberPreTranslators(d, numGeomDims): ...

class RectilinearDatasetTranslator:
    @staticmethod
    def applyConventionToStruct(d) -> None: ...
    @staticmethod
    def createStructMemberPreTranslators(d): ...

class UnstructuredDatasetTranslator:
    @staticmethod
    def applyConventionToStruct(d) -> None: ...
    @staticmethod
    def createStructMemberPreTranslators(d): ...

class PutTranslator:
    @staticmethod
    def translateStruct(handle, value): ...
    @staticmethod
    def translateList(handle, values): ...
    @staticmethod
    def translate(handle, value): ...
    @staticmethod
    def createStructMemberPreTranslators(value): ...
    @staticmethod
    def putStructMembers(handle, value): ...
    @staticmethod
    def putListMembers(handle, value): ...

class GetTranslator:
    @staticmethod
    def translateString(strVal): ...
    @staticmethod
    def recalculateSize(size, elements): ...
    @staticmethod
    def translate(handle, d, element): ...
    @staticmethod
    def applyLumDatasetConventions(d) -> None: ...
    @staticmethod
    def getStructMembers(handle, value): ...
    @staticmethod
    def getListMembers(handle, value): ...

def removePromptLineNo(strval): ...
def appCallWithConstructor(self, funcName, *args, **kwargs): ...
def appCall(self, name, *args): ...
def lumTypes(argList): ...

class SimObjectResults:
    def __init__(self, parent) -> None: ...
    def __dir__(self): ...
    def __getitem__(self, name): ...
    def __getattr__(self, name): ...
    def __setattr__(self, name, value): ...

class GetSetHelper(dict):
    def __init__(self, owner, name, **kwargs) -> None: ...
    def __getitem__(self, key): ...
    def __setitem__(self, key, val) -> None: ...
    def __getattr__(self, key): ...
    def __setattr__(self, key, val): ...

class SimObjectId:
    name: Incomplete
    index: Incomplete
    def __init__(self, id) -> None: ...

class SimObject:
    results: Incomplete
    def __init__(self, parent, id) -> None: ...
    def build_nested(self, properties): ...
    def __dir__(self): ...
    def __getitem__(self, key): ...
    def __setitem__(self, key, item) -> None: ...
    def __getattr__(self, name): ...
    def __setattr__(self, name, value): ...
    def getParent(self): ...
    def getChildren(self): ...

class Lumerical:
    keepCADOpened: Incomplete
    handle: Incomplete
    syncUserFunctionsFlag: bool
    userFunctions: Incomplete
    def __init__(
        self, product, filename, key, hide, serverArgs, remoteArgs, **kwargs
    ) -> None: ...
    def __extractKeepCADOpenedArgument__(self, serverArgs): ...
    def __del__(self) -> None: ...
    def __enter__(self) -> "Lumerical": ...
    def __exit__(
        self,
        type: type[BaseException] | None,
        value: BaseException | None,
        traceback: types.TracebackType | None,
    ) -> None: ...
    def __getattr__(self, name): ...
    def __open__(
        self,
        iapi,
        product,
        key: Incomplete | None = None,
        hide: bool = False,
        serverArgs={},
        remoteArgs={},
    ): ...
    def close(self) -> None: ...
    def eval(self, code) -> None: ...
    def getv(self, varname): ...
    def putv(self, varname, value) -> None: ...
    def getObjectById(self, id): ...
    def getObjectBySelection(self): ...
    def getAllSelectedObjects(self): ...

class INTERCONNECT(Lumerical):
    def __init__(
        self,
        filename: Incomplete | None = None,
        key: Incomplete | None = None,
        hide: bool = False,
        serverArgs={},
        remoteArgs={},
        **kwargs,
    ) -> None: ...

class DEVICE(Lumerical):
    def __init__(
        self,
        filename: Incomplete | None = None,
        key: Incomplete | None = None,
        hide: bool = False,
        serverArgs={},
        remoteArgs={},
        **kwargs,
    ) -> None: ...

class FDTD(Lumerical):
    def __init__(
        self,
        filename: Incomplete | None = None,
        key: Incomplete | None = None,
        hide: bool = False,
        serverArgs={},
        remoteArgs={},
        **kwargs,
    ) -> None: ...
    def __enter__(self) -> "FDTD": ...
    def getfdtdindex(
        self, materialName, f_range, f_min, f_max, verifyConn: bool = False
    ) -> np.ndarray: ...
    def stackrt(self, materialIndex, thickness, f_range, verifyConn: bool = False): ...

class MODE(Lumerical):
    def __init__(
        self,
        filename: Incomplete | None = None,
        key: Incomplete | None = None,
        hide: bool = False,
        serverArgs={},
        remoteArgs={},
        **kwargs,
    ) -> None: ...

# class Lumerical(object):
#     def __init__(self, product, filename, key, hide, serverArgs, remoteArgs, **kwargs):
#         """Keyword Arguments:
#                 script: A single string containing a script filename, or a collection of strings
#                         that are filenames. Preffered types are list and tuple, dicts are not
#                         supported. These scripts will run after the project specified by the
#                         project keyword is opened. If no project is specified, they will run
#                         in a new blank project.

#                 project: A single string containing a project filename. This project will be
#                          opened before any scripts specified by the script keyword are run.
#         """
#          # this is to keep backward compatibility with applications that need a CAD running until the 
#          # Python interpreter shuts down
#         self.keepCADOpened = self.__extractKeepCADOpenedArgument__(serverArgs)
#         iapi = initLib(remoteArgs)
#         handle = self.__open__(iapi, product, key, hide, serverArgs, remoteArgs)
#         self.handle = LumApiSession(iapi, handle)

#         self.syncUserFunctionsFlag = False
#         self.userFunctions = set() # variable to keep track of all added user methods

#         # get a list of commands from script interpreter and register them
#         # an error here is a constructor failure to populate class methods
#         try:
#             self.eval('api29538 = getcommands;')
#             commands = self.getv('api29538').split("\n")
#             commands = [x for x in commands if len(x) > 0 and x[0].isalpha()]
#             self.eval('clear(api29538);')
#         except:
#             close(self.handle)
#             raise

#         try:
#             with biopen(INTEROPLIBDIR + '/docs.json') as docFile:
#                 docs = json.load(docFile)
#         except:
#             docs = {}

#         # add methods to class corresponding to Lumerical script
#         # use lambdas to create closures on the name argument.
#         keywordsLumerical = ['for', 'if', 'else', 'exit', 'break', 'del', 'eval', 'try', 'catch', 'assert', 'end',
#                              'true', 'false', 'isnull']
#         deprecatedScriptCommands = ['addbc', 'addcontact', 'addeigenmode', 'addpropagator', 'deleteallbc', 'deletebc',
#                                     'getasapdata', 'getbc', 'getcompositionfraction', 'getcontact', 'getglobal',
#                                     'importdoping', 'lum2mat', 'monitors', 'new2d', 'new3d', 'newmode', 'removepropertydependency',
#                                     'setbc', 'setcompositionfraction', 'setcontact', 'setglobal', 'setsolver', 'setparallel',
#                                     'showdata', 'skewness', 'sources', 'structures']
#         functionsToExclude = keywordsLumerical + deprecatedScriptCommands

#         addScriptCommands =    ['add2drect', 'add2dpoly', 'addabsorbing', 'addanalysisgroup', 'addanalysisprop',
#                                 'addanalysisresult', 'addbandstructuremonitor', 'addbulkgen', 'addchargemesh',
#                                 'addchargemonitor', 'addchargesolver', 'addcircle', 'addconvectionbc',
#                                 'addctmaterialproperty', 'addcustom', 'adddeltachargesource', 'addelectricalcontact',
#                                 'addelement', 'addemabsorptionmonitor', 'addemfieldmonitor', 'addemfieldtimemonitor',
#                                 'addemmaterialproperty', 'adddevice', 'adddgtdmesh', 'adddgtdsolver', 'adddiffusion',
#                                 'addimplant', 'adddipole', 'adddope', 'addeffectiveindex', 'addefieldmonitor',
#                                 'addelectricalcontact', 'addelement', 'addeme', 'addemeindex', 'addemeport',
#                                 'addemeprofile', 'addfde', 'addfdtd', 'addfeemsolver', 'addfeemmesh', 'addgaussian',
#                                 'addgridattribute', 'addgroup', 'addheatfluxbc', 'addheatfluxmonitor', 'addheatmesh',
#                                 'addheatsolver', 'addhtmaterialproperty', 'addimport', 'addimportdope',
#                                 'addimportedsource', 'addimportgen', 'addimportheat', 'addimporttemperature',
#                                 'addindex', 'addjfluxmonitor', 'addlayer', 'addlayerbuilder',
#                                 'addimportnk', 'addmesh', 'addmode', 'addmodeexpansion', 'addmodelmaterial',
#                                 'addmodesource', 'addmovie', 'addobject', 'addparameter', 'addpath', 'addpec',
#                                 'addperiodic', 'addplane', 'addplanarsolid', 'addpmc', 'addpml', 'addpoly', 'addpower',
#                                 'addprofile', 'addproperty', 'addpyramid', 'addradiationbc', 'addrect',
#                                 'addring', 'addsimulationregion', 'addsphere', 'addstructuregroup', 'addsurface',
#                                 'addsurfacerecombinationbc', 'addtemperaturebc', 'addtemperaturemonitor', 'addtfsf',
#                                 'addthermalinsulatingbc', 'addthermalpowerbc', 'addtime', 'addtriangle',
#                                 'adduniformheat', 'adduserprop', 'addvarfdtd', 'addvoltagebc', 'addwaveguide']

#         for name in [n for n in commands if n not in functionsToExclude]:
#             if name in addScriptCommands:
#                 method = (lambda x: lambda self, *args, **kwargs:
#                 appCallWithConstructor(self, x, args, **kwargs))(name)
#             else:
#                 method = (lambda x: lambda self, *args: appCall(self, x, args))(name)
#             method.__name__ = str(name)
#             try:
#                 method.__doc__ = docs[name]['text'] + "\n" + docs[name]['link']
#             except:
#                 pass
#             setattr(Lumerical, name, method)

#         # change the working directory to match Python program
#         # load or run any file provided as argument
#         # an error here is a constructor failure due to invalid user argument
#         try:
#             if REMOTE_MODULE_ON is False: # we are not on remote mode
#                 self.cd(os.getcwd())
#             if filename is not None:
#                 if filename.endswith('.lsf'):
#                     self.feval(filename)
#                 elif filename.endswith('.lsfx'):
#                     self.eval(filename[:-5] + ';')
#                 else:
#                     self.load(filename)

#             if kwargs is not None:
#                 if 'project' in kwargs:
#                     self.load(kwargs['project'])
#                 if 'script' in kwargs:
#                     if type(kwargs['script']) is not str:
#                         for script in kwargs['script']:
#                             if script.endswith('.lsfx'):
#                                 self.eval(script[:-5] + ';')
#                             else:
#                                 self.feval(script)
#                     else:
#                         if kwargs['script'].endswith('.lsfx'):
#                             self.eval(kwargs['script'][:-5] + ';')
#                         else:
#                             self.feval(kwargs['script'])
#         except:
#             close(self.handle)
#             raise