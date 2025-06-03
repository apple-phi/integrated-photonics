import json

class LumEncoder(json.JSONEncoder):
    TYPE_NODE_NAME: str
    DATA_NODE_NAME: str
    SIZE_NODE_NAME: str
    COMPLEX_NODE_NAME: str
    TYPE_NODE_NAME_MATRIX: str
    TYPE_NODE_NAME_CELL: str
    def default(self, obj): ...

class LumDecoder(json.JSONDecoder):
    TYPE_NODE_NAME: str
    DATA_NODE_NAME: str
    SIZE_NODE_NAME: str
    COMPLEX_NODE_NAME: str
    TYPE_NODE_NAME_MATRIX: str
    TYPE_NODE_NAME_CELL: str
    def __init__(self, *args, **kwargs) -> None: ...
    def object_hook(self, dct): ...
