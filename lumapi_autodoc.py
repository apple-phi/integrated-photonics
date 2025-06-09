#!/usr/bin/env python3
"""
Auto-generate type stubs for lumapi from documentation.

TODO for future development:
- Inspect the link for each method to get more detailed parameter information
- Parse actual parameter types from documentation instead of relying on name heuristics
- Extract more precise return types from examples and descriptions
- Improve parameter name extraction to handle complex syntax (quotes, nested structures)
- Add proper union types for parameters that accept multiple types
- Consider using the Lumerical-Python conversion rules more systematically:
  * Strings remain strings
  * Real numbers become float
  * Complex numbers become 1x1 numpy arrays
  * Numpy arrays become matrices
  * Lists become cell arrays
  * Dicts become structures
- Add validation of parameter names against the actual API
- Consider generating separate overloads for different Lumerical products (FDTD, MODE, etc.)
- Parse parameter descriptions from documentation tables for better docstrings
"""
import json
import os
import re
import keyword
from pathlib import Path
from typing import Any, List
import sys

root = Path("C:\\Program Files\\Lumerical")

# Recursive search for a python/lumapi.py file in the specified directory
for dirpath, dirnames, filenames in root.walk():
    if "lumapi.py" in filenames:
        sys.path.append(str(dirpath))
        print("Found lumapi at: '", dirpath / "lumapi.py'")
        break
else:
    raise ImportError("lumapi.py not found in the specified directory tree.")

import lumapi  # type: ignore[import-untyped]
import numpy as np
import matplotlib.pyplot as plt

docs = json.load(open(Path(lumapi.INTEROPLIBDIR) / "docs.json"))  # type: ignore[no-untyped-call]

# Regexes for pulling out the ASCII "Syntax" tables + return hints
table_block_re = re.compile(r"^\+[-+]+\+\s*\n(.*?)^\+[-+]+\+", re.MULTILINE | re.DOTALL)
sig_re = re.compile(r"o\.(?P<name>\w+)\s*\(\s*(?P<params>.*?)\s*\)", re.DOTALL)
returns_re = re.compile(r"Returns the ([^.]+)\.", re.IGNORECASE)


def make_safe(name: str) -> str:
    """If name is a Python keyword, append an underscore."""
    return name + "_" if keyword.iskeyword(name) else name


def fix_escaped_underscores(text: str) -> str:
    """
    Replace all '\\_' with '_' in `text`, but for any ASCII tables
    (lines using '+' and '|' as delimiters), re-pad each column so
    that the table stays aligned.

    Returns the transformed text.
    """
    lines = text.splitlines()
    out_lines: List[str] = []
    i = 0
    sep_re = re.compile(r"^\+(-+\+)+\s*$")  # e.g. "+----+-------+"

    while i < len(lines):
        line = lines[i]
        # If this line is a separator, we assume a table starts here
        if sep_re.match(line):
            # find end of table
            j = i
            while j < len(lines) and (
                lines[j].startswith("+") or lines[j].startswith("|")
            ):
                j += 1
            table_block = lines[i:j]
            out_lines.extend(_rebuild_table(table_block))
            i = j
        else:
            # normal text: just do a global replace of "\\_" → _
            out_lines.append(line.replace(r"\_", "_"))
            i += 1

    return "\n".join(out_lines)


def _rebuild_table(block: List[str]) -> List[str]:
    """
    Given the lines of one ASCII table, first replace "\\_" → _,
    then recompute and re-pad all interior rows to the column widths
    defined by the first separator line.
    """
    # 1) replace all "\\_" → _
    block = [ln.replace(r"\_", "_") for ln in block]

    # 2) parse first separator line to get column widths
    sep = block[0]
    # find all runs of '-' between '+'; their lengths are the column widths
    widths = [len(run) for run in re.findall(r"\+(-+)", sep)]

    new_block: List[str] = []
    for ln in block:
        if ln.startswith("+"):
            # keep separator lines exactly as is
            new_block.append(ln)
        elif ln.startswith("|"):
            # split into cells; this yields ["", cell1, cell2, …, ""]
            parts = ln.split("|")
            # inner cells are parts[1:-1]
            cells = parts[1:-1]
            # rebuild each cell, preserving leading spaces but left-aligning content
            new_cells = []
            for cell, w in zip(cells, widths):
                # count leading spaces
                lead = len(cell) - len(cell.lstrip(" "))
                text = cell.strip()
                # pad: leading spaces + content + remaining spaces
                pad = w - lead - len(text)
                new_cells.append(" " * lead + text + " " * pad)
            # reassemble
            new_block.append("|" + "|".join(new_cells) + "|")
        else:
            # shouldn't happen in a well-formed table, but pass through
            new_block.append(ln)
    return new_block


def infer_param_type(param_name: str) -> str:
    """
    Infer a better type hint for a parameter based on its name.
    Based on Lumerical-Python conversion rules.
    """
    param_lower = param_name.lower()

    # String parameters (often quoted in docs)
    if any(
        word in param_lower
        for word in [
            "name",
            "str",
            "filename",
            "file",
            "path",
            "mode",
            "type",
            "label",
            "title",
        ]
    ):
        return "str"

    # Matrix/array parameters
    if any(
        word in param_lower
        for word in ["matrix", "data", "array", "field", "tet", "vtx", "mesh", "tri"]
    ):
        return "np.ndarray"

    # Numeric parameters
    if any(
        word in param_lower
        for word in [
            "min",
            "max",
            "value",
            "num",
            "count",
            "size",
            "length",
            "width",
            "height",
            "thickness",
            "option",
            "index",
        ]
    ):
        return "float"

    # Boolean-like parameters
    if any(word in param_lower for word in ["enable", "disable", "flag", "bool"]):
        return "bool"

    # Coordinate parameters
    if param_lower in ["x", "y", "z", "dx", "dy", "dz"]:
        return "float"

    if param in ["x_span", "y_span", "z_span"]:
        return "float"

    if param in ["T"]:
        return "str"

    # Common single-letter variables that are usually numeric
    if len(param_name) == 1 and param_name.lower() in "abcdefghijklmnopqrstuvwxyz":
        return "float"

    # Default to Any for unclear cases
    return "Any"


pre = r"""
# This file is auto-generated from the lumapi documentation.
from _typeshed import Incomplete
from collections.abc import Generator
from contextlib import contextmanager
from ctypes import Structure, Union
from typing import Any, Optional, List, Tuple, Dict, overload
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
# class Any(Structure): ...

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
    # def eval(self, code) -> None: ...
    def getv(self, varname): ...
    def putv(self, varname, value) -> None: ...
    def getObjectById(self, id): ...
    def getObjectBySelection(self): ...
    def getAllSelectedObjects(self): ...
"""

post = r"""
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
"""


class TableParser:
    """
    A class to parse text-based tables with variable row heights.
    This version uses a split-based strategy for robustness.
    """

    def __init__(self, table_string: str):
        """
        Initializes the parser and processes the table string.

        Args:
            table_string: A string containing the entire table to be parsed.
        """
        self.raw_table = table_string
        self.header = []
        self.data = []
        self._parse()

    def _parse(self):
        """The main method to perform the parsing logic."""
        lines = self.raw_table.strip().splitlines()
        if not lines:
            return

        # Determine the number of columns from the header row (the first data line)
        # We split by '|' and subtract 2 for the empty strings outside the table
        try:
            num_columns = len(lines[1].split("|")) - 2
        except IndexError:
            return  # Table is too short

        if num_columns <= 0:
            return

        parsed_rows = []
        current_row_parts = [""] * num_columns

        # Start from line 1 (the first row with |), skipping the top border
        for line in lines[1:]:
            # Check if the line is a separator, marking the end of a row
            if line.strip().startswith("+"):
                # Process the completed row
                # Use ' '.join(cell.split()) to normalize internal whitespace
                cleaned_row = [" ".join(cell.split()) for cell in current_row_parts]
                if any(cleaned_row):
                    parsed_rows.append(cleaned_row)

                # Reset for the next row
                current_row_parts = [""] * num_columns

            # Check if the line contains data
            elif "|" in line:
                parts = line.split("|")
                # The actual content is in parts[1], parts[2], etc.
                if len(parts) == num_columns + 2:
                    for i in range(num_columns):
                        # The cell content is the (i+1)th element from the split
                        cell_content = parts[i + 1].strip()
                        if cell_content:
                            # Append content with a space if cell already has text
                            if current_row_parts[i]:
                                current_row_parts[i] += " " + cell_content
                            else:
                                current_row_parts[i] = cell_content

        # Assign header and data from the parsed rows
        if parsed_rows:
            self.header = parsed_rows[0]
            data_rows = parsed_rows[1:]
            self.data = [dict(zip(self.header, row)) for row in data_rows]

    def to_dicts(self) -> list[dict]:
        """Returns the parsed table data as a list of dictionaries."""
        return self.data


stub_lines = str(pre).splitlines()

for cmd, info in docs.items():
    text = info["text"]
    link = info["link"]

    signatures = []  # List of tuples: (param_list, original_signature)
    return_type = "Any"
    all_raw_signatures = []

    # Find tables by looking for table boundary patterns
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        # Look for table start (line with +---+ pattern)
        if line.startswith("+") and "-" in line and line.endswith("+"):
            # Found table start, find table end
            table_start = i
            table_end = i

            # Find the last line of this table
            j = i + 1
            while j < len(lines):
                next_line = lines[j].strip()
                if (
                    next_line.startswith("+")
                    and "-" in next_line
                    and next_line.endswith("+")
                ) or next_line.startswith("|"):
                    table_end = j
                elif (
                    next_line.startswith("+")
                    and "-" in next_line
                    and next_line.endswith("+")
                ):
                    # This is the final border of the table
                    table_end = j
                    break
                elif not next_line.strip():
                    # Empty line might end the table, but check next line
                    if j + 1 < len(lines) and not (
                        lines[j + 1].strip().startswith(("+", "|"))
                    ):
                        break
                else:
                    # Non-table line
                    break
                j += 1

            # Extract the complete table text
            table_text = "\n".join(lines[table_start : table_end + 1])

            # Parse with TableParser - this is our single source of truth
            try:
                parser = TableParser(table_text)
                parsed_data = parser.to_dicts()

                if parsed_data:
                    # Look specifically in "Syntax" column for function signatures
                    for row_dict in parsed_data:
                        syntax_content = row_dict.get("Syntax", "")
                        if syntax_content:
                            # Look for function signatures like o.cmd(...)
                            for m in sig_re.finditer(syntax_content):
                                name, params = m.group("name"), m.group("params")
                                if name != cmd:
                                    continue

                                # collapse whitespace, strip trailing commas
                                params = re.sub(r"\s+", " ", params).strip().rstrip(",")
                                all_raw_signatures.append(params)
            except Exception as e:
                print(f"Warning: TableParser failed for {cmd}: {e}")
                # Don't fall back to regex - trust TableParser or skip
                pass

            i = table_end + 1
        else:
            i += 1  # Now process all found signatures and extract clean parameter lists
    for params in all_raw_signatures:
        if params:
            # Extract actual parameter names using a much simpler approach
            param_list = []
            # Split by comma and clean each parameter
            raw_params = [p.strip() for p in params.split(",")]

            for param in raw_params:
                # Remove quotes if present
                param = re.sub(r'^["\']|["\']$', "", param)
                # Clean up escaped underscores
                param = param.replace(r"\_", "_")

                # Only accept simple identifiers (no complex expressions)
                if re.match(r"^[a-zA-Z_][a-zA-Z0-9_]*$", param) and len(param) <= 15:
                    clean_param = make_safe(param)
                    param_list.append(clean_param)
                # Skip parameters that look like complex expressions, descriptions, etc.

            # Only add signature if we extracted meaningful parameters
            if param_list:
                signatures.append((param_list, f"o.{cmd}({params})"))
        else:
            # No parameters - only add if we don't already have a no-param signature
            if not any(sig[0] == [] for sig in signatures):
                signatures.append(([], f"o.{cmd}())"))

    # Remove duplicates by comparing parameter lists, but also check for parameter name conflicts
    seen_signatures = set()
    unique_signatures = []
    for param_list, orig_sig in signatures:
        # Check for duplicate parameter names within this signature
        if len(param_list) != len(set(param_list)):
            # Skip signatures with duplicate parameter names
            continue

        # Create a signature key based on parameter count and names
        sig_key = (len(param_list), tuple(param_list))
        if sig_key not in seen_signatures:
            seen_signatures.add(sig_key)
            unique_signatures.append((param_list, orig_sig))

    # If no signatures were extracted, create a generic fallback
    if not unique_signatures:
        unique_signatures.append(([], f"o.{cmd}()"))

    # 3) Heuristic on return type
    rmatch = returns_re.search(text)
    if rmatch:
        desc = rmatch.group(1).lower()
        lowered = text.lower()
        if "does not return any" in lowered:
            return_type = "None"
        elif any(k in desc for k in ("matrix", "value", "number", "array")):
            return_type = "float"
        elif any(k in desc for k in ("string", "text")):
            return_type = "str"

    safe_cmd = make_safe(cmd)

    # 4) Build the signature lines
    if len(unique_signatures) > 1:
        # Multiple overloads needed
        for i, (param_list, orig_sig) in enumerate(unique_signatures):
            stub_lines.append("    @overload")

            if param_list:
                safe_params = [make_safe(p) for p in param_list]
                # Use better type hints based on parameter names
                params_annot = ", ".join(
                    f"{s}: {infer_param_type(s)}" for s in safe_params
                )
                sig = f"def {safe_cmd}(self, {params_annot}, **kwargs: Any) -> {return_type}: ..."
            else:
                sig = f"def {safe_cmd}(self, **kwargs: Any) -> {return_type}: ..."

            stub_lines.append(f"    {sig}")
        # Add the actual implementation signature
        stub_lines.append(
            f"    def {safe_cmd}(self, *args: Any, **kwargs: Any) -> {return_type}:"
        )
    elif len(unique_signatures) == 1:
        # Single signature
        param_list, orig_sig = unique_signatures[0]
        if param_list:
            safe_params = [make_safe(p) for p in param_list]
            # Use better type hints based on parameter names
            params_annot = ", ".join(f"{s}: {infer_param_type(s)}" for s in safe_params)
            sig = (
                f"def {safe_cmd}(self, {params_annot}, **kwargs: Any) -> {return_type}:"
            )
        else:
            sig = f"def {safe_cmd}(self, **kwargs: Any) -> {return_type}:"
        stub_lines.append(f"    {sig}")
    else:
        # No signatures found in tables, use generic
        sig = f"def {safe_cmd}(self, *args: Any, **kwargs: Any) -> {return_type}:"
        stub_lines.append(f"    {sig}")

    # Docstring: break out the ASCII prose, but reformat “See Also” as NumPy-doc
    lines = text.strip().splitlines()
    # find See Also start
    see_idx = None
    for i, L in enumerate(lines):
        if L.strip().startswith("See Also"):
            see_idx = i
            break

    # split into intro vs see_also block
    intro = lines[:see_idx] if see_idx is not None else lines
    see_block = lines[see_idx + 1 :] if see_idx is not None else []
    see_block = [L.replace("()", "") for L in see_block if L.strip()]

    stub_lines.append('        """')

    # dump the intro lines
    intro = fix_escaped_underscores("\n".join(intro)).splitlines()
    for line in intro:
        stub_lines.append(f"        {line}")

    # If there was a See Also, reformat it
    if see_idx is not None:
        # extract all function names (split on commas, strip)
        names = []
        for L in see_block:
            for part in L.split(","):
                name = part.strip()
                if name:
                    names.append(name)
        # NumPy-doc See Also
        stub_lines.append("        ")
        stub_lines.append("        See Also")
        stub_lines.append("        --------\n")
        # render each as a Sphinx method reference
        meths = ", ".join(f":meth:`Lumerical.{n}`" for n in names)
        stub_lines.append(f"        {meths}")

    stub_lines.append(
        f"""\n        Link
        ----
        {link}
"""
    )
    if unique_signatures:
        # Show all the original signatures we found
        sig_examples = [orig_sig for _, orig_sig in unique_signatures]
        stub_lines.append(
            f"""
        Note
        ----
            Signature autogen'd from: {', '.join(f'`{fix_escaped_underscores(sig)}`' for sig in sig_examples)}"""
        )
    else:
        stub_lines.append(
            f"""
        Note
        ----
            Signature autogen'd from: `o.{cmd}(…)`"""
        )
    stub_lines.append('        """')
    stub_lines.append("")

# Add the postamble
stub_lines.extend(str(post).splitlines())

os.makedirs("types", exist_ok=True)
(Path("types") / "lumapi.pyi").write_text("\n".join(stub_lines), encoding="utf-8")
print("Stub generated to ./types/lumapi.pyi")
