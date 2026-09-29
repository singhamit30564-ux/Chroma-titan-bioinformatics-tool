"""Chroma Titan -- a Galaxy-inspired bioinformatics tool suite.

The package mirrors the tool panel of https://usegalaxy.org/ : every tool is a
small, testable Python function that is *declared* with its inputs, outputs and
a runnable example.  The declarations are picked up by :mod:`chroma_titan.registry`
and rendered by the Streamlit front end (``app.py``) or the CLI
(``python -m chroma_titan``).

Layout::

    chroma_titan/
      registry.py     tool metadata, search, dispatch
      panel.py        usegalaxy.org panel sections/groups
      datasets.py     catalogue, loaders and previews for the bundled example data
      core/           reusable engines (io, sequences, intervals, stats, ...)
      tools/          the tools themselves, one module per panel area
"""

from chroma_titan.registry import Input, Tool, all_tools, get_tool, run_tool, search

__version__ = "6.0.0"
__all__ = [
    "Input",
    "Tool",
    "all_tools",
    "get_tool",
    "run_tool",
    "search",
    "__version__",
]
