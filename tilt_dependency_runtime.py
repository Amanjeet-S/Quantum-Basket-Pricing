"""Optional dependency discovery for the retained tilt synthesis experiment.

Normal installed packages are supported. QBP_SYNTHESIS_VENDOR can point to
an existing dependency target, or repository-local synthesis_vendor is used
if present. This helper neither installs nor copies third-party libraries.
"""
from pathlib import Path
import os
import sys


def add_optional_vendor():
    root=Path(__file__).resolve().parent
    configured=os.environ.get('QBP_SYNTHESIS_VENDOR')
    vendor=Path(configured).expanduser().resolve() if configured else root/'synthesis_vendor'
    if configured and not vendor.is_dir():
        raise RuntimeError('QBP_SYNTHESIS_VENDOR must name an existing dependency directory')
    if vendor.is_dir():
        value=str(vendor)
        if value not in sys.path:sys.path.insert(0,value)
        return vendor
    return None


def mpmath_runtime():
    add_optional_vendor()
    try:
        import mpmath as mp
    except ImportError as error:
        raise RuntimeError('Tilt generation requires mpmath1.4.1. Use an existing installed package or set QBP_SYNTHESIS_VENDOR; optional versions are pinned in requirements-tilt-synthesis.txt.') from error
    if mp.__version__!='1.4.1':
        raise RuntimeError('Retained tilt generation requires mpmath1.4.1 for reproducible source/version evidence')
    return mp
