"""
Compatibility functions (currently for Python 2.6).
"""


try:
    # Python 2.7 and up.
    from collections import Counter
except ImportError:
    from _counter import Counter
