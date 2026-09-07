# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


_LAZY = {
    "contour_map": "_contour_map",
    "parallel_categories": "_parallel_categories",
}


def __getattr__(name):
    if name in _LAZY:
        from importlib import import_module
        mod = import_module(f".{_LAZY[name]}", __name__)
        value = getattr(mod, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = sorted(_LAZY)
