# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com



_LAZY = {
    "p2text": "_p2text",
    "complementary_color": "_complementary_color",
    "bezier_curve_S": "_bezier_curve_S",
    "grid_average": "_grid_average",
    "text_similarity": "_text_similarity",
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
