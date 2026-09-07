from ._fit_models import *
from ._cell_annotation import *

_LAZY = {
    "cox": "_cox",
    "group4": "_group4",
    "signature_score": "_signature_score",
    "pathway_enrichment": "_pathway_enrichment",
}


def __getattr__(name):
    if name in _LAZY:
        from importlib import import_module
        mod = import_module(f".{_LAZY[name]}", __name__)
        value = getattr(mod, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["anno_cells","anno_level_cells","cox","extract_cell_types_recursive","feature_ranking","feature_selection","fit_models","get_all_subclasses","get_max_level","group4","leave_one_out","pathway_enrichment","performance_evaluation","rfe_features","signature_score","sklearn_models"]
