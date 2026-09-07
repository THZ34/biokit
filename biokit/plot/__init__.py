


_LAZY = {
    "cumulative_bar": "baseplot",
    "cumulative_barh": "baseplot",
    "parallel_categories": "complexplot",
    "testbox": "_box",
    "Circos": "_circos",
    "pathway_circos": "_circos",
    "pathway_circos_style2": "_circos",
    "lr_circos": "_circos",
    "scicolors": "_colors",
    "show_scicolors": "_colors",
    "scicmap": "_colors",
    "show_scicmap": "_colors",
    "forest_plot": "_cox",
    "crosstab_plot": "_crosstab",
    "EvolutionTree": "_evolution_tree",
    "flow_chart": "_fit_models",
    "heatmap": "_heatmap",
    "kaplan_meier": "_kaplan_meier",
    "lassocv": "_lasso",
    "rocplot": "_roc_plot",
    "rocplots": "_roc_plot",
    "mutation_timescape": "_timescape",
    "timescape": "_timescape",
    "oncoplot": "_variants_landscape",
    "volcano_plot": "_volcano_plot",
    "metascape_dotplot": "_dot",
    "draw_neural_net_double": "_neuronetwork",
    "draw_neural_net_single": "_neuronetwork",
    "create_fig": "_fig",
    "radarplot": "_radar",
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
