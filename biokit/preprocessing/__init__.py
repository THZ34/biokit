from ._gene_trans import *


_LAZY = {
    "km_best_cutoff": "_kaplan_meier",
    "read_vcf": "_read_vcf",
    "refactor_json": "_refactor_json",
    "sampleinfo_stat": "_sample_statistic",
    "read_aachange": "_variants_landscape",
    "sort_mutation": "_variants_landscape",
    "vcf_to_mutation": "_vcf_to_mutation",
}


def __getattr__(name):
    if name in _LAZY:
        from importlib import import_module
        mod = import_module(f".{_LAZY[name]}", __name__)
        value = getattr(mod, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["detect_version","download_ensembl","download_gencode","genename_version_convert","get_ensembl_genename_df","get_gencode_genename_df","get_genename_df","km_best_cutoff","read_aachange","read_vcf","refactor_json","sampleinfo_stat","sort_mutation","vcf_to_mutation"]
