# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


def p2text(p, cutoff: dict = None):
    """
    将p值转换为显著性符号(*,**,***,ns)
    :param p:
    :param cutoff:
    :return:
    """
    if not isinstance(cutoff, dict):
        cutoff = {0.05: '*', 0.01: '**', 0.001: '***',0.0001: '****'}
    for cutoff_value in sorted(list(cutoff.keys())):
        if p <= cutoff_value:
            return cutoff[cutoff_value]
    return 'ns'
