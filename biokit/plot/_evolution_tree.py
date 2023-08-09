# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


class EvolutionTree(object):
    def __init__(self, tree):
        self.tree = tree
        self.__level_nodes_y = {}

    def get_all_keys(self, tree):
        keys = []
        for key in tree.keys():
            keys.append(key)
            if len(tree[key]) > 0:
                keys.extend(self.get_all_keys(tree[key]))
        return keys

    def get_nodes_per_level(self, tree=None, level=0, nodes_per_level=None):
        if tree is None:
            tree = self.tree
        if not nodes_per_level:
            nodes_per_level = []
        if len(nodes_per_level) <= level:
            nodes_per_level.append(0)
        nodes_per_level[level] += 1

        if isinstance(tree, dict):
            for key in tree.keys():
                self.get_nodes_per_level(tree[key], level + 1, nodes_per_level)
        elif isinstance(tree, list):
            for item in tree:
                self.get_nodes_per_level(item, level + 1, nodes_per_level)
        return nodes_per_level

    def get_tree_depth(self, tree=None):
        if tree is None:
            tree = self.tree
        if isinstance(tree, dict):
            if len(tree) == 0:
                return 1
            else:
                depths = []
                for key in tree.keys():
                    depths.append(self.get_tree_depth(tree[key]) + 1)
                return max(depths)
        elif isinstance(tree, list):
            depths = []
            for item in tree:
                depths.append(self.get_tree_depth(item))
            return max(depths)
        else:
            return 0

    def treeplot(self, tree=None, tree_name='root', ax=None, level=0, color_dict=None):
        level_nodes_y = self.__level_nodes_y
        if tree is None:
            tree = self.tree
        if not ax:
            nodes_per_level = self.get_nodes_per_level()
            fig, ax = plt.subplots(figsize=(len(nodes_per_level) * 2, max(nodes_per_level) * 2))
        if not color_dict:
            color_dict = dict(zip(self.get_all_keys(tree), sns.hls_palette(len(self.get_all_keys(tree)))))
        color_dict['root'] = 'grey'

        # 绘制进化树
        childs_y = []
        if len(tree) == 0:  # 如果自身是空字典, 根据level_nodes_y的坐标占位判断自身y坐标并返回
            self_y = level_nodes_y.get(level, 0)  # 自身占据y
            if self_y <= level_nodes_y.get(level - 1, 0):  # 如果自身y坐标小于上一级的y坐标, 则把上一级的y坐标设为自身y坐标
                self_y = level_nodes_y.get(level - 1, 0) + 1
            level_nodes_y[level] = self_y + 1  # 更新level_nodes_y, 表示本级已经有节点占据了y坐标, 并把下一个点的坐标设为y+1
            # 绘制节点
            ax.scatter(level, self_y, s=100, c=color_dict[tree_name], zorder=1,edgecolors='k',linewidths=1)
            ax.text(x=level, y=self_y + 0.3, s=tree_name, ha='center', va='center', fontsize=10)
            return self_y, ax
        else:
            # 如果自身不是叶节点, 根据子节点的y坐标计算自身y坐标
            for key in tree.keys():
                childs_y.append(
                    self.treeplot(tree[key], tree_name=key, ax=ax, level=level + 1, color_dict=color_dict)[0])
            self_y = np.mean(childs_y)  # 自身占据子节点y坐标的中间
            level_nodes_y[level] = self_y + 1  # 更新level_nodes_y, 表示本级已经有节点占据了y坐标, 并把下一个点的坐标设为y+1
            # 绘制节点
            ax.scatter(level, self_y, s=100, c=color_dict[tree_name], zorder=1,edgecolors='k',linewidths=1)
            ax.text(x=level, y=self_y + 0.3, s=tree_name, ha='center', va='center', fontsize=10)
            # 绘制从本级到下一级全部节点的连线
            for child_y in childs_y:
                ax.plot([level, level + 0.5, level + 0.5, level + 1], [self_y, self_y, child_y, child_y], c='k',
                        linewidth=1, zorder=0)
            if tree_name == 'root':
                # 如果已经递归到制根节点，对ax进行设置, 并清除level_nodes_y
                ax.axis('off')
                ax.set_ylim(-0.5, max(level_nodes_y.values()) + 0.5)
                self.__level_nodes_y = {}
            return self_y, ax

