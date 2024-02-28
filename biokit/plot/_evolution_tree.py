# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import matplotlib.collections as mpcollections
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


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
            ax.scatter(level, self_y, s=100, c=color_dict[tree_name], zorder=1, edgecolors='k', linewidths=1)
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
            ax.scatter(level, self_y, s=100, c=color_dict[tree_name], zorder=1, edgecolors='k', linewidths=1)
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


# %%

def draw(tree, ax=None, linecolor='black', point_color_dict=None):
    """Plot the given tree using matplotlib (or pylab).

    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.

    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).

    Example using the pyplot options 'axhspan' and 'axvline'::

        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})

    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).

    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_confidence : bool
            Whether to display confidence values, if present on the tree.
        ax : matplotlib/pylab ax
            If a valid matplotlib.ax.ax instance, the phylogram is plotted
            in that ax. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_confidence`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.

    """

    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths

    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = {tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))}

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights

    def draw_clade_lines(use_linecollection=False, orientation="horizontal", y_here=0, x_start=0, x_here=0, y_bot=0,
                         y_top=0, color="black", lw=".1", ):
        """Create a line with or without a line collection object.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if not use_linecollection and orientation == "horizontal":
            ax.hlines(y_here, x_start, x_here, color=color, lw=lw, zorder=0)
        elif use_linecollection and orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection([[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw, zorder=0))
        elif not use_linecollection and orientation == "vertical":
            ax.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection([[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw, zorder=0)
            )

    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(use_linecollection=True, orientation="horizontal", y_here=y_here, x_start=x_start,
                         x_here=x_here, color=linecolor, lw=lw)
        # Add node/taxon labels
        label = [i for i in [clade.name, clade.confidence] if i]
        label = None if len(label) == 0 else str(label[0])

        print(label)
        if label:
            ax.text(x_here, y_here, " %s" % label, ha='center', va="center", color='black', zorder=2,
                    fontsize=12)
            ax.scatter(x_here, y_here, s=400, color=point_color_dict[label], zorder=1)

        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(use_linecollection=True, orientation="vertical", x_here=x_here, y_bot=y_bot, y_top=y_top,
                             color=linecolor, lw=lw, )
        # Draw descendents
        for child in clade:
            draw_clade(child, x_here, color, lw)

    if not ax:
        ax = plt.subplots()
    if not point_color_dict:
        all_nodes = tree.get_terminals() + tree.get_nonterminals()
        labels = []

    horizontal_linecollections = []
    vertical_linecollections = []
    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    draw_clade(tree.root, 0, 'k', plt.rcParams["lines.linewidth"])

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        ax.add_collection(i)
    for i in vertical_linecollections:
        ax.add_collection(i)
    # Aesthetics
    try:
        name = tree.name
    except AttributeError:
        pass
    else:
        if name:
            ax.set_title(name)
    ax.set_xlabel("branch length")
    ax.set_ylabel("taxa")
    # Add margins around the tree to prevent overlapping the ax
    xmax = max(x_posns.values())
    ax.set_xlim(-0.05 * xmax, 1.25 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    ax.set_ylim(max(y_posns.values()) + 0.8, 0.2)
