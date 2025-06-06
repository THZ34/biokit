from math import pi

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class Circos(object):
    def __init__(self, bottom=5, layer_spacing=0.5, background: pd.DataFrame = None, ax=None):
        """初始化一个circos类"""

        self.background = None
        self.background_scale = None
        self.background_ring = []
        self.fig = None
        self.ymax = 0
        self.rings = {}
        self.bottom = bottom
        self.layer_spacing = layer_spacing
        self.kinds = ['bar', 'barh', 'scatter', 'line', 'bezier', 'bezierarea', 'text']
        self.gene_ref = None
        self.genes = None
        if isinstance(background, pd.DataFrame):
            self.set_base(background)
        if not ax:
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(projection='polar')
            ax.set_yticks([])
            ax.set_xticks([])
            ax.grid(False)
            ax.spines['polar'].set_visible(False)
        self.ax = ax

    def add_ring(self):
        """向 self.rings 中添加一个环

        """
        ring = []
        if not self.rings:
            current_ring_num = 1
        else:
            current_ring_num = max(list(self.rings.keys())) + 1
        self.rings[current_ring_num] = ring
        return ring

    def add_space(self, n):
        if not self.rings:
            current_ring_num = 1
        else:
            current_ring_num = max(list(self.rings.keys())) + 1
        for i in range(1, n + 1):
            self.rings[current_ring_num + i] = []

    def add_layer(self, data, value_col, kind, ring, color_dict=None, height=None, linepoints=None,
                  size=None, alpha=None, marker=None, rotation=None, vmax=None, width=None, hue_color=None,
                  hue_label=None, smooth=None, n_smooth=None, **kwargs):
        """向环中添加一个图层

        :param smooth:
        :param hue_label:
        :param hue_color:
        :param width:
        :param vmax:
        :param rotation:
        :param marker:
        :param alpha:
        :param size:
        :param data: pd.DataFrame
        :param value_col: any
        :param kind: str
        :param ring: list
        :param color_dict: dict
        :param height: int or float
        :param linepoints: int
        :return: layer

        example layer:

        {'kind': 'barh',
        'value': 'type',
        'color': {'gneg': 'white', 'gpos25': 'lightgrey', 'gpos50': 'grey', 'gpos75': 'dimgrey', 'gpos100': 'black',
        'acen': 'lightcoral', 'gvar': 'pink', 'stalk': 'red'},
        'data':        chr     start       end    name   color    length
        0     chr1  0.000000  0.004244  p36.33    gneg  0.004244
        1     chr1  0.004244  0.009964  p36.32  gpos25  0.005720
        2     chr1  0.009964  0.013285  p36.31    gneg  0.003321
        3     chr1  0.013285  0.016975  p36.23  gpos25  0.003690
        4     chr1  0.016975  0.023433  p36.22    gneg  0.006458
        ..     ...       ...       ...     ...     ...       ...
        857  chr22  5.790448  5.796722   q13.1    gneg  0.006274
        858  chr22  5.796722  5.802626   q13.2  gpos50  0.005904
        859  chr22  5.802626  5.810376  q13.31    gneg  0.007750
        860  chr22  5.810376  5.812221  q13.32  gpos50  0.001845
        861  chr22  5.812221  5.815735  q13.33    gneg  0.003514

        """
        data = data.copy()
        # 检查画图类型是否正确
        if kind not in self.kinds:
            raise ValueError(f'kind allowed {self.kinds}')

        layer = {'kind': kind, 'value_col': value_col, 'alpha': alpha}

        if kind == 'barh':
            layer['data'] = self.coordinating(data)
            layer['height'] = height or 0.8
            values = data[value_col].unique()
            layer['color'] = color_dict or dict(zip(values, sns.hls_palette(len(values))))

        elif kind == 'bar':
            if not isinstance(value_col, list):
                value_col = [value_col]
            layer['data'] = self.coordinating(data)
            layer['value_cols'] = value_col
            layer['width'] = width or 0.8
            layer['vmax'] = vmax or data[value_col].max().max() * 1.25
            layer['hue_color'] = hue_color or sns.hls_palette(len(value_col))
            layer['hue_label'] = hue_label or value_col

        elif kind == 'line':
            layer['data'] = self.coordinating(data)
            layer['linewidth'] = size or 5
            layer['smooth'] = smooth or False
            layer['n_smooth'] = n_smooth or 100
            layer['color'] = color_dict or 'black'

        elif kind == 'bezier' or kind == 'bezierarea':
            values = data[value_col].unique()
            layer['color'] = color_dict or dict(zip(values, sns.hls_palette(len(values))))

            temp_df1 = data[['chr1', 'start1', 'end1']].copy()
            temp_df1.columns = ['chr', 'start', 'end']
            data[['chr1', 'start1', 'end1']] = self.coordinating(temp_df1)

            temp_df2 = data[['chr2', 'start2', 'end2']].copy()
            temp_df2.columns = ['chr', 'start', 'end']
            data[['chr2', 'start2', 'end2']] = self.coordinating(temp_df2)

            layer['data'] = data
            layer['linewidth'] = size or 5
            layer['linepoints'] = linepoints or 50
            layer['kwargs'] = kwargs

        elif kind == 'text':
            layer['data'] = self.coordinating(data)
            layer['fontsize'] = size or 10
            # layer['ha'] = ha or 'center'
            # layer['va'] = va or 'center'
            layer['rotation'] = rotation or 90

        elif kind == 'scatter':
            layer['data'] = self.coordinating(data)
            layer['scattersize'] = size or 10
            layer['marker'] = marker or 'o'
            layer['height'] = height or 0
            values = data[value_col].unique()
            layer['color'] = color_dict or dict(zip(values, sns.hls_palette(len(values))))

        ring.append(layer)
        return layer

    def set_base(self, chr_df, interval_proportion=0.1, keep_space=0.5):
        """设定基础坐标层，其他层的横坐标会基于这一层进行矫正

        :param keep_space: 留空的范围
        :param interval_proportion:
        :param chr_df:
        :return: None
        """
        chr_df['length'] = chr_df['end']
        ring_length = chr_df['end'].sum() * (1 + interval_proportion) * (1 + keep_space) / (2 * pi)
        chr_interval = chr_df['end'].sum() * interval_proportion / chr_df.shape[0]
        chr_df['chr_ring_start'] = [(chr_df['end'][:i].sum() + chr_interval * i) for i in range(chr_df.shape[0])]
        chr_df['chr_ring_start'] = chr_df['chr_ring_start'] / ring_length
        chr_df['length'] = chr_df['length'] / ring_length
        self.background = chr_df
        self.background_scale = ring_length  # 比例尺，将染色体距离缩放到2 * pi

    def init_hg19_base_layer(self):
        """默认基层为hg19染色体及核型

        :return:
        """
        from biokit.data import load_hg19_karyo
        karyoband_df = load_hg19_karyo()
        karyoband_df.columns = ['chr', 'start', 'end', 'name', 'type']
        chr_order = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
                     'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
                     'chr22', 'chrX', 'chrY']
        # 染色体
        chr_df = []
        for chrom in chr_order:
            chr_df.append([chrom, 0, karyoband_df[karyoband_df['chr'] == chrom]['end'].max()])
        chr_df = pd.DataFrame(chr_df, columns=['chr', 'start', 'end'])
        self.set_base(chr_df, 0.1, 0)

        # 条带
        karyoband_df = self.coordinating(karyoband_df)
        karyoband_df['length'] = karyoband_df['end'] - karyoband_df['start']
        color_dict = dict(zip(['gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100', 'acen', 'gvar', 'stalk'],
                              ['white', 'lightgrey', 'grey', 'dimgrey', 'black', 'lightcoral', 'pink', 'red']))
        band_layer = {'kind': 'barh', 'value_col': 'type', 'color': color_dict, 'data': karyoband_df.copy()}
        self.background_ring.append(band_layer)

    def init_gene_base_layer(self, genes):
        from biokit.data import load_hg19_ref
        ref_df = load_hg19_ref()
        print(set(genes) - set(ref_df['name2']))

        # 选择每个基因的最长转录本
        temp_refs = []
        for gene in genes:
            temp_ref = ref_df[ref_df['name2'] == gene]
            temp_ref['txLength'] = temp_ref['txEnd'] - temp_ref['txStart']
            temp_ref.sort_values(by='txLength', ascending=False, inplace=True)
            temp_ref = ref_df.loc[temp_ref.index[0]]
            temp_refs.append(temp_ref)

        ref_df = pd.DataFrame(temp_refs)
        ref_df = ref_df[['chrom', 'txStart', 'txEnd', 'name2', 'exonStart', 'exonEnd']]
        ref_df.columns = ['chr', 'start', 'end', 'gene', 'exonstart', 'exonend']

        chr_df = ref_df[['gene', 'start', 'end']].copy()
        chr_df.rename({'gene': 'chr'}, axis=1, inplace=True)
        chr_df['end'] -= chr_df['start']
        chr_df['start'] = 0
        self.set_base(chr_df, 0.1, 0)

        gene_ref = ref_df[['chr', 'gene', 'start', 'end']].copy()
        gene_ref.index = gene_ref['gene']
        self.gene_ref = gene_ref

    def gene_coordinating(self, df):
        """基于self.gene_ref中的基因坐标，将新输入的染色体坐标映射到基因上"""
        df = df.copy()
        gene_ref = self.gene_ref

        if 'gene' not in df.columns:  # 根据坐标判断基因
            genes = []
            for chrom, start, end in df[['chr', 'start', 'end']].to_numpy():
                gene = gene_ref[(gene_ref['chr'] == chrom) &
                                (gene_ref['start'] <= start) &
                                (gene_ref['end'] >= end)].index[0]
                genes.append(gene)
        df['gene'] = genes
        df['start'] -= gene_ref.loc[df['gene']]['start']
        df['end'] -= gene_ref.loc[df['gene']]['start']
        df['chr'] = df['gene']
        return df

    def coordinating(self, df):
        """基于self.background中的染色体坐标，将新输入的染色体坐标映射到环上
        """
        df = df.copy()
        chr_df = self.background
        chr_ring_start_dict = dict(zip(chr_df['chr'], chr_df['chr_ring_start']))
        background_scale = self.background_scale
        df[['start', 'end']] = df[['start', 'end']] / background_scale
        df['start'] = [start + chr_ring_start_dict[chrom] for chrom, start in df[['chr', 'start']].to_numpy()]
        df['end'] = [start + chr_ring_start_dict[chrom] for chrom, start in df[['chr', 'end']].to_numpy()]
        return df

    def plot_base(self, color_dict=None, alpha=1, y=None, fontproperties=None, plot_chrom=True, text_chrom=True,
                  chrname_extend=0, plot_band=True, fontsize_scale=False, ):
        ax = self.ax
        rings = self.rings
        bottom = self.bottom
        if not y:
            y = bottom + len(rings)
        background = self.background
        chr_df = background

        if not color_dict:
            color_dict = dict(zip(chr_df['chr'], ['grey'] * chr_df.shape[0]))

        if plot_chrom:
            ax.barh(y=y, height=0.8, left=chr_df['chr_ring_start'], width=chr_df['length'],
                    color=[color_dict[chrom] for chrom in chr_df['chr']],
                    edgecolor='grey', alpha=alpha)

        if not fontproperties:
            fontproperties = {}
        fontsize = fontproperties.get('fontsize', 20)
        fontweight = fontproperties.get('fontweight', 'bold')
        # 染色体名
        if text_chrom:
            for chrom, start, length in chr_df[['chr', 'chr_ring_start', 'length']].to_numpy():
                ax.text(x=start + length / 2, y=y + 2 + chrname_extend, s=chrom, va='center',
                        ha=fontproperties.get('ha', 'center'),
                        rotation=360 * (start + length / 2) / (2 * pi) - 90 + fontproperties.get('rotation', 0),
                        fontsize=fontsize * np.log(1 + length) / np.log(1.45) if fontsize_scale else fontsize,
                        fontweight=fontweight, rotation_mode='anchor')
        # 染色体条带
        if plot_band:
            for layer in self.background_ring:
                self.plot_layer(layer, y)

    def plot_rings(self, rings=None):
        """循环绘制每个环"""
        if not rings:
            rings = self.rings
        ymax = self.ymax
        for ring_num in rings:
            print(f'Plot ring {ring_num}')
            y = ymax - ring_num
            ring = rings[ring_num]
            # 逐层绘制
            for layer in ring:
                print(f'Plot ring {ring_num} layer {layer["kind"]}')
                self.plot_layer(layer, y)

    def plot_layer(self, layer, y):
        """画图
        """
        data = layer['data']
        kind = layer['kind']
        value_col = layer['value_col']
        color_dict = layer.get('color', 'deepskyblue')
        height = layer.get('height', 0.8)
        linewidth = layer.get('linewidth', 50)
        linepoints = layer.get('linepoints', 50)
        fontsize = layer.get('fontsize', 50)
        alpha = layer.get('alpha', 1)
        scattersize = layer.get('scattersize', 5)
        marker = layer.get('marker', 'O')
        rotation = layer.get('rotation', 0)
        smooth = layer.get('smooth', False)
        n_smooth = layer.get('n_smooth', 100)
        kwargs = layer.get('kwargs', {})

        if kind == 'barh':
            self.barhplot(data, value_col, height, y, color_dict)
        elif kind == 'bar':
            width = layer.get('width', 0.8)
            vmax = layer.get('vmax', data[value_col].max().max() * 1.25)
            hue_color = layer.get('hue_color', sns.hls_palette(len(value_col)))
            hue_label = layer.get('hue_label', value_col)
            self.barplot(data, value_col, y, width, vmax, hue_color, hue_label)
        elif kind == 'bezier':
            self.bezierplot(data, value_col, y, linewidth, linepoints, color_dict, alpha)
        elif kind == 'bezierarea':
            self.bezierareaplot(data, value_col, y, linewidth, linepoints, color_dict, alpha, **kwargs)
        elif kind == 'text':
            self.textanno(data, value_col, y, fontsize, rotation)
        elif kind == 'scatter':
            self.scatter(data, value_col, y, height, scattersize, marker, color_dict)
        elif kind == 'line':
            self.lineplot(data, value_col, y, linewidth, smooth, n_smooth, color_dict)

    def draw(self):
        rings = self.rings
        bottom = self.bottom
        ymax = bottom + len(rings)
        self.ymax = ymax
        self.ax.set_ylim(0, ymax + 1)
        self.plot_rings()

    def show(self):
        self.draw()
        plt.show()

    def barhplot(self, data, value, height, y, color_dict, **kwargs):
        ax = self.ax
        for color in color_dict:
            temp_data = data[data[value] == color]
            ax.barh(y=y, left=temp_data['start'], width=temp_data['length'], height=height,
                    color=color_dict[color])

    def bezierplot(self, data, value, y, linewidth, linepoint, color_dict, alpha, **kwargs):
        """

        :param alpha:
        :param data:
        :param value:
        :param y:
        :param linewidth:
        :param linepoint:
        :param color_dict:
        :return:
        """
        ax = self.ax

        data = data.copy()
        data['y'] = y
        data['p0_x'] = data[['start1', 'end1']].mean(1)
        data['p1_x'] = data[['start2', 'end2']].mean(1)
        data['vertice_x'] = data[['p0_x', 'p1_x']].mean(1)
        data['vertice_y'] = 0
        data[['p0_x', 'p0_y']] = Circos.coordinates_transformation(data[['p0_x', 'y']].to_numpy(), to='angle')
        data[['p1_x', 'p1_y']] = Circos.coordinates_transformation(data[['p1_x', 'y']].to_numpy(), to='angle')
        data[['vertice_x', 'vertice_y']] = Circos.coordinates_transformation(
            data[['vertice_x', 'vertice_y']].to_numpy(), to='angle')

        for p0_x, p0_y, p1_x, p1_y, vertice_x, vertice_y, color in data[
            ['p0_x', 'p0_y', 'p1_x', 'p1_y', 'vertice_x', 'vertice_y', value]].to_numpy():
            color = color_dict[color]
            p0 = np.array([p0_x, p0_y])
            p1 = np.array([p1_x, p1_y])
            vertice = np.array([vertice_x, vertice_y])
            curve = Circos.bezier_curve(p0, p1, vertices=vertice, n=linepoint)
            curve = Circos.coordinates_transformation(curve, to='polar')
            ax.plot(curve[:, 0], curve[:, 1], color=color, linewidth=linewidth, alpha=alpha)

    def bezierareaplot(self, data, value, y, linewidth, linepoints, color_dict, alpha, **kwargs):
        ax = self.ax
        for start1, end1, start2, end2, color in data[['start1', 'end1', 'start2', 'end2', value]].to_numpy():
            color = color_dict[color]
            p1 = np.array([start1, y])
            p2 = np.array([end1, y])
            p3 = np.array([start2, y])
            p4 = np.array([end2, y])
            area = Circos.bezier_area(p2, p3, p4, p1, n=linepoints)
            ax.fill(area[:, 0], area[:, 1], color=color, alpha=alpha, linewidth=linewidth, edgecolor='white', )

    def textanno(self, data, value, y, fontsize, rotation, min_interval=0.01, **kwargs):
        """标注文字"""
        ax = self.ax
        data = data.copy()
        data['center'] = data[['start', 'end']].mean(1)
        adjust_xs = []
        last_x = data['center'].max() - 2 * pi
        data.sort_values(by='center', inplace=True)
        for x in data['center']:
            if x - last_x < min_interval:
                new_x = max(x, last_x) + min_interval
            else:
                new_x = x
            adjust_xs.append(new_x)
            last_x = new_x

        data['adjust_x'] = adjust_xs
        for chrom, adjust_x, text in data[['chr', 'adjust_x', value]].to_numpy():
            ax.text(x=adjust_x, y=y + 2, s=f'{text:>{len(text) * 4}}', rotation=360 * adjust_x / (2 * pi) + rotation,
                    fontsize=fontsize, ha='center', va='center')
        self.foldline(y=y, adjust_xs=adjust_xs, true_xs=data['center'])

    def foldline(self, y, adjust_xs, true_xs, height=2, **kwargs):
        # 文字注释占2-3个环的宽度，外层环显示调整过坐标的文字以避免互相覆盖，内层环通过折线连接到文字的实际坐标
        ax = self.ax
        for adjust_x, true_x in zip(adjust_xs, true_xs):
            ax.plot([true_x, true_x, adjust_x, adjust_x], [y, y + 0.15 * height, y + 0.85 * height, y + height],
                    c='black')

    def lineplot(self, data, value_col, y, linewidth, smooth, n_smooth, color_dict, ylim=4, **kwargs):
        data['x'] = data[['start', 'end']].mean(1)
        data.sort_values(by='x', ascending=True, inplace=True)
        # 两个区间的中点可能相同导致报错
        dup_x = data['x'].value_counts()[data['x'].value_counts() > 1].index
        if dup_x.size > 0:
            for x in dup_x:
                xdup_data = data[data['x'].isin(dup_x)]
                xdup_data['position'] = xdup_data['chr'] + ':' + xdup_data['start'].astype(str) + '-' + xdup_data[
                    'end'].astype(str)
                print(', '.join(xdup_data['position'].to_list()) + ' 区间中点相同')

        data.drop_duplicates(subset=['x'], inplace=True)

        data['y'] = y + (data[value_col] - data[value_col].min()) / (data[value_col].max() - data[value_col].min())

        if smooth:
            from scipy.interpolate import make_interp_spline
            smooth_x = np.linspace(data['x'].min(), data['x'].max(), n_smooth)
            smooth_y = make_interp_spline(data['x'].to_numpy(), data['y'], k=1)(smooth_x, )
        else:
            smooth_x = data['x']
            smooth_y = data['y']

        # smooth_y 标准化到ylim 之间
        if ylim:
            smooth_y = smooth_y - y
            smooth_y = smooth_y - smooth_y.min()
            smooth_y = smooth_y / smooth_y.max()
            smooth_y = smooth_y * ylim - ylim / 2
            smooth_y = smooth_y + y

        # 首尾衔接
        smooth_x = smooth_x.tolist()
        smooth_x.append(smooth_x[0])
        smooth_y = smooth_y.tolist()
        smooth_y.append(smooth_y[0])

        ax = self.ax
        ax.plot(smooth_x, smooth_y, linewidth=linewidth, color=color_dict)

    def scatter(self, data, value_col, y, height, size, style, color_dict, ylim=3, **kwargs):
        ax = self.ax
        if height == 0:
            data['y'] = y
        elif height in data.columns:
            data['y'] = y + (data[height] - data[height].min()) / (data[height].max() - data[height].min()) * ylim

        for value in color_dict.keys():
            temp_data = data[data[value_col] == value]
            ax.scatter(x=temp_data[['start', 'end']].mean(1).to_list(), y=temp_data['y'], marker=style,
                       c=color_dict[value], s=size)

    def annotate_rings(self, xs, ys, texts, fontsize=15):
        ax = self.ax
        for x, y, text in zip(xs, ys, texts):
            ax.text(x=x, y=y, s=text, ha='center', va='center', rotation=self.circos_rotation(x),
                    fontsize=fontsize)

    def barplot(self, data, value_cols, y, width=0.8, vmax=None, hue_color=None, hue_label=None, **kwargs):
        data = data.copy()
        if not vmax:
            vmax = data[value_cols].max().max() * 1.25

        if not hue_color:
            hue_color = sns.hls_palette(n_colors=len(value_cols))

        if not hue_label:
            hue_label = value_cols

        ax = self.ax
        bar_table = pd.DataFrame()
        for gene_number in data.index:
            start, end = data.loc[gene_number, ['start', 'end']]
            block_length = end - start
            bar_width = width * block_length / len(value_cols)
            bar_x = np.linspace(start, end, len(value_cols) * 2 + 1)[1::2]
            bar_height = data.loc[gene_number][value_cols] / vmax
            temp_df = pd.DataFrame([bar_x, [bar_width] * len(value_cols), bar_height, hue_label],
                                   index=['x', 'width', 'height', 'label']).T
            bar_table = pd.concat([bar_table, temp_df], axis=0)

        for label, color in zip(hue_label, hue_color):
            bar_x = bar_table[bar_table['label'] == label]['x']
            bar_height = bar_table[bar_table['label'] == label]['height']
            bar_width = bar_table[bar_table['label'] == label]['width']
            ax.bar(x=bar_x, height=bar_height, width=bar_width, color=color, label=label, bottom=y)
        self.bar_table = bar_table

    @staticmethod
    def bezier_curve(p0, p1, vertices=np.array([0, 0]), n=20):
        """贝塞尔曲线


        :param vertices: 顶点
        :param p0: 曲线端点1
        :param p1: 曲线端点2
        :param n: 曲线段数，决定曲线平滑度
        :return: np.array 曲线坐标
        """
        line_coordinates = [p0]
        for t in range(1, n):
            t0 = t * (vertices - p0) / n + p0
            t1 = t * (p1 - vertices) / n + vertices
            q0 = t * (t1 - t0) / n + t0
            line_coordinates.append(q0)
        line_coordinates.append(p1)
        return np.array(line_coordinates)

    @staticmethod
    def circos_rotation(x):
        return (x - pi / 2) / (2 * pi) * 360

    @staticmethod
    def coordinate_correction(points):
        """极坐标角度矫正至 0 ~ 2pi

        :param points: np.array 输入点坐标
        :return: np.array 矫正后坐标
        """
        points[:, 0] = points[:, 0] % (2 * pi)
        return points

    @staticmethod
    def fill_endpoints(point1, point2):
        angle1, radius1 = point1
        angle2, radius2 = point2
        if angle2 > angle1:
            interval = np.linspace(point1, point2, num=100)
        else:
            interval_1 = np.linspace(point1, np.array([2 * pi, point1[1]]), num=100)
            interval_2 = np.linspace(np.array([0, point2[1]]), point2, num=100)
            interval = np.concatenate([interval_1, interval_2], axis=0)
        return interval

    @staticmethod
    def coordinates_transformation(points, to='polar'):
        """坐标系转换, 极坐标系使用0~2pi角度，而非0-360

        :param points: np.array 二位点坐标
        :param to: 'angle' 或 'polar' ， 直角坐标或极坐标
        :return: np.array 转换后的坐标
        """
        points = np.array(points)
        new_coordinates = []
        if to == 'polar':
            for point in points:
                x, y = point
                radius = np.sqrt(x ** 2 + y ** 2)
                angle = np.arctan(y / x)
                new_coordinates.append([angle, radius])
            new_coordinates = np.array(new_coordinates)
            new_coordinates[points[:, 0] < 0, 0] += pi
        elif to == 'angle':
            for point in points:
                angle, radius = point
                x = np.cos(angle) * radius
                y = np.sin(angle) * radius
                new_coordinates.append([x, y])
            new_coordinates = np.array(new_coordinates)
        else:
            raise ValueError("to='angle' or to='polar'")

        return new_coordinates

    @staticmethod
    def bezier_area(p1, p2, p3, p4, n=100):
        """两条贝塞尔曲线连接包围的区域
        输入点坐标为极坐标
        连接顺序 p1 - p2 - p3 - p4 - p1

        :param n: 曲线分段数，影响曲线平滑度
        :param p1: p1,p2连接一条曲线
        :param p2:
        :param p3: p3,p4连接一条曲线
        :param p4:
        :return: np.array 区域边框
        """
        # 矫正极坐标 angle 至 0 ~ 2pi
        p1, p2, p3, p4 = Circos.coordinate_correction(np.array([p1, p2, p3, p4]))
        # 点转换直角坐标
        transformation = Circos.coordinates_transformation(np.array([p1, p2, p3, p4]), to='angle')
        angle_p1, angle_p2, angle_p3, angle_p4 = transformation
        # 曲线转换极坐标
        bezier1 = Circos.coordinates_transformation(Circos.bezier_curve(angle_p1, angle_p2, n=n), to='polar')
        bezier2 = Circos.coordinates_transformation(Circos.bezier_curve(angle_p3, angle_p4, n=n), to='polar')

        # 曲线间圆弧
        interval1 = Circos.fill_endpoints(p4, p1)
        interval2 = Circos.fill_endpoints(p2, p3)
        area = np.concatenate([interval1, bezier1, interval2, bezier2], axis=0)
        return area


# %%
def pathway_circos(pathway_gene_dict, color_dict=None, pathway_gene_ratio=(1, 1), connect_from='pathway', cmap=None,
                   gene_order=None, pathway_order=None, alpha=1):
    """
    Pathway Circos
    :param cmap:
    :param pathway_gene_dict:
    :param color_dict:
    :param pathway_gene_ratio:
    :param connect_from:
    :return:
    """
    pathway_names = list(pathway_gene_dict.keys())
    all_genes = sum(list(pathway_gene_dict.values()), [])

    pathway_size = pd.Series({k: len(v) for k, v in pathway_gene_dict.items()})
    pathway_size.sort_values(inplace=True, ascending=False)
    pathway_df = pd.DataFrame()
    pathway_df['end'] = pathway_size
    pathway_df['chr'] = pathway_size.index
    pathway_df['start'] = 0

    gene_counts = pd.value_counts(all_genes)
    gene_df = pd.DataFrame()
    gene_df['end'] = gene_counts
    gene_df['chr'] = gene_counts.index
    gene_df['start'] = 0
    gene_df[['start', 'end']] = gene_df[['start', 'end']] * pathway_gene_ratio[1] / pathway_gene_ratio[0]

    if not cmap:
        cmap = 'rainbow'
    if not color_dict:
        color_dict = dict(zip(pathway_names, sns.color_palette(cmap, len(pathway_names)))).update(
            dict(zip(all_genes, sns.color_palette(cmap, len(all_genes)))))

    if gene_order:
        gene_df = gene_df.loc[gene_order]
    if pathway_order:
        pathway_df = pathway_df.loc[pathway_order]

    base_df = pd.concat([pathway_df, gene_df], axis=0)
    bezier_df = []
    gene_start_dict = dict.fromkeys(gene_counts.index, 0)
    for pathway, genes in pathway_gene_dict.items():
        pathway_start = 0
        for gene in gene_counts.index:
            gene_start = gene_start_dict[gene]
            if gene in genes:
                bezier_df.append([pathway, pathway_start, pathway_start + 1, gene, gene_start, gene_start + 1])
                pathway_start += 1
                gene_start_dict[gene] += 1
    bezier_df = pd.DataFrame(bezier_df,
                             columns=['pathway', 'pathway_start', 'pathway_end', 'gene', 'gene_start', 'gene_end'])
    if connect_from == 'pathway':
        bezier_df.rename(dict(zip(['pathway', 'pathway_start', 'pathway_end', 'gene', 'gene_start', 'gene_end'],
                                  ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])), axis=1, inplace=True)

    elif connect_from == 'gene':
        bezier_df.rename(dict(zip(['pathway', 'pathway_start', 'pathway_end', 'gene', 'gene_start', 'gene_end'],
                                  ['chr2', 'start2', 'end2', 'chr1', 'start1', 'end1'])), axis=1, inplace=True)

    bezier_df['value'] = bezier_df['chr1']

    # 画图
    circos = Circos(bottom=20)
    circos.set_base(base_df, keep_space=0)
    ring = circos.add_ring()
    circos.add_layer(data=bezier_df, kind='bezierarea', value_col='value', ring=ring, alpha=alpha, linepoints=200,
                     color_dict=color_dict, size=0.5)
    circos.plot_base(fontproperties={'fontsize': 10, 'ha': 'left', 'rotation': 180}, text_chrom=True,
                     color_dict=color_dict)
    circos.ax.set_theta_offset(np.pi / 2)
    circos.draw()

    ax = circos.ax
    texts = ax.texts
    # pathway统一转移到图左侧,保持高度不变,向右对齐,通过annotation线连接到原本的坐标点
    for text in texts:
        if text.get_text() in pathway_df.index:
            theta, radius = text.get_position()
            x, y = circos.coordinates_transformation([[theta - np.pi / 2, radius]], to='angle')[0]
            new_x = 26
            new_theta, new_radius = circos.coordinates_transformation([[new_x, y]], to='polar')[0]
            new_theta = new_theta + np.pi / 2
            # text.set_position([new_radius, new_theta])
            # text.set_rotation(0)
            # text.set_ha('right')
            print(text.get_text())
            print(x, y)
            ax.annotate(text.get_text(), xy=(theta, 22), xytext=(new_theta, new_radius),
                        # ax.annotate(text.get_text(), xy=(new_x, y), xytext=(new_theta, new_radius),
                        arrowprops=dict(arrowstyle='-', color='grey'), fontsize=12, fontweight='bold', ha='right',
                        va='center')
            text.set_text('')
    ax.set_ylim(0, 27)
    return circos


# %% 通路circos style2
def pathway_circos_style2(pathway_gene_dict, center_name='Pathway', pathway_radius=0.7, pathway_node_base_size=0.03,
                pathway_node_adapt_size=0.1, pathway_line_base_width=1, pathway_line_adapt_width=3, gene_radius=1.2,
                cmap='Set2', ax=None, dpi=100):
    """
    环形分叉树
    :param pathway_gene_dict:
    :param center_name:
    :param pathway_radius:
    :param pathway_node_base_size:
    :param pathway_node_adapt_size:
    :param pathway_line_base_width:
    :param pathway_line_adapt_width:
    :param gene_radius:
    :param cmap:
    :param ax:
    :param dpi:
    :return:
    """
    import matplotlib.patches as mpatches
    import math
    from biokit.data_convert import bezier_curve_S

    if not ax:
        fig, ax = plt.subplots(figsize=(12, 12), dpi=dpi)

    # 
    pathway_gene_df = pd.DataFrame([(p, g) for p, gl in pathway_gene_dict.items() for g in gl],
                                   columns=['pathway', 'gene'])

    # ---------- 角度分配 ----------
    total_genes = len(pathway_gene_df)
    angle_step = 360 / total_genes
    pathway_gene_df['angle'] = [i * angle_step for i in range(total_genes)]

    # ---------- 计算颜色 ----------
    pathways = pathway_gene_df['pathway'].unique()
    color_dict = dict(zip(pathways, sns.color_palette(cmap, len(pathways))))

    # ---------- 计算中心-通路连线粗细,以及通路node大小 ----------
    pathway_gene_counts = pathway_gene_df.groupby('pathway').size()
    pathway_gene_counts = pathway_gene_counts / pathway_gene_counts.max()

    # ---------- 画图 ----------
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.axis('off')

    # 中心点
    ax.add_patch(mpatches.Circle((0, 0), 0.20, color="#f7cccc", alpha=0.6))
    ax.text(0, 0, center_name, ha='center', va='center',
            fontsize=20, fontweight='bold', color="#D62728")

    # 通路节点和基因节点
    control_point_radius = (pathway_radius + gene_radius) / 2
    for pathway in pathways:
        pathway_df = pathway_gene_df[pathway_gene_df['pathway'] == pathway]

        # 通路节点
        pathway_angle = math.radians(pathway_df['angle'].mean())
        pathway_node_x, pathway_node_y = pathway_radius * math.cos(pathway_angle), pathway_radius * math.sin(
            pathway_angle)
        pathway_node_radius = pathway_node_base_size + pathway_node_adapt_size * pathway_gene_counts[pathway]
        pathway_line_width = pathway_line_base_width + pathway_line_adapt_width * pathway_gene_counts[pathway]
        ax.add_patch(
            mpatches.Circle((pathway_node_x, pathway_node_y), radius=pathway_node_radius,
                            color=color_dict[pathway], alpha=0.6))
        ax.text(pathway_node_x, pathway_node_y, pathway, ha='center', va='center',
                fontsize=16, fontweight='bold', color=color_dict[pathway])
        ax.plot([0, pathway_node_x], [0, pathway_node_y], color=color_dict[pathway],
                linewidth=pathway_line_width, alpha=0.8)

        # 贝塞尔曲线控制点1
        control_point1 = (control_point_radius * math.cos(pathway_angle),
                          control_point_radius * math.sin(pathway_angle))

        # 基因节点 & 曲线
        for _, row in pathway_df.iterrows():
            gene_angle_rad = math.radians(row['angle'])
            # 基因node位置
            gene_x, gene_y = gene_radius * math.cos(gene_angle_rad), gene_radius * math.sin(gene_angle_rad)
            # 贝塞尔曲线控制点2
            control_point2 = (control_point_radius * math.cos(gene_angle_rad),
                              control_point_radius * math.sin(gene_angle_rad))
            # 曲线坐标
            ctrl_pts = np.array([[pathway_node_x, pathway_node_y], control_point1, control_point2, [gene_x, gene_y]])
            line = bezier_curve_S(ctrl_pts)
            ax.plot(line[:, 0], line[:, 1], color=color_dict[pathway], linewidth=1.5, alpha=0.8)

            # 基因文本
            gene_angle = row['angle']
            ratation = gene_angle + 180 if 90 <= gene_angle <= 270 else gene_angle
            text_ha = 'right' if 90 <= gene_angle <= 270 else 'left'
            ax.text(gene_x, gene_y, row['gene'], ha=text_ha, va='center',
                    fontsize=10, color=color_dict[pathway],
                    rotation=ratation, rotation_mode='anchor')

    ax.set_title("Pathway–Gene Network (Curved Line Links)", fontsize=22, fontweight='bold')
    return ax
