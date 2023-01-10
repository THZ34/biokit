from math import pi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from biokit.data import load_hg19_karyo
from adjustText import adjust_text


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


def coordinate_correction(points):
    """极坐标角度矫正至 0 ~ 2pi

    :param points: np.array 输入点坐标
    :return: np.array 矫正后坐标
    """
    points[:, 0] = points[:, 0] % (2 * pi)
    return points


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


def coordinates_transformation(points, to='polar'):
    """坐标系转换, 极坐标系使用0~2pi角度，而非0-360

    :param points: np.array 二位点坐标
    :param to: 'angle' 或 'polar' ， 直角坐标或极坐标
    :return: np.array 转换后的坐标
    """
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
    p1, p2, p3, p4 = coordinate_correction(np.array([p1, p2, p3, p4]))
    # 点转换直角坐标
    transformation = coordinates_transformation(np.array([p1, p2, p3, p4]), to='angle')
    angle_p1, angle_p2, angle_p3, angle_p4 = transformation
    # 曲线转换极坐标
    bezier1 = coordinates_transformation(bezier_curve(angle_p1, angle_p2, n=n), to='polar')
    bezier2 = coordinates_transformation(bezier_curve(angle_p3, angle_p4, n=n), to='polar')

    # 曲线间圆弧
    interval1 = fill_endpoints(p4, p1)
    interval2 = fill_endpoints(p2, p3)
    area = np.concatenate([interval1, bezier1, interval2, bezier2], axis=0)
    return area


class Circos(object):
    def __init__(self, bottom=5, layer_spacing=0.5, background: pd.DataFrame = None):
        """初始化一个circos类"""

        self.background = None
        self.background_scale = None
        self.fig = None
        self.ax = None
        self.ymax = 0
        self.rings = {}
        self.bottom = bottom
        self.layer_spacing = layer_spacing
        self.kinds = ['bar', 'barh', 'scatter', 'line', 'bezier', 'bezierarea', 'text']
        if background:
            self.set_base(background)
        else:
            self.init_default_base_layer()

    def add_ring(self):
        """向 self.rings 中添加一个环

        """
        ring = []
        current_ring_num = max(list(self.rings.keys())) + 1
        self.rings[current_ring_num] = ring
        return ring

    def add_layer(self, data, value_col, kind, ring, color_dict=None, linewidth=None, height=None, linepoints=None,
                  fontsize=None, cmap='viridis'):
        """向环中添加一个图层

        :param fontsize: matplotlib fontsize
        :param data: pd.DataFrame
        :param value_col: any
        :param kind: str
        :param ring: list
        :param color_dict: dict
        :param linewidth: int or float
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

        layer = {'kind': kind, 'value': value_col}
        if kind == 'barh':
            layer['data'] = self.coordinating(data)
            layer['height'] = height or 0.8
            values = data[value_col].unique()
            layer['color'] = color_dict or dict(zip(values, sns.hls_palette(len(values))))

        elif kind == 'bar':
            layer['data'] = self.coordinating(data)
            pass

        elif kind == 'line':
            layer['data'] = self.coordinating(data)
            layer['linewidth'] = linewidth or 5

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
            layer['linewidth'] = linewidth or 5
            layer['linepoints'] = linepoints or 50

        elif kind == 'text':
            layer['data'] = self.coordinating(data)
            layer['fontsize'] = fontsize or 10

        ring.append(layer)
        return layer

    def set_base(self, chr_df):
        """设定基础坐标层，其他层的横坐标会基于这一层进行矫正

        :param chr_df:
        :return: None
        """
        chr_df['length'] = chr_df['end']
        interval_proportion = 0.1
        ring_length = chr_df['end'].sum() * (1 + interval_proportion) / (2 * pi)
        chr_interval = chr_df['end'].sum() * interval_proportion / chr_df.shape[0]
        chr_df['chr_ring_start'] = [(chr_df['end'][:i].sum() + chr_interval * i) for i in range(chr_df.shape[0])]
        chr_df['chr_ring_start'] = chr_df['chr_ring_start'] / ring_length
        chr_df['length'] = chr_df['length'] / ring_length
        self.background = chr_df
        self.background_scale = ring_length  # 比例尺，将染色体距离缩放到2 * pi

    def init_default_base_layer(self):
        """默认基层为hg19染色体及核型

        :return:
        """
        chr_order = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                     'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',
                     'chrY']
        karyoband_df = load_hg19_karyo()
        karyoband_df.columns = ['chr', 'start', 'end', 'name', 'type']

        # 染色体
        chr_df = []
        for chrom in chr_order:
            chr_df.append([chrom, 0, karyoband_df[karyoband_df['chr'] == chrom]['end'].max()])
        chr_df = pd.DataFrame(chr_df, columns=['chr', 'start', 'end'])
        self.set_base(chr_df)

        # 条带
        karyoband_df = self.coordinating(karyoband_df)
        karyoband_df['length'] = karyoband_df['end'] - karyoband_df['start']
        color_dict = dict(zip(['gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100', 'acen', 'gvar', 'stalk'],
                              ['white', 'lightgrey', 'grey', 'dimgrey', 'black', 'lightcoral', 'pink', 'red']))
        band_layer = {'kind': 'barh', 'value': 'type', 'color': color_dict, 'data': karyoband_df.copy()}
        self.rings[0] = [band_layer]

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

    def plot_base(self):
        ax = self.ax
        ymax = self.ymax
        background = self.background

        # 染色体
        chr_df = background
        ax.barh(y=ymax, height=0.8, left=chr_df['chr_ring_start'], width=chr_df['length'], color='white',
                edgecolor='grey')
        for chrom, start, length in chr_df[['chr', 'chr_ring_start', 'length']].to_numpy():
            ax.text(x=start + length / 2, y=ymax + 1, s=chrom, va='center', ha='center',
                    rotation=360 * (start + length / 2) / (2 * pi) + 90)

    def plot_rings(self):
        """循环绘制每个环"""
        rings = self.rings
        ymax = self.ymax
        for ring_num in rings:
            y = ymax - ring_num
            ring = rings[ring_num]
            # 逐层绘制
            for layer in ring:
                self.plot_layer(layer, y)

    def plot_layer(self, layer, y):
        """画图
        """
        data = layer['data']
        kind = layer['kind']
        value = layer['value']
        color_dict = layer.get('color', 'deepskyblue')
        height = layer.get('height', 0.8)
        linewidth = layer.get('linewidth', 5)
        linepoints = layer.get('linepoints', 50)
        fontsize = layer.get('fontsize', 50)

        if kind == 'barh':
            self.barhplot(data, value, height, y, color_dict)
        elif kind == 'bar':
            pass
        elif kind == 'bezier':
            self.bezierplot(data, value, y, linewidth, linepoints, color_dict)
        elif kind == 'bezierarea':
            self.bezierplot(data, value, y, linewidth, linepoints, color_dict)
        elif kind == 'text':
            self.textanno(data, value, y, fontsize)

    def show(self):
        rings = self.rings
        bottom = self.bottom
        ymax = bottom + len(rings)
        self.ymax = ymax

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection='polar')
        ax.set_yticks([])
        ax.set_xticks([])
        ax.grid(False)
        ax.spines['polar'].set_visible(False)
        ax.set_ylim(0, ymax + 1)
        self.ax = ax

        self.plot_base()
        self.plot_rings()

    def barhplot(self, data, value, height, y, color_dict):
        ax = self.ax
        for color in color_dict:
            temp_data = data[data[value] == color]
            ax.barh(y=y, left=temp_data['start'], width=temp_data['length'], height=height,
                    color=color_dict[color])

    def bezierplot(self, data, value, y, linewidth, linepoint, color_dict):
        """

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
        data[['p0_x', 'p0_y']] = coordinates_transformation(data[['p0_x', 'y']].to_numpy(), to='angle')
        data[['p1_x', 'p1_y']] = coordinates_transformation(data[['p1_x', 'y']].to_numpy(), to='angle')
        data[['vertice_x', 'vertice_y']] = coordinates_transformation(data[['vertice_x', 'vertice_y']].to_numpy(),
                                                                      to='angle')

        for p0_x, p0_y, p1_x, p1_y, vertice_x, vertice_y, color in data[
            ['p0_x', 'p0_y', 'p1_x', 'p1_y', 'vertice_x', 'vertice_y', value]].to_numpy():
            color = color_dict[color]
            p0 = np.array([p0_x, p0_y])
            p1 = np.array([p1_x, p1_y])
            vertice = np.array([vertice_x, vertice_y])
            curve = bezier_curve(p0, p1, vertices=vertice, n=linepoint)
            curve = coordinates_transformation(curve, to='polar')
            ax.plot(curve[:, 0], curve[:, 1], color=color, linewidth=linewidth, alpha=0.3)

    def bezierareaplot(self, data, value, y, linewidth, linepoints, color_dict):
        ax = self.ax
        for start1, end1, start2, end2, color in data[['start1', 'end1', 'start2', 'end2', value]].to_numpy():
            color = color_dict[color]
            p1 = np.array([start1, y])
            p2 = np.array([end1, y])
            p3 = np.array([start2, y])
            p4 = np.array([end2, y])
            area = bezier_area(p1, p2, p3, p4, n=linepoints)
            ax.fill(area[:, 0], area[:, 1], color=color, alpha=0.3, linewidth=linewidth)

    def textanno(self, data, value, y, fontsize):
        """标注文字"""
        ax = self.ax

        # ax.scatter(x=data[['start', 'end']].mean(1), y=[y] * data.shape[0], color='red', s=5)
        new_texts = []
        for chrom, start, end, text in data[['chr', 'start', 'end', value]].to_numpy():
            # new_texts.append(ax.text(x=(start + end) / 2, y=y, s=text, va='center', ha='center',
            #                          rotation=360 * (start + end) / (4 * pi) + 90, fontsize=fontsize))
            new_texts.append(ax.text(x=(start + end) / 2, y=y, s=text, rotation=360 * (start + end) / (4 * pi) + 90,
                                     fontsize=fontsize, va='center', ha='center'))
        adjust_text(new_texts, only_move={'text': 'xy'}, va='center', ha='center')
