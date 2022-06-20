import matplotlib.pyplot as plt
import numpy as np
from math import pi
import seaborn as sns


def bezier_curve(p0, p1, vertices, n=20):
    """贝塞尔曲线

    :param vertices: 顶点
    :param p0: 曲线端点1
    :param p1: 曲线端点2
    :param n: 曲线段数，决定曲线平滑度
    :return:
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
    """坐标系转换

    :param points: np.array 二位点坐标
    :param to: 'angle' 或 'polar' ， 直角坐标或极坐标
    :return: np.array 转换后的坐标
    """
    new_coordinates = []
    if to == 'polar':
        for point in points:
            x, y = point
            radius = np.sqrt(x ** 2 + y ** 2)
            angle = np.arctan(x / y)
            new_coordinates.append([angle, radius])
        new_coordinates = np.array(new_coordinates)
        new_coordinates[points[:, 1] < 0, 0] += pi
    elif to == 'angle':
        for point in points:
            angle, radius = point
            x = np.sin(angle) * radius
            y = np.cos(angle) * radius
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
    def __init__(self, bottom=4, layer_spacing=0.5):
        self.layers = {}
        self.bottom = bottom
        self.layer_spacing = layer_spacing

        self.fig = plt.figure(figsize=(10, 10))
        self.ax = self.fig.add_subplot(projection='polar')

    # def add_layer(self, df, kind, height):
    #     """向 self.layers 中添加一层
    #
    #     :param df: pd.DataFrame
    #     :param variables: 一个或多个变量
    #     :param kind: box,boxen,violin,bar
    #     :return: None
    #     """
    #     current_layers_num = len(self.layers)
    #     self.layers[current_layers_num + 1] = [df, kind, height]

    def plot_df(self, df, layers):
        """绘制整个数据框

        :param df: pd.DataFrame
        :param layers: Dict
        :return: None
        """
        colors = dict(zip(range(len(layers)), sns.hls_palette(len(layers))))
        x = 2 * pi * np.arange(df.shape[0]) / df.shape[0]
        for layer_num in range(len(layers)):
            var, kind, height = layers[layer_num]
            color = colors[layer_num]
            y = height * df[var] / df[var].max()
            self.plot(x, y, kind, color, layer_height=height)
            self.bottom += (height + self.layer_spacing)

    def plot(self, x, y, kind, color, layer_height=None):
        plot_functions = {'box': self.box,
                          'violin': self.violin,
                          'boxen': self.boxen,
                          'bar': self.bar,
                          'heatmap': self.heatmap}
        plot_func = plot_functions[kind]
        if kind == 'bar':
            plot_func(x=x, height=y, color=color)
        elif kind == 'heatmap':
            plot_func(x=x, y=y, height=layer_height)
        else:
            pass

    def bar(self, x, height, color='blue', width=0.5):
        ax = self.ax
        ax.bar(x, height, bottom=self.bottom, color=color, width=width / (2 * pi))

    def heatmap(self, x, y, height, width=0.95, cmap=sns.color_palette("flare", as_cmap=True)):
        ax = self.ax
        if not height:
            height = 1
        vmin = min(y)
        vmax = max(y)
        colors = [int(256 * i / (vmax - vmin)) for i in y]
        colors = cmap(colors)
        ax.bar(x=x, height=[height] * len(x), bottom=self.bottom, width=width / (2 * pi), color=colors)

    def show(self):
        ax = self.ax
        ax.set_yticks([])
        ax.grid(False)
        ax.set_xticks([])
        ax.set_ylim(0, self.bottom + 1)

    def box(self):
        pass

    def violin(self):
        pass

    def boxen(self):
        pass

    def barh(self, x, width, height, color):
        ax = self.ax
        ax.barh(y=self.bottom, left=x, width=width, height=height, color=color)

    def plot_bezier_area(self, p1, p2, p3, p4, color='blue', alpha=0.3):
        ax = self.ax
        area = bezier_area(p1, p2, p3, p4)
        ax.fill(area[:, 0], area[:, 1], color=color, alpha=alpha)
