# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
# %%
import matplotlib.pyplot as plt

rgb_color_style = {'style1': [(0, 0, 0), (19, 33, 60), (252, 163, 17), (229, 229, 229)],
                   'style2': [(244, 241, 222), (223, 122, 94), (60, 64, 91), (130, 178, 154), (242, 204, 142)],
                   'style3': [(38, 70, 83), (42, 157, 142), (233, 196, 107), (243, 162, 97), (230, 111, 81)],
                   'style4': [(246, 111, 105), (254, 179, 174), (255, 244, 242), (21, 151, 165), (14, 96, 107),
                              (255, 194, 75)],
                   'style5': [(144, 201, 231), (33, 158, 188), (19, 103, 131), (2, 48, 74), (254, 183, 5),
                              (255, 158, 2), (250, 134, 0)],
                   'style6': [(115, 186, 214), (13, 76, 109), (3, 50, 80), (2, 38, 62), (239, 65, 67), (191, 30, 46),
                              (196, 50, 63)],
                   'style7': [(231, 56, 71), (240, 250, 239), (168, 218, 219), (69, 123, 157), (29, 53, 87)],
                   'style8': [(183, 181, 160), (68, 117, 122), (69, 42, 61), (212, 76, 60), (221, 108, 76),
                              (229, 133, 93), (238, 213, 183)],
                   'style9': [(78, 171, 144), (142, 182, 156), (237, 221, 195), (238, 191, 109), (217, 79, 51),
                              (131, 64, 38)],
                   'style10': [(219, 49, 36), (252, 140, 90), (255, 223, 146), (230, 241, 243), (144, 190, 224),
                               (75, 16, 178)],
                   'style11': [(203, 153, 126), (221, 190, 169), (253, 232, 213), (184, 183, 163), (165, 165, 141),
                               (107, 112, 92)],
                   'style12': [(251, 240, 195), (84, 104, 111), (229, 123, 127), (158, 49, 80), (135, 187, 164)],
                   'style13': [(144, 201, 230), (33, 158, 188), (2, 48, 71), (255, 183, 3), (251, 132, 2)],
                   'style14': [(1, 7, 19), (1, 36, 76), (1, 53, 101), (255, 195, 0), (251, 132, 2)],
                   'style15': [(120, 0, 1), (193, 18, 33), (254, 240, 213), (0, 47, 73), (102, 155, 187)],
                   'style16': [(38, 70, 83), (40, 114, 113), (42, 157, 140), (138, 176, 125), (233, 196, 107),
                               (243, 162, 97), (230, 111, 81)],
                   'style17': [(68, 4, 90), (65, 62, 133), (48, 104, 141), (31, 146, 139), (53, 183, 119),
                               (145, 213, 66), (248, 230, 32)],
                   'style18': [(0, 0, 0), (58, 0, 100), (122, 27, 109), (189, 55, 82), (237, 104, 37), (251, 180, 26)],
                   'style19': [(0, 0, 0), (65, 14, 115), (140, 42, 129), (223, 74, 104), (252, 154, 107),
                               (252, 248, 187)],
                   'style20': [(22, 6, 138), (99, 0, 169), (158, 24, 157), (204, 73, 117), (236, 120, 83),
                               (253, 179, 46)],
                   'style21': [(164, 5, 69), (244, 111, 68), (253, 217, 133), (233, 254, 161), (127, 203, 164),
                               (75, 101, 175)]}

hex_color_dict = {}
for key in rgb_color_style:
    hex_color_dict[key] = []
    for rgb in rgb_color_style[key]:
        hex_color_dict[key].append('#%02x%02x%02x' % rgb)


def scicolors(n_colors, style='style1'):
    if n_colors <= len(hex_color_dict[style]):
        fold = n_colors // len(hex_color_dict[style])
        colors = hex_color_dict[style] * (fold + 1)
    elif n_colors == len(hex_color_dict[style]):
        colors = hex_color_dict[style]
    else:
        step = len(hex_color_dict[style]) // n_colors
        colors = hex_color_dict[style][:step:]
    return colors[:n_colors]


def show_scicolors():
    max_n_colors = max([len(i) for i in hex_color_dict.values()])
    fig, ax = plt.subplots(figsize=(max_n_colors, len(hex_color_dict)))
    for i, key in enumerate(hex_color_dict):
        for j, color in enumerate(hex_color_dict[key]):
            ax.barh(y=i, width=0.8, height=0.8, left=j, color=color)
    ax.set_xticks([])
    ax.set_yticks(range(len(hex_color_dict)))
    ax.set_yticklabels(hex_color_dict.keys(), fontsize=16)
    ax.spines.clear()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()
