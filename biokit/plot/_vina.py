# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import os

# pymol作图优化模板


"""
reinitialize
cd Y:/data/project/AE0004
load vina_best_poses/955_1201531_best.pdbqt
load 5VHG_final.pdbqt
bg white
hide everything
set_color color_A, [0.4,0.76,0.65]
set_color color_B, [0.99,0.55,0.38]
set_color hb_col,  [0.20,0.40,0.85]
set_color sb_col,  [0.85,0.25,0.70]

pseudoatom label_A, 955_1201531_best
label label_A, "Dehy"
set label_color, black, label_A
set label_size, 24, label_A
set label_position, [0,0,30], label_A

pseudoatom label_B, 5VHG_final
label label_B, "Beclin1"
set label_color, black, label_B
set label_size, 24, label_B
set label_position, [0,0,30], label_B

set depth_cue, off
set fog, off
set orthoscopic, off
set field_of_view, 70

select iface_A, byres 955_1201531_best within 4 of Beclin1
select iface_B, byres Beclin1 within 4 of 955_1201531_best
show sticks, iface_A or iface_B
set stick_radius, 0.15, iface_A or iface_B
color color_A, iface_A
color color_B, iface_B
show spheres, iface_A or iface_B
set sphere_scale, 0.2, iface_A or iface_B
set sphere_transparency, 0.0, iface_A or iface_B
set stick_transparency, 0.0, iface_A or iface_B

show surface, all
color color_A, 955_1201531_best
color color_B, Beclin1
set transparency, 0.5, all

remove solvent
find_pairs iface_A, iface_B, mode=1, cut=3.5, angle=35
select hb_res, byres find_pairs and iface_A or iface_B
color hb_col, hb_res
distance hbonds, iface_A, iface_B, mode=2, cutoff=3.5, angle=35
hide labels, hbonds
set dash_width, 2.0
set dash_color, hb_col

select salt_A, (iface_A and (resn ASP+GLU)) and (name OE1+OE2+OD1+OD2) within 4 of (iface_B and (resn LYS+ARG) and name NZ+NH1+NH2)
select salt_B, (iface_B and (resn ASP+GLU)) and (name OE1+OE2+OD1+OD2) within 4 of (iface_A and (resn LYS+ARG) and name NZ+NH1+NH2)
select salt_res, byres salt_A or salt_B
color sb_col, salt_res
distance saltbridges, salt_A, salt_B, mode=2, cutoff=4.0
hide labels, saltbridges
set dash_width, 2.0
set dash_color, sb_col

center iface_A or iface_B
zoom iface_A or iface_B, 10
clip slab, 40, iface_A or iface_B
set ambient, 0.18
set direct, 0.82
set antialias, 2
set ray_shadow, off
set two_sided_lighting, on

set ray_trace_mode, 1
set ray_texture, 1
ray 2400, 2400
png Be_interface.png, dpi=300


"""



####################################
# %% 高饱和卡通风格
"""

# 重置 PyMOL 会话
reinitialize

# 设置工作目录
cd Y:/data/project/AE0004

# 加载分子
load vina_best_poses/955_1201531_best.pdbqt
load 5VHG_final.pdbqt

# 设置背景颜色为白色
bg white

# 设置 5VHG_final 的颜色为青色
set_color protein_blue, [0.212, 0.549, 0.741]   # #368CBD
color protein_blue, 5VHG_final

set_color ligand_C, [0.9, 0.5, 0.1]
color ligand_C, 955_1201531_best

# 将 5VHG_final 设置为表面风格
show surface, 5VHG_final

# 提高渲染精度
set ray_trace_mode, 3        
set ray_trace_frames, 3      
set ray_opaque_background, 0 

# 更改表面的材质为半透明玻璃风格
set transparency, 0.3, 5VHG_final  
set shininess, 100, 5VHG_final     
set reflect, 0.7, 5VHG_final      
set specular, 0.5, 5VHG_final     
set ambient, 0, 5VHG_final        
set direct, 0, 5VHG_final         

# 渲染图像
ray
png beclin1_dehy.png,dpi=1200


select ligand, 955_1201531_best
select pocket, byres (5VHG_final within 4 of ligand)

# 显示 pocket 残基
show sticks, pocket
color gray70, pocket and elem C

# 显示 ligand 为 sticks
show sticks, ligand
set stick_radius, 0.2, ligand

# 氢键
dist hbonds, (ligand), (pocket), mode=2
set dash_color, marine
set dash_width, 2

# 盐桥（带电原子之间）
select pos_res, pocket and resn ARG+LYS+HIS
select neg_res, pocket and resn ASP+GLU
dist saltbridge, (pos_res), (neg_res), cutoff=4
set dash_color, magenta
set dash_width, 2

# 标注结合残基
label pocket and name CA, "%s%s" % (resn,resi)
set label_color, black
set label_size, 20

# 聚焦到结合位点
zoom ligand, 8




"""