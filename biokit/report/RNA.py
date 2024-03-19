from jinja2 import Template

# 差异基因筛选
diff_gene_filtering = """
<h2>2.4.1. 差异基因筛选</h2>
<p>在这里插入差异基因筛选的内容。</p>
"""

# 火山图
volcano_plot = """
<h2>2.4.2. 火山图</h2>
<p>在这里插入火山图的内容。</p>
"""

# 差异基因表达水平聚类分析
gene_expression_clustering = """
<h2>2.4.3. 差异基因表达水平聚类分析</h2>
<p>在这里插入差异基因表达水平聚类分析的内容。</p>
"""

# 差异基因表达水平雷达图
gene_expression_radar_plot = """
<h2>2.4.4. 差异基因表达水平雷达图</h2>
<p>在这里插入差异基因表达水平雷达图的内容。</p>
"""

# 差异基因GO富集分析
go_enrichment_analysis = """
<h2>2.4.5. 差异基因GO富集分析</h2>
<p>在这里插入差异基因GO富集分析的内容。</p>
"""

# 差异基因KEGG富集分析
kegg_enrichment_analysis = """
<h2>2.4.6. 差异基因KEGG富集分析</h2>
<p>在这里插入差异基因KEGG富集分析的内容。</p>
"""

# 差异转录因子家族分布
tf_family_distribution = """
<h2>2.6.1. 差异转录因子家族分布</h2>
<p>在这里插入差异转录因子家族分布的内容。</p>
"""

# 差异转录因子（家族）靶基因统计
tf_target_gene_statistics = """
<h2>2.6.2. 差异转录因子（家族）靶基因统计</h2>
<p>在这里插入差异转录因子（家族）靶基因统计的内容。</p>
"""

# 构建完整的HTML报告
template_html = """
<html>
<head>
<title>RNA-seq表达差异分析报告</title>
</head>
<body>
<h1>Analysis Report</h1>

{diff_gene_filtering}

{volcano_plot}

{gene_expression_clustering}

{gene_expression_radar_plot}

{go_enrichment_analysis}

{kegg_enrichment_analysis}

{tf_family_distribution}

{tf_target_gene_statistics}

</body>
</html>
"""

# 使用模板引擎渲染HTML报告
template = Template(template_html)


final_html = template.render(
diff_gene_filtering=diff_gene_filtering,
volcano_plot=volcano_plot,
gene_expression_clustering=gene_expression_clustering,
gene_expression_radar_plot=gene_expression_radar_plot,
go_enrichment_analysis=go_enrichment_analysis,
kegg_enrichment_analysis=kegg_enrichment_analysis,
tf_family_distribution=tf_family_distribution,
tf_target_gene_statistics=tf_target_gene_statistics)
