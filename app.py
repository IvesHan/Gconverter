import streamlit as st
import pandas as pd
import mygene
import io
import math
from gprofiler import GProfiler
import plotly.express as px
import plotly.graph_objects as go

# --- 页面基础设置 ---
st.set_page_config(page_title="BioInfo Tool Pro v3.0", layout="wide", page_icon="🧬")

st.title("🧬 基因组学多功能工具 (v3.0 - 可视化大师版)")
st.markdown("""
**核心功能：** ID 转换 | 注释 | 富集分析 | **交互式可视化 (Plotly)**
""")

# --- 全局物种映射 ---
species_map = {
    "Human (Homo sapiens)": (9606, 'hsapiens'),
    "Mouse (Mus musculus)": (10090, 'mmusculus'),
    "Rat (Rattus norvegicus)": (10116, 'rnorvegicus')
}

# --- 侧边栏：全局物种 ---
st.sidebar.header("🛠️ 全局设置")
selected_species_key = st.sidebar.selectbox("选择物种:", options=list(species_map.keys()))
species_id, gprofiler_organism_code = species_map[selected_species_key]

# --- 辅助函数 ---
def clean_cell_data(cell):
    if isinstance(cell, list):
        cleaned_list = [str(item) if not isinstance(item, dict) else f"{item.get('chr','N/A')}:{item.get('start','N/A')}" for item in cell]
        return "; ".join(cleaned_list)
    return str(cell) if isinstance(cell, dict) else cell

# --- 主体 Tabs ---
tab1, tab2 = st.tabs(["1. 基因 ID 转换与注释", "2. 富集分析与可视化"])

# =================================================================================
# Tab 1: 基因 ID 转换 (保持原样，功能已很完善)
# =================================================================================
with tab1:
    st.header("功能一：ID 转换与详细注释")
    
    # ... (此处代码逻辑与 v2.1 相同，为节省篇幅，只保留核心逻辑结构) ...
    # 为了完整性，这里保留最简化的输入输出逻辑
    
    field_mapping = {
        '基因全名 (Name)': 'name', '别名 (Alias)': 'alias', 
        '功能简介 (Summary)': 'summary', '基因类型 (Type)': 'type_of_gene', 
        '染色体位置 (Pos)': 'genomic_pos'
    }
    
    col_t1_1, col_t1_2 = st.columns([1, 1])
    with col_t1_2:
        add_info = st.multiselect("添加额外注释:", options=list(field_mapping.keys()), default=['基因全名 (Name)'])
    
    query_fields = ['symbol', 'entrezgene', 'ensembl.gene'] + [field_mapping[i] for i in add_info]
    
    raw_text = st.text_area("输入基因 ID (每行一个):", height=100, key="t1_text")
    if st.button("🚀 开始转换", key="t1_btn"):
        if raw_text:
            gene_list = [x.strip() for x in raw_text.split('\n') if x.strip()]
            mg = mygene.MyGeneInfo()
            res = mg.querymany(gene_list, scopes='symbol,entrezgene,ensembl.gene,alias', fields=query_fields, species=species_id, as_dataframe=True)
            df_res = res.reset_index()
            for col in df_res.columns: df_res[col] = df_res[col].apply(clean_cell_data)
            st.dataframe(df_res)
            
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df_res.to_excel(writer, index=False)
            st.download_button("📥 下载 Excel", output.getvalue(), "gene_conversion.xlsx")

# =================================================================================
# Tab 2: 富集分析与可视化 (重点升级)
# =================================================================================
with tab2:
    st.header("功能二：富集分析 & 交互式绘图")
    
    # 1. 输入区
    col_in1, col_in2 = st.columns([1, 2])
    with col_in1:
        raw_text_enrich = st.text_area("粘贴差异基因列表:", height=150, placeholder="TP53\nEGFR...", key="t2_text")
        
    with col_in2:
        st.markdown("**参数配置**")
        enrich_sources = st.multiselect("富集来源:", ['KEGG', 'GO:BP', 'GO:CC', 'GO:MF', 'Reactome'], default=['KEGG', 'GO:BP'])
        p_threshold = st.slider("P-value 阈值:", 0.01, 0.1, 0.05)
        run_enrich = st.button("📈 运行分析", type="primary")

    # 2. 分析与绘图逻辑
    if run_enrich and raw_text_enrich:
        gene_list = [x.strip() for x in raw_text_enrich.split('\n') if x.strip()]
        
        try:
            gp = GProfiler(return_style='pandas')
            results = gp.profile(organism=gprofiler_organism_code, query=gene_list, sources=enrich_sources, user_threshold=p_threshold, no_iea=True)
            
            if results.empty:
                st.warning("未发现显著富集结果。")
            else:
                # --- 数据预处理 ---
                # 计算 -log10(pvalue) 用于绘图
                results['neg_log10_p'] = results['p_value'].apply(lambda x: -math.log10(x))
                # 截取过长的通路名称
                results['short_name'] = results['name'].apply(lambda x: x[:50] + '...' if len(x)>50 else x)
                
                # 保留 Top N 用于绘图 (默认 20)
                results = results.sort_values('p_value', ascending=True) # 越小越显著
                
                st.success(f"发现 {len(results)} 条显著通路！")
                
                # --- 核心升级：可视化控制面板 ---
                st.divider()
                st.subheader("📊 可视化工作室 (Visualization Studio)")
                
                # 布局：左边是控制面板，右边是图
                viz_col1, viz_col2 = st.columns([1, 3])
                
                with viz_col1:
                    with st.container(border=True):
                        st.markdown("### 🎨 绘图设置")
                        
                        # 图表类型
                        plot_type = st.selectbox("图表类型:", ["气泡图 (Dot Plot)", "柱状图 (Bar Chart)"])
                        
                        # 数据筛选
                        top_n = st.slider("展示 Top N 通路:", 5, 50, 20)
                        
                        # 颜色与主题
                        color_scale = st.selectbox("配色方案:", ["Tealgrn", "Viridis", "Plasma", "Bluered", "Portland"])
                        
                        # 字体与尺寸
                        base_font_size = st.slider("字体大小:", 8, 20, 12)
                        plot_height = st.slider("图表高度:", 400, 1000, 600)
                        
                        # 排序逻辑
                        sort_by = st.selectbox("排序依据:", ["显著性 (P-value)", "基因数量 (Count)"])
                
                # 准备绘图数据 (取 Top N)
                plot_data = results.head(top_n).copy()
                # 为了让画图时显著的在上面，通常需要反转顺序
                plot_data = plot_data.sort_values('p_value', ascending=False)
                
                with viz_col2:
                    fig = None
                    
                    # --- 气泡图逻辑 (Dot Plot) ---
                    if plot_type == "气泡图 (Dot Plot)":
                        fig = px.scatter(
                            plot_data, 
                            x="intersection_size",  # X轴：富集基因数
                            y="short_name",         # Y轴：通路名称
                            size="intersection_size", # 气泡大小
                            color="neg_log10_p",    # 颜色：显著性
                            hover_data=["p_value", "source"],
                            color_continuous_scale=color_scale,
                            labels={"intersection_size": "Gene Count", "short_name": "Pathway", "neg_log10_p": "-log10(P-value)"},
                            title=f"Top {top_n} Enriched Pathways"
                        )
                    
                    # --- 柱状图逻辑 (Bar Chart) ---
                    elif plot_type == "柱状图 (Bar Chart)":
                        fig = px.bar(
                            plot_data,
                            x="neg_log10_p",        # X轴：显著性
                            y="short_name",         # Y轴：通路名称
                            color="intersection_size", # 颜色：基因数
                            orientation='h',        # 水平柱状图
                            color_continuous_scale=color_scale,
                            labels={"neg_log10_p": "-log10(P-value)", "short_name": "Pathway", "intersection_size": "Count"},
                            title=f"Top {top_n} Enriched Pathways"
                        )
                    
                    # --- 统一美化设置 ---
                    fig.update_layout(
                        height=plot_height,
                        font=dict(size=base_font_size),
                        xaxis=dict(showgrid=True, gridcolor='lightgrey'),
                        yaxis=dict(showgrid=True, gridcolor='lightgrey'),
                        plot_bgcolor='white'
                    )
                    
                    # 展示图表
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # --- 完美导出 ---
                    st.markdown("### 📥 导出结果")
                    col_exp1, col_exp2 = st.columns(2)
                    
                    # 1. 导出 Excel
                    output_excel = io.BytesIO()
                    with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
                        results.drop(columns=['neg_log10_p', 'short_name']).to_excel(writer, index=False)
                    col_exp1.download_button("下载数据表 (Excel)", output_excel.getvalue(), "enrichment_results.xlsx")
                    
                    # 2. 导出 HTML (保留交互)
                    buffer_html = io.StringIO()
                    fig.write_html(buffer_html)
                    col_exp2.download_button("下载交互式图表 (HTML)", buffer_html.getvalue().encode(), "enrichment_plot.html", help="下载后用浏览器打开，可缩放、保存为SVG/PNG")

        except Exception as e:
            st.error(f"分析出错: {e}")