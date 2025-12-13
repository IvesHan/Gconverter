import streamlit as st
import pandas as pd
import mygene
import io
import math
from gprofiler import GProfiler
import plotly.express as px
import plotly.graph_objects as go

# --- 页面基础设置 ---
st.set_page_config(page_title="BioInfo Tool Pro v3.1", layout="wide", page_icon="🧬")

st.title("🧬 基因组学多功能工具 (v3.1 - 修复版)")
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
    """清理 Excel 单元格数据，将列表转为字符串"""
    if isinstance(cell, list):
        cleaned_list = [str(item) if not isinstance(item, dict) else f"{item.get('chr','N/A')}:{item.get('start','N/A')}" for item in cell]
        return "; ".join(cleaned_list)
    return str(cell) if isinstance(cell, dict) else cell

# --- 主体 Tabs ---
tab1, tab2 = st.tabs(["1. 基因 ID 转换与注释", "2. 富集分析与可视化"])

# =================================================================================
# Tab 1: 基因 ID 转换与注释
# =================================================================================
with tab1:
    st.header("功能一：ID 转换与详细注释")
    
    field_mapping = {
        '基因全名 (Name)': 'name', '别名 (Alias)': 'alias', 
        '功能简介 (Summary)': 'summary', '基因类型 (Type)': 'type_of_gene', 
        '染色体位置 (Pos)': 'genomic_pos'
    }
    
    col_t1_1, col_t1_2 = st.columns([1, 1])
    with col_t1_2:
        add_info = st.multiselect("添加额外注释:", options=list(field_mapping.keys()), default=['基因全名 (Name)'])
    
    query_fields = ['symbol', 'entrezgene', 'ensembl.gene'] + [field_mapping[i] for i in add_info]
    
    # 输入区域
    input_method = st.radio("输入方式:", ("直接粘贴文本", "上传 Excel/CSV 文件"), key="t1_method")
    
    gene_list = []
    df_input = None
    col_name = "Input_ID"

    if input_method == "直接粘贴文本":
        raw_text = st.text_area("输入基因 ID (每行一个):", height=100, key="t1_text")
        if raw_text:
            gene_list = [x.strip() for x in raw_text.split('\n') if x.strip()]
            df_input = pd.DataFrame({col_name: gene_list})
    else:
        uploaded_file = st.file_uploader("上传文件", type=['xlsx', 'csv'], key="t1_file")
        if uploaded_file:
            if uploaded_file.name.endswith('.csv'): df_input = pd.read_csv(uploaded_file)
            else: df_input = pd.read_excel(uploaded_file)
            col_name = st.selectbox("选择 ID 列:", df_input.columns)
            gene_list = df_input[col_name].dropna().astype(str).tolist()

    if st.button("🚀 开始转换", key="t1_btn"):
        if not gene_list:
            st.warning("请输入基因 ID")
        else:
            with st.spinner("查询中..."):
                try:
                    mg = mygene.MyGeneInfo()
                    res = mg.querymany(gene_list, scopes='symbol,entrezgene,ensembl.gene,alias', fields=query_fields, species=species_id, as_dataframe=True)
                    df_res = res.reset_index()
                    
                    # 清洗数据
                    for col in df_res.columns: df_res[col] = df_res[col].apply(clean_cell_data)
                    
                    # 合并结果
                    if input_method == "直接粘贴文本":
                        final_df = df_res
                    else:
                        final_df = pd.merge(df_input, df_res, left_on=col_name, right_on='query', how='left')

                    st.dataframe(final_df)
                    
                    output = io.BytesIO()
                    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                        final_df.to_excel(writer, index=False)
                    st.download_button("📥 下载 Excel", output.getvalue(), "gene_conversion.xlsx")
                except Exception as e:
                    st.error(f"Error: {e}")

# =================================================================================
# Tab 2: 富集分析与可视化 (已修复 BUG)
# =================================================================================
with tab2:
    st.header("功能二：富集分析 & 交互式绘图")
    
    # 1. 输入区
    col_in1, col_in2 = st.columns([1, 2])
    with col_in1:
        raw_text_enrich = st.text_area("粘贴差异基因列表 (Symbol/Ensembl):", height=150, placeholder="TP53\nEGFR...", key="t2_text")
        
    with col_in2:
        st.markdown("**参数配置**")
        enrich_sources = st.multiselect("富集来源:", ['KEGG', 'GO:BP', 'GO:CC', 'GO:MF', 'Reactome'], default=['KEGG', 'GO:BP'])
        p_threshold = st.slider("P-value 阈值:", 0.01, 0.1, 0.05)
        run_enrich = st.button("📈 运行分析", type="primary")

    # 2. 分析与绘图逻辑
    if run_enrich and raw_text_enrich:
        gene_list = [x.strip() for x in raw_text_enrich.split('\n') if x.strip()]
        
        try:
            with st.spinner("正在进行富集分析 (g:Profiler)..."):
                # --- 修复点：初始化方式更新 ---
                gp = GProfiler(user_agent='streamlit_app')
                
                # --- 修复点：手动获取结果并转 DataFrame ---
                raw_results = gp.profile(organism=gprofiler_organism_code, query=gene_list, sources=enrich_sources, user_threshold=p_threshold, no_iea=True)
                
                # 检查是否有结果
                if not raw_results or 'result' not in raw_results or not raw_results['result']:
                    st.warning("⚠️ 未发现显著富集结果 (No significant pathways found)。")
                else:
                    # 将字典列表转换为 DataFrame
                    results = pd.DataFrame(raw_results['result'])
                    
                    # --- 数据预处理 ---
                    # 计算 -log10(pvalue)
                    results['neg_log10_p'] = results['p_value'].apply(lambda x: -math.log10(x))
                    # 截取过长名称
                    results['short_name'] = results['name'].apply(lambda x: x[:50] + '...' if len(x)>50 else x)
                    # 处理 intersections (列表转字符串)
                    if 'intersections' in results.columns:
                        results['intersections'] = results['intersections'].apply(lambda x: "; ".join(x) if isinstance(x, list) else str(x))
                    
                    # 排序
                    results = results.sort_values('p_value', ascending=True) 
                    
                    st.success(f"发现 {len(results)} 条显著通路！")
                    
                    # --- 可视化控制面板 ---
                    st.divider()
                    st.subheader("📊 可视化工作室")
                    
                    viz_col1, viz_col2 = st.columns([1, 3])
                    
                    with viz_col1:
                        with st.container(border=True):
                            st.markdown("### 🎨 绘图设置")
                            plot_type = st.selectbox("图表类型:", ["气泡图 (Dot Plot)", "柱状图 (Bar Chart)"])
                            top_n = st.slider("展示 Top N:", 5, 50, 20)
                            color_scale = st.selectbox("配色方案:", ["Tealgrn", "Viridis", "Plasma", "Bluered"])
                            base_font_size = st.slider("字体大小:", 8, 20, 12)
                            plot_height = st.slider("图表高度:", 400, 1000, 600)
                    
                    # 准备绘图数据 (Top N)
                    plot_data = results.head(top_n).copy()
                    plot_data = plot_data.sort_values('p_value', ascending=False) # 反转顺序使显著的在图上方
                    
                    with viz_col2:
                        fig = None
                        if plot_type == "气泡图 (Dot Plot)":
                            fig = px.scatter(
                                plot_data, x="intersection_size", y="short_name", size="intersection_size", 
                                color="neg_log10_p", hover_data=["p_value", "source"], color_continuous_scale=color_scale,
                                labels={"intersection_size": "Gene Count", "short_name": "Pathway", "neg_log10_p": "-log10(P-value)"},
                                title=f"Top {top_n} Enriched Pathways"
                            )
                        elif plot_type == "柱状图 (Bar Chart)":
                            fig = px.bar(
                                plot_data, x="neg_log10_p", y="short_name", color="intersection_size", orientation='h',
                                color_continuous_scale=color_scale,
                                labels={"neg_log10_p": "-log10(P-value)", "short_name": "Pathway", "intersection_size": "Count"},
                                title=f"Top {top_n} Enriched Pathways"
                            )
                        
                        fig.update_layout(height=plot_height, font=dict(size=base_font_size), plot_bgcolor='white')
                        st.plotly_chart(fig, use_container_width=True)
                        
                        # --- 导出 ---
                        st.markdown("### 📥 导出数据")
                        col_exp1, col_exp2 = st.columns(2)
                        
                        output_excel = io.BytesIO()
                        with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
                            # 剔除绘图辅助列再导出
                            export_cols = [c for c in results.columns if c not in ['neg_log10_p', 'short_name']]
                            results[export_cols].to_excel(writer, index=False)
                        col_exp1.download_button("下载数据表 (Excel)", output_excel.getvalue(), "enrichment_results.xlsx")
                        
                        buffer_html = io.StringIO()
                        fig.write_html(buffer_html)
                        col_exp2.download_button("下载交互式图表 (HTML)", buffer_html.getvalue().encode(), "enrichment_plot.html")

        except Exception as e:
            st.error(f"分析出错: {e}")
            st.markdown("建议：检查输入基因是否正确，或者网络连接是否正常。")
