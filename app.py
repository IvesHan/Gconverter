import streamlit as st
import pandas as pd
import mygene
import io
import math
from gprofiler import GProfiler
import plotly.express as px
import plotly.graph_objects as go

# --- 页面基础设置 ---
st.set_page_config(page_title="BioInfo Tool Pro v3.2", layout="wide", page_icon="🧬")

st.title("🧬 基因组学多功能工具 (v3.2 - 宽松修正版)")
st.markdown("""
**R 语言对齐版：** 调整了统计算法 (FDR) 和注释范围，确保结果与 R/Web 端一致。
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
tab1, tab2 = st.tabs(["1. 基因 ID 转换与注释", "2. 富集分析与可视化 (Debug)"])

# =================================================================================
# Tab 1: 基因 ID 转换
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
        if not gene_list: st.warning("请输入基因 ID")
        else:
            with st.spinner("查询中..."):
                try:
                    mg = mygene.MyGeneInfo()
                    res = mg.querymany(gene_list, scopes='symbol,entrezgene,ensembl.gene,alias', fields=query_fields, species=species_id, as_dataframe=True)
                    df_res = res.reset_index()
                    for col in df_res.columns: df_res[col] = df_res[col].apply(clean_cell_data)
                    if input_method == "直接粘贴文本": final_df = df_res
                    else: final_df = pd.merge(df_input, df_res, left_on=col_name, right_on='query', how='left')
                    st.dataframe(final_df)
                    output = io.BytesIO()
                    with pd.ExcelWriter(output, engine='xlsxwriter') as writer: final_df.to_excel(writer, index=False)
                    st.download_button("📥 下载 Excel", output.getvalue(), "gene_conversion.xlsx")
                except Exception as e: st.error(f"Error: {e}")

# =================================================================================
# Tab 2: 富集分析 (深度修正版)
# =================================================================================
with tab2:
    st.header("功能二：富集分析 & 交互式绘图")
    st.markdown("如果之前跑不出结果，请尝试勾选 **'包含电子注释 (IEA)'** 并将矫正方法改为 **'FDR'**。")
    
    col_in1, col_in2 = st.columns([1, 2])
    with col_in1:
        raw_text_enrich = st.text_area("粘贴基因列表 (Symbol/Ensembl/Entrez):", height=150, placeholder="TP53\nEGFR\nCD4...", key="t2_text")
        
    with col_in2:
        with st.container(border=True):
            st.subheader("⚙️ 关键参数 (对齐 R 语言)")
            enrich_sources = st.multiselect("数据库:", ['KEGG', 'GO:BP', 'GO:CC', 'GO:MF', 'Reactome', 'WP'], default=['KEGG', 'GO:BP'])
            
            # --- 关键修正：参数调整 ---
            col_p1, col_p2 = st.columns(2)
            with col_p1:
                p_threshold = st.slider("P-value 阈值:", 0.01, 1.0, 0.05, help="R 默认通常为 0.05")
            with col_p2:
                correction_method = st.selectbox("矫正算法 (Correction):", 
                                                 ["fdr", "g_SCS", "bonferroni"], 
                                                 index=0, 
                                                 help="推荐使用 FDR (即 BH 方法)，这与 clusterProfiler 结果一致。g_SCS 是 g:Profiler 独有的严格算法。")
            
            exclude_iea = st.checkbox("排除电子注释 (Exclude IEA)", value=False, help="❌ 取消勾选（即包含 IEA）能获得更多 GO 结果。R 默认通常包含。")
            
            run_enrich = st.button("📈 运行分析 (Run Analysis)", type="primary")

    if run_enrich and raw_text_enrich:
        gene_list = [x.strip() for x in raw_text_enrich.split('\n') if x.strip()]
        
        try:
            with st.spinner("正在分析..."):
                gp = GProfiler(user_agent='streamlit_app_v3.2')
                
                # --- 1. 执行富集分析 ---
                raw_results = gp.profile(
                    organism=gprofiler_organism_code, 
                    query=gene_list, 
                    sources=enrich_sources, 
                    user_threshold=p_threshold, 
                    no_iea=exclude_iea,  # 这里的 False 意味着包含更多结果
                    significance_threshold_method=correction_method # 显式指定算法
                )
                
                # --- 2. 调试信息 (如果没结果，看看到底识别了几个基因) ---
                if not raw_results:
                    st.error("❌ API 返回为空。可能是网络问题或所有基因都未被识别。")
                elif not raw_results['result']:
                    st.warning(f"⚠️ 分析完成，但在当前阈值 (P<{p_threshold}, {correction_method}) 下未发现显著通路。")
                    
                    # 尝试调用 mapping 接口看看基因识别情况
                    st.info("🔎 正在诊断 ID 识别情况...")
                    try:
                        # 这是一个 hack，利用 convert 接口看识别率
                        mg_check = mygene.MyGeneInfo()
                        check_res = mg_check.querymany(gene_list[:10], scopes='symbol,entrezgene,ensembl.gene', fields='symbol', species=species_id)
                        found = sum([1 for x in check_res if 'notfound' not in x])
                        st.write(f"ID 诊断: 输入的前 10 个 ID 中，MyGene 成功识别了 {found} 个。如果这个数字很低，请检查物种选择或 ID 拼写。")
                    except:
                        pass
                else:
                    # --- 3. 处理成功结果 ---
                    results = pd.DataFrame(raw_results['result'])
                    
                    # 数据清洗
                    results['neg_log10_p'] = results['p_value'].apply(lambda x: -math.log10(x))
                    results['short_name'] = results['name'].apply(lambda x: x[:50] + '...' if len(x)>50 else x)
                    if 'intersections' in results.columns:
                        results['intersections'] = results['intersections'].apply(lambda x: "; ".join(x) if isinstance(x, list) else str(x))
                    
                    results = results.sort_values('p_value', ascending=True)
                    
                    st.success(f"✅ 成功发现 {len(results)} 条显著通路！(Top 50 shown below)")
                    
                    # --- 可视化 ---
                    st.divider()
                    viz_col1, viz_col2 = st.columns([1, 3])
                    
                    with viz_col1:
                        with st.container(border=True):
                            st.markdown("### 🎨 绘图控制")
                            plot_type = st.selectbox("图表类型:", ["气泡图 (Dot Plot)", "柱状图 (Bar Chart)"])
                            top_n = st.slider("展示数量:", 5, 50, 20)
                            color_scale = st.selectbox("配色:", ["Tealgrn", "Viridis", "Plasma", "Bluered"])
                    
                    plot_data = results.head(top_n).copy().sort_values('p_value', ascending=False)
                    
                    with viz_col2:
                        if plot_type == "气泡图 (Dot Plot)":
                            fig = px.scatter(
                                plot_data, x="intersection_size", y="short_name", size="intersection_size", 
                                color="neg_log10_p", hover_data=["p_value", "source"], color_continuous_scale=color_scale,
                                labels={"intersection_size": "Count", "short_name": "Pathway", "neg_log10_p": "-log10(P)"},
                                title=f"Top {top_n} Enriched Pathways ({correction_method})"
                            )
                        else:
                            fig = px.bar(
                                plot_data, x="neg_log10_p", y="short_name", color="intersection_size", orientation='h',
                                color_continuous_scale=color_scale,
                                labels={"neg_log10_p": "-log10(P)", "short_name": "Pathway", "intersection_size": "Count"},
                                title=f"Top {top_n} Enriched Pathways ({correction_method})"
                            )
                        
                        fig.update_layout(height=600, plot_bgcolor='white')
                        st.plotly_chart(fig, use_container_width=True)
                        
                        # 导出
                        col_e1, col_e2 = st.columns(2)
                        output_excel = io.BytesIO()
                        with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
                            results.drop(columns=['neg_log10_p', 'short_name']).to_excel(writer, index=False)
                        col_e1.download_button("📥 下载 Excel", output_excel.getvalue(), "enrichment.xlsx")
                        
                        buf = io.StringIO()
                        fig.write_html(buf)
                        col_e2.download_button("📥 下载 HTML 图表", buf.getvalue().encode(), "plot.html")

        except Exception as e:
            st.error(f"Error: {e}")
