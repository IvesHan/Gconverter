import streamlit as st
import pandas as pd
import mygene
import io
import math
import requests # <--- 弃用 gprofiler 库，直接用原生请求
import plotly.express as px
import plotly.graph_objects as go

# --- 页面基础设置 ---
st.set_page_config(page_title="BioInfo Tool Pro v3.7", layout="wide", page_icon="🧬")

st.title("🧬 基因组学多功能工具 (v3.7 - 核武版)")
st.markdown("""
**修复策略：** 弃用第三方 Python 库，直接调用 g:Profiler 官方 API 接口。强制指定 `numeric_ns='ENTREZGENE_ACC'`，确保纯数字 ID 被 100% 识别。
""")

# --- 全局物种映射 ---
species_map = {
    "Human (Homo sapiens)": (9606, 'hsapiens'),
    "Mouse (Mus musculus)": (10090, 'mmusculus'),
    "Rat (Rattus norvegicus)": (10116, 'rnorvegicus')
}

st.sidebar.header("🛠️ 全局设置")
selected_species_key = st.sidebar.selectbox("选择物种:", options=list(species_map.keys()))
species_id, gprofiler_organism_code = species_map[selected_species_key]

# --- 辅助函数 ---
def clean_cell_data(cell):
    if isinstance(cell, list):
        cleaned_list = [str(item) if not isinstance(item, dict) else f"{item.get('chr','N/A')}:{item.get('start','N/A')}" for item in cell]
        return "; ".join(cleaned_list)
    return str(cell) if isinstance(cell, dict) else cell

# --- Tabs ---
tab1, tab2 = st.tabs(["1. 基因 ID 转换", "2. 富集分析 (Direct API)"])

# =================================================================================
# Tab 1: 基因 ID 转换 (保持不变，避免引入新错误)
# =================================================================================
with tab1:
    # --- (此处代码逻辑与 v3.2 相同，为节省篇幅，省略但请在实际文件中保留) ---
    st.markdown("ID 转换功能代码保持不变...")
    
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
# Tab 2: 富集分析 (修复 TypeError)
# =================================================================================
with tab2:
    st.header("功能二：富集分析 & 可视化")
    st.markdown("✅ **技术说明：** Direct API 模式 + 数据类型强制转换。")
    
    col_in1, col_in2 = st.columns([1, 2])
    with col_in1:
        raw_text_enrich = st.text_area("粘贴基因列表:", height=200, placeholder="TP53\nEGFR...", key="t2_text")
        
    with col_in2:
        with st.container(border=True):
            st.subheader("⚙️ 参数设置")
            enrich_sources = st.multiselect("数据库:", ['KEGG', 'GO:BP', 'GO:CC', 'GO:MF', 'Reactome'], default=['KEGG', 'GO:BP'])
            
            col_p1, col_p2 = st.columns(2)
            with col_p1:
                p_threshold = st.slider("P-value 阈值:", 0.01, 1.0, 0.05)
            with col_p2:
                correction_map = {"fdr": "fdr", "bonferroni": "bonferroni", "g_SCS": "g_SCS"}
                correction_method = st.selectbox("矫正算法:", list(correction_map.keys()), index=0)
            
            exclude_iea = st.checkbox("排除电子注释 (建议不勾选)", value=False)
            
            run_enrich = st.button("📈 运行分析", type="primary")

    if run_enrich and raw_text_enrich:
        raw_gene_list = [x.strip() for x in raw_text_enrich.split('\n') if x.strip()]
        
        # --- 步骤 1: 转换为 Entrez ID ---
        with st.spinner("第一步: MyGene 获取 Entrez ID..."):
            try:
                mg = mygene.MyGeneInfo()
                map_res = mg.querymany(raw_gene_list, scopes='symbol,entrezgene,ensembl.gene,alias', fields='entrezgene', species=species_id)
                
                converted_ids = []
                for item in map_res:
                    if 'entrezgene' in item:
                        converted_ids.append(str(item['entrezgene']))
                converted_ids = list(set(converted_ids))
                
            except Exception as e:
                st.error(f"ID 转换失败: {e}")
                st.stop()

        # 调试面板
        with st.expander(f"🔍 ID 转换日志 (获取到 {len(converted_ids)} 个唯一 Entrez ID)", expanded=True):
            st.text(f"发送给 API 的 ID: {converted_ids[:10]} ...")
            if not converted_ids:
                st.error("❌ 无法获取 Entrez ID，请检查物种或输入。")
                st.stop()

        # --- 步骤 2: 原生 API 调用 ---
        with st.spinner("第二步: 调用 g:Profiler 官方 API..."):
            try:
                payload = {
                    'organism': gprofiler_organism_code,
                    'query': converted_ids,
                    'sources': enrich_sources,
                    'user_threshold': p_threshold,
                    'no_iea': exclude_iea,
                    'significance_threshold_method': correction_method,
                    'numeric_ns': 'ENTREZGENE_ACC'
                }
                
                response = requests.post(
                    'https://biit.cs.ut.ee/gprofiler/api/gost/profile/', 
                    json=payload
                )
                
                if response.status_code != 200:
                    st.error(f"服务器连接失败 (Status {response.status_code}): {response.text}")
                    st.stop()
                    
                raw_results = response.json()
                
            except Exception as e:
                st.error(f"API 通讯错误: {e}")
                st.stop()

        # --- 步骤 3: 结果处理 (关键修复点) ---
        if 'result' not in raw_results or not raw_results['result']:
            st.warning(f"⚠️ 分析完成，未发现显著通路。")
        else:
            # 成功！
            results = pd.DataFrame(raw_results['result'])
            
            # 1. 计算 -log10 P值
            results['neg_log10_p'] = results['p_value'].apply(lambda x: -math.log10(x))
            
            # 2. 截断过长名称
            results['short_name'] = results['name'].apply(lambda x: x[:50] + '...' if len(x)>50 else x)
            
            # 3. 修复 intersections 列 (TypeError 修复点)
            if 'intersections' in results.columns:
                # 这里的 map(str, x) 是核心修复：把数字列表转为字符串列表
                results['intersections'] = results['intersections'].apply(
                    lambda x: "; ".join(map(str, x)) if isinstance(x, list) else str(x)
                )
            
            # 4. 排序与筛选列
            display_df = results[[
                'source', 'native', 'name', 'p_value', 'intersection_size', 'term_size', 'intersections', 'neg_log10_p', 'short_name'
            ]].sort_values('p_value')
            
            st.success(f"✅ 成功发现 {len(display_df)} 条通路！(分析完成)")
            
            # --- 可视化 ---
            st.divider()
            viz_col1, viz_col2 = st.columns([1, 3])
            
            with viz_col1:
                with st.container(border=True):
                    st.markdown("### 🎨 绘图控制")
                    plot_type = st.selectbox("图表类型:", ["气泡图 (Dot Plot)", "柱状图 (Bar Chart)"])
                    top_n = st.slider("展示数量:", 5, 50, 20)
                    color_scale = st.selectbox("配色:", ["Tealgrn", "Viridis", "Plasma", "Bluered"])
            
            plot_data = display_df.head(top_n).copy().sort_values('p_value', ascending=False)
            
            with viz_col2:
                if plot_type == "气泡图 (Dot Plot)":
                    fig = px.scatter(
                        plot_data, x="intersection_size", y="short_name", size="intersection_size", 
                        color="neg_log10_p", hover_data=["p_value", "source"], color_continuous_scale=color_scale,
                        labels={"intersection_size": "Count", "short_name": "Pathway", "neg_log10_p": "-log10(P)"},
                        title=f"Top {top_n} Enriched Pathways"
                    )
                else:
                    fig = px.bar(
                        plot_data, x="neg_log10_p", y="short_name", color="intersection_size", orientation='h',
                        color_continuous_scale=color_scale,
                        labels={"neg_log10_p": "-log10(P)", "short_name": "Pathway", "intersection_size": "Count"},
                        title=f"Top {top_n} Enriched Pathways"
                    )
                
                fig.update_layout(height=600, plot_bgcolor='white')
                st.plotly_chart(fig, use_container_width=True)
                
                col_e1, col_e2 = st.columns(2)
                output_excel = io.BytesIO()
                with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
                    display_df.drop(columns=['neg_log10_p', 'short_name']).to_excel(writer, index=False)
                col_e1.download_button("📥 下载 Excel", output_excel.getvalue(), "enrichment.xlsx")
                
                buf = io.StringIO()
                fig.write_html(buf)
                col_e2.download_button("📥 下载 HTML", buf.getvalue().encode(), "plot.html")
