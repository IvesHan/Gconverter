import streamlit as st
import pandas as pd
import mygene
import io
import math
import requests
import plotly.express as px
import plotly.graph_objects as go

# --- 1. 页面配置 ---
st.set_page_config(page_title="Omics Analysis Tool", layout="wide", page_icon="🔬")
st.title("🔬 Omics Data Assistant (v5.0 - Custom Filter)")

# --- 2. 全局物种映射 ---
species_map = {
    "Human (Homo sapiens)": (9606, 'hsapiens'),
    "Mouse (Mus musculus)": (10090, 'mmusculus'),
    "Rat (Rattus norvegicus)": (10116, 'rnorvegicus')
}

st.sidebar.header("Settings")
selected_species_key = st.sidebar.selectbox("Species:", options=list(species_map.keys()))
species_id, gprofiler_organism_code = species_map[selected_species_key]

# --- 3. 辅助函数 ---
def clean_cell_data(cell):
    if isinstance(cell, list):
        cleaned_list = [str(item) if not isinstance(item, dict) else f"{item.get('chr','N/A')}:{item.get('start','N/A')}" for item in cell]
        return "; ".join(cleaned_list)
    return str(cell) if isinstance(cell, dict) else cell

def calculate_jaccard(list1, list2):
    """计算两个列表的 Jaccard 相似系数"""
    s1 = set(list1)
    s2 = set(list2)
    if not s1 or not s2: return 0.0
    return len(s1 & s2) / len(s1 | s2)

def simplify_results(df, threshold=0.7):
    """
    仿照 clusterProfiler 的 simplify 逻辑：
    1. 按 p值排序
    2. 遍历，如果当前通路与已保留通路的基因重叠度 > 阈值，则丢弃
    """
    if df.empty: return df
    
    # 确保按显著性排序
    df = df.sort_values('p_value', ascending=True)
    keep_indices = []
    
    # 获取所有基因列表 (注意：这里需要原始的 intersection ID 列表)
    # 我们假设 df['intersections_raw'] 存储了原始 ID 列表
    genes_list = df['intersections_raw'].tolist()
    
    for i in range(len(df)):
        current_genes = genes_list[i]
        is_redundant = False
        
        # 与已经保留的通路进行比较
        for kept_idx in keep_indices:
            kept_genes = genes_list[kept_idx]
            sim = calculate_jaccard(current_genes, kept_genes)
            if sim > threshold:
                is_redundant = True
                break
        
        if not is_redundant:
            keep_indices.append(i)
            
    return df.iloc[keep_indices]

# --- 4. 页面布局 ---
tab1, tab2 = st.tabs(["ID Conversion", "Enrichment & Visualization"])

# =================================================================================
# Tab 1: ID 转换 (保持不变)
# =================================================================================
with tab1:
    st.header("ID Conversion")
    c1, c2 = st.columns([1, 1])
    with c1:
        input_method = st.radio("Input Method:", ("Paste Text", "Upload File"))
    with c2:
        target_fields = st.multiselect("Additional Info:", ['name', 'alias', 'summary', 'type_of_gene', 'genomic_pos'], default=['name'])

    df_input = None
    gene_list = []
    col_name = "Input_ID"

    if input_method == "Paste Text":
        raw_text = st.text_area("Gene List (One per line):", height=150)
        if raw_text:
            gene_list = [x.strip() for x in raw_text.split('\n') if x.strip()]
            df_input = pd.DataFrame({col_name: gene_list})
    else:
        uploaded_file = st.file_uploader("Upload Excel/CSV", type=['xlsx', 'csv'])
        if uploaded_file:
            if uploaded_file.name.endswith('.csv'): df_input = pd.read_csv(uploaded_file)
            else: df_input = pd.read_excel(uploaded_file)
            col_name = st.selectbox("Select ID Column:", df_input.columns)
            gene_list = df_input[col_name].dropna().astype(str).tolist()

    if st.button("Convert IDs"):
        if not gene_list:
            st.warning("Please input gene IDs.")
        else:
            with st.spinner("Processing..."):
                try:
                    mg = mygene.MyGeneInfo()
                    fields = ['symbol', 'entrezgene', 'ensembl.gene'] + target_fields
                    res = mg.querymany(gene_list, scopes='symbol,entrezgene,ensembl.gene,alias', fields=fields, species=species_id, as_dataframe=True)
                    df_res = res.reset_index()
                    for col in df_res.columns: df_res[col] = df_res[col].apply(clean_cell_data)
                    
                    if input_method == "Paste Text": final_df = df_res
                    else: final_df = pd.merge(df_input, df_res, left_on=col_name, right_on='query', how='left')
                    
                    st.dataframe(final_df)
                    output = io.BytesIO()
                    with pd.ExcelWriter(output, engine='xlsxwriter') as writer: final_df.to_excel(writer, index=False)
                    st.download_button("Download Result", output.getvalue(), "conversion_result.xlsx")
                except Exception as e:
                    st.error(f"Error: {e}")

# =================================================================================
# Tab 2: 富集分析 (新增 Filter 和 Manual Select)
# =================================================================================
with tab2:
    st.header("Enrichment Analysis")
    
    # 2.1 参数区
    with st.expander("Step 1: Analysis Parameters", expanded=True):
        col_in1, col_in2 = st.columns([1, 2])
        with col_in1:
            raw_text_enrich = st.text_area("Paste Gene List:", height=150, placeholder="TP53\nEGFR...")
        
        with col_in2:
            enrich_sources = st.multiselect("Databases:", ['KEGG', 'GO:BP', 'GO:CC', 'GO:MF', 'Reactome'], default=['KEGG', 'GO:BP'])
            c_p1, c_p2, c_p3 = st.columns(3)
            p_threshold = c_p1.slider("P-value Cutoff:", 0.01, 1.0, 0.05)
            correction_method = c_p2.selectbox("Correction:", ["fdr", "bonferroni", "g_SCS"], index=0)
            exclude_iea = c_p3.checkbox("No IEA", value=False)
            run_enrich = st.button("Run Analysis", type="primary")

    # 2.2 API 分析逻辑
    if run_enrich and raw_text_enrich:
        raw_gene_list = [x.strip() for x in raw_text_enrich.split('\n') if x.strip()]
        
        # 清除缓存
        for key in ['raw_results', 'filtered_results']:
            if key in st.session_state: del st.session_state[key]

        with st.spinner("Talking to g:Profiler..."):
            try:
                # ID Mapping
                mg = mygene.MyGeneInfo()
                map_res = mg.querymany(raw_gene_list, scopes='symbol,entrezgene,ensembl.gene,alias', fields='entrezgene,symbol', species=species_id)
                converted_ids = []
                entrez_to_symbol = {}
                for item in map_res:
                    if 'entrezgene' in item:
                        eid = str(item['entrezgene'])
                        sym = item.get('symbol', item.get('query', eid))
                        converted_ids.append(eid)
                        entrez_to_symbol[eid] = sym
                unique_converted_ids = list(set(converted_ids))

                if not unique_converted_ids:
                    st.error("No valid IDs identified.")
                else:
                    # API Call
                    payload = {
                        'organism': gprofiler_organism_code,
                        'query': unique_converted_ids,
                        'sources': enrich_sources,
                        'user_threshold': p_threshold,
                        'no_iea': exclude_iea,
                        'significance_threshold_method': correction_method,
                        'numeric_ns': 'ENTREZGENE_ACC'
                    }
                    response = requests.post('https://biit.cs.ut.ee/gprofiler/api/gost/profile/', json=payload)
                    raw_results = response.json()

                    if 'result' in raw_results and raw_results['result']:
                        results = pd.DataFrame(raw_results['result'])
                        results['neg_log10_p'] = results['p_value'].apply(lambda x: -math.log10(x))
                        
                        # 解析 Genes
                        def decode_intersections(inter_list):
                            if not isinstance(inter_list, list): return ""
                            hit_genes = []
                            for idx, evidences in enumerate(inter_list):
                                if evidences: 
                                    if idx < len(unique_converted_ids):
                                        hit_genes.append(entrez_to_symbol.get(unique_converted_ids[idx], unique_converted_ids[idx]))
                            return "; ".join(hit_genes)

                        # 保存原始的 intersections 列表用于 Jaccard 计算 (intersections_raw)
                        # API 返回的 intersections 实际上是列表的列表 (证据代码)，我们需要根据这个位置反推 ID
                        # 为了简化 Jaccard，我们只关心“哪些基因ID在这个通路里”
                        # 构造一个真实的 Entrez ID 列表给每一行
                        
                        def get_gene_ids_list(inter_list):
                            ids = []
                            for idx, evidences in enumerate(inter_list):
                                if evidences and idx < len(unique_converted_ids):
                                    ids.append(unique_converted_ids[idx])
                            return ids

                        if 'intersections' in results.columns:
                            results['hit_genes'] = results['intersections'].apply(decode_intersections)
                            results['intersections_raw'] = results['intersections'].apply(get_gene_ids_list)
                            # 显示用的列
                            results['intersections'] = results['hit_genes']

                        # 存入 Session State (Raw Data)
                        st.session_state['raw_results'] = results.sort_values('p_value')
                        st.success(f"Analysis Done. Found {len(results)} pathways.")
                    else:
                        st.warning("No significant pathways found.")

            except Exception as e:
                st.error(f"Error: {e}")

    # 2.3 过滤与手动选择区域 (独立处理)
    if 'raw_results' in st.session_state:
        df_raw = st.session_state['raw_results']
        
        st.divider()
        st.subheader("Step 2: Filter & Select")
        
        # --- 新功能 1: 相似性过滤 (Simplify) ---
        col_f1, col_f2 = st.columns([1, 2])
        with col_f1:
            st.markdown("##### 1. Reduce Redundancy (Simplify)")
            use_simplify = st.checkbox("Apply Similarity Filter", value=False, help="Remove pathways that share too many genes with more significant pathways.")
            sim_threshold = st.slider("Similarity Threshold (Jaccard):", 0.1, 1.0, 0.7, help="Lower value = More strict filtering (fewer pathways).")
            
        # 计算过滤后的数据
        if use_simplify:
            df_processed = simplify_results(df_raw, threshold=sim_threshold)
            st.info(f"Filtered: {len(df_raw)} -> {len(df_processed)} pathways.")
        else:
            df_processed = df_raw.copy()

        # --- 新功能 2: 手动选择 (Data Editor) ---
        with col_f2:
            st.markdown("##### 2. Manual Selection")
            st.markdown("Check the box to include specific pathways in the plot. Uncheck all to auto-select Top N.")
        
        # 构建一个带勾选框的 DataFrame
        df_display = df_processed[['source', 'name', 'p_value', 'intersection_size', 'hit_genes']].copy()
        df_display.insert(0, "Select", False) # 插入第一列为复选框
        
        # 使用 data_editor 让用户操作
        edited_df = st.data_editor(
            df_display,
            column_config={
                "Select": st.column_config.CheckboxColumn(
                    "Plot?",
                    help="Select pathways to visualize",
                    default=False,
                ),
                "p_value": st.column_config.NumberColumn(format="%.2e")
            },
            disabled=["source", "name", "p_value", "intersection_size", "hit_genes"],
            hide_index=True,
            height=300
        )
        
        # 获取用户选中的行
        selected_indices = edited_df[edited_df["Select"]].index
        
        if len(selected_indices) > 0:
            # 如果用户手动选了，就只画选中的
            df_final_plot = df_processed.iloc[selected_indices].copy()
            auto_top_n = False
        else:
            # 如果用户没选，就走默认 Top N
            df_final_plot = df_processed.copy()
            auto_top_n = True

        # --- Visualization Studio ---
        st.divider()
        st.subheader("Step 3: Visualization Studio")
        
        viz_c1, viz_c2 = st.columns([1, 3])
        with viz_c1:
            with st.container(border=True):
                st.markdown("**Chart Settings**")
                plot_type = st.selectbox("Chart Type:", ["Dot Plot", "Bar Chart"])
                
                # 如果是自动 Top N，才显示滑块
                if auto_top_n:
                    top_n = st.slider("Top N Pathways:", 5, 50, 20)
                    plot_data = df_final_plot.head(top_n).sort_values('p_value', ascending=False)
                    plot_title = f"Top {top_n} Enriched Pathways"
                else:
                    st.success(f"Plotting {len(df_final_plot)} manually selected pathways.")
                    plot_data = df_final_plot.sort_values('p_value', ascending=False)
                    plot_title = "Manually Selected Pathways"

                color_scale = st.selectbox("Color Theme:", ["Tealgrn", "Viridis", "Plasma", "Bluered", "Sunset"])

        with viz_c2:
            # 缩短名字用于绘图
            plot_data['short_name'] = plot_data['name'].apply(lambda x: x[:50] + '...' if len(x)>50 else x)
            
            if plot_type == "Dot Plot":
                fig = px.scatter(
                    plot_data, x="intersection_size", y="short_name", size="intersection_size", 
                    color="neg_log10_p", hover_data=["p_value", "source", "hit_genes"],
                    color_continuous_scale=color_scale,
                    labels={"intersection_size": "Count", "short_name": "Pathway", "neg_log10_p": "-log10(P)"},
                    title=plot_title
                )
            else:
                fig = px.bar(
                    plot_data, x="neg_log10_p", y="short_name", color="intersection_size", orientation='h',
                    color_continuous_scale=color_scale, hover_data=["hit_genes"],
                    labels={"neg_log10_p": "-log10(P)", "short_name": "Pathway", "intersection_size": "Count"},
                    title=plot_title
                )
            
            fig.update_layout(height=max(500, len(plot_data)*25), plot_bgcolor='white', font=dict(family="Arial", size=12))
            st.plotly_chart(fig, use_container_width=True)

        # Export
        st.markdown("### Export")
        e1, e2, e3 = st.columns(3)
        
        out_df = df_final_plot.drop(columns=['neg_log10_p', 'intersections_raw'], errors='ignore')
        output_excel = io.BytesIO()
        with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
            out_df.to_excel(writer, index=False)
        e1.download_button("📥 Data (Excel)", output_excel.getvalue(), "enrichment.xlsx")
        
        buf_html = io.StringIO()
        fig.write_html(buf_html)
        e2.download_button("📥 Plot (HTML)", buf_html.getvalue().encode(), "plot.html")

        try:
            pdf_bytes = fig.to_image(format="pdf", engine="kaleido")
            e3.download_button("📥 Plot (PDF)", pdf_bytes, "enrichment_plot.pdf", "application/pdf")
        except:
            e3.error("PDF export requires 'kaleido'.")
