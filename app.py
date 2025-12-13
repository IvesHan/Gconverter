import streamlit as st
import pandas as pd
import mygene
import io
import math
import requests
import plotly.express as px
import plotly.graph_objects as go

# --- 1. 页面配置 (低调专业风格) ---
st.set_page_config(page_title="Omics Analysis Tool", layout="wide", page_icon="🔬")

st.title("🔬 Omics Data Assistant")
# st.markdown("Data processing and visualization tool.") # 甚至可以不写描述，保持极简

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
    """Excel 导出格式清洗"""
    if isinstance(cell, list):
        cleaned_list = [str(item) if not isinstance(item, dict) else f"{item.get('chr','N/A')}:{item.get('start','N/A')}" for item in cell]
        return "; ".join(cleaned_list)
    return str(cell) if isinstance(cell, dict) else cell

# --- 4. 页面布局 ---
tab1, tab2 = st.tabs(["ID Conversion", "Enrichment & Visualization"])

# =================================================================================
# Tab 1: ID 转换 (保持功能，界面简化)
# =================================================================================
with tab1:
    st.header("ID Conversion")
    
    # 简单的两列布局
    c1, c2 = st.columns([1, 1])
    with c1:
        input_method = st.radio("Input Method:", ("Paste Text", "Upload File"))
    with c2:
        target_fields = st.multiselect("Additional Info:", 
                                     ['name', 'alias', 'summary', 'type_of_gene', 'genomic_pos'], 
                                     default=['name'])

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
# Tab 2: 富集分析 (加入 Session State 状态保持)
# =================================================================================
with tab2:
    st.header("Enrichment Analysis")
    
    # 2.1 输入与参数区
    with st.expander("Analysis Parameters", expanded=True):
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

    # 2.2 分析逻辑 (只在点击按钮时触发 API)
    if run_enrich and raw_text_enrich:
        raw_gene_list = [x.strip() for x in raw_text_enrich.split('\n') if x.strip()]
        
        # 清除旧结果
        if 'enrich_data' in st.session_state:
            del st.session_state['enrich_data']

        with st.spinner("Analyzing..."):
            try:
                # Step 1: MyGene ID Mapping
                mg = mygene.MyGeneInfo()
                map_res = mg.querymany(raw_gene_list, scopes='symbol,entrezgene,ensembl.gene,alias', fields='entrezgene', species=species_id)
                converted_ids = []
                for item in map_res:
                    if 'entrezgene' in item: converted_ids.append(str(item['entrezgene']))
                converted_ids = list(set(converted_ids))

                if not converted_ids:
                    st.error("No valid IDs identified.")
                else:
                    # Step 2: Direct API Call
                    payload = {
                        'organism': gprofiler_organism_code,
                        'query': converted_ids,
                        'sources': enrich_sources,
                        'user_threshold': p_threshold,
                        'no_iea': exclude_iea,
                        'significance_threshold_method': correction_method,
                        'numeric_ns': 'ENTREZGENE_ACC'
                    }
                    response = requests.post('https://biit.cs.ut.ee/gprofiler/api/gost/profile/', json=payload)
                    raw_results = response.json()

                    # Step 3: Data Processing
                    if 'result' in raw_results and raw_results['result']:
                        results = pd.DataFrame(raw_results['result'])
                        results['neg_log10_p'] = results['p_value'].apply(lambda x: -math.log10(x))
                        results['short_name'] = results['name'].apply(lambda x: x[:50] + '...' if len(x)>50 else x)
                        if 'intersections' in results.columns:
                            # 修复 list join 报错的关键代码
                            results['intersections'] = results['intersections'].apply(lambda x: "; ".join(map(str, x)) if isinstance(x, list) else str(x))
                        
                        # 保存到 Session State (持久化存储)
                        st.session_state['enrich_data'] = results.sort_values('p_value')
                        st.success(f"Analysis Done. Found {len(results)} pathways.")
                    else:
                        st.warning("No significant pathways found.")

            except Exception as e:
                st.error(f"Error: {e}")

    # 2.3 可视化区域 (独立于按钮，只要 Session State 里有数据就显示)
    if 'enrich_data' in st.session_state:
        df = st.session_state['enrich_data']
        
        st.divider()
        st.subheader("Visualization Studio")
        
        # 绘图控制面板
        viz_c1, viz_c2 = st.columns([1, 3])
        
        with viz_c1:
            with st.container(border=True):
                st.markdown("**Chart Settings**")
                # 这些控件的改动会触发页面刷新，但因为数据在 session_state 里，图表会立刻重绘
                plot_type = st.selectbox("Chart Type:", ["Dot Plot", "Bar Chart"])
                top_n = st.slider("Top N Pathways:", 5, 50, 20)
                color_scale = st.selectbox("Color Theme:", ["Tealgrn", "Viridis", "Plasma", "Bluered", "Sunset"])
                
                # 数据筛选
                plot_data = df.head(top_n).copy().sort_values('p_value', ascending=False)

        with viz_c2:
            if plot_type == "Dot Plot":
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
            
            fig.update_layout(height=600, plot_bgcolor='white', font=dict(family="Arial", size=12))
            st.plotly_chart(fig, use_container_width=True)

        # 导出区域
        st.markdown("### Export")
        e1, e2 = st.columns(2)
        
        # Excel 导出
        out_df = df.drop(columns=['neg_log10_p', 'short_name'], errors='ignore')
        output_excel = io.BytesIO()
        with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
            out_df.to_excel(writer, index=False)
        e1.download_button("📥 Download Data (Excel)", output_excel.getvalue(), "enrichment_results.xlsx")
        
        # 图片导出
        buf = io.StringIO()
        fig.write_html(buf)
        e2.download_button("📥 Download Plot (HTML)", buf.getvalue().encode(), "plot.html")
