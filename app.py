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
st.title("🔬 Omics Data Assistant (v5.1 - Layout Control)")

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
    s1 = set(list1)
    s2 = set(list2)
    if not s1 or not s2: return 0.0
    return len(s1 & s2) / len(s1 | s2)

def simplify_results(df, threshold=0.7):
    if df.empty: return df
    df = df.sort_values('p_value', ascending=True)
    keep_indices = []
    genes_list = df['intersections_raw'].tolist()
    
    for i in range(len(df)):
        current_genes = genes_list[i]
        is_redundant = False
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
# Tab 1: ID 转换
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
# Tab 2: 富集分析
# =================================================================================
with tab2:
    st.header("Enrichment Analysis")
    
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

    if run_enrich and raw_text_enrich:
        raw_gene_list = [x.strip() for x in raw_text_enrich.split('\n') if x.strip()]
        for key in ['raw_results', 'filtered_results']:
            if key in st.session_state: del st.session_state[key]

        with st.spinner("Talking to g:Profiler..."):
            try:
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
                        
                        def decode_intersections(inter_list):
                            if not isinstance(inter_list, list): return ""
                            hit_genes = []
                            for idx, evidences in enumerate(inter_list):
                                if evidences: 
                                    if idx < len(unique_converted_ids):
                                        hit_genes.append(entrez_to_symbol.get(unique_converted_ids[idx], unique_converted_ids[idx]))
                            return "; ".join(hit_genes)

                        def get_gene_ids_list(inter_list):
                            ids = []
                            for idx, evidences in enumerate(inter_list):
                                if evidences and idx < len(unique_converted_ids):
                                    ids.append(unique_converted_ids[idx])
                            return ids

                        if 'intersections' in results.columns:
                            results['hit_genes'] = results['intersections'].apply(decode_intersections)
                            results['intersections_raw'] = results['intersections'].apply(get_gene_ids_list)
                            results['intersections'] = results['hit_genes']

                        st.session_state['raw_results'] = results.sort_values('p_value')
                        st.success(f"Analysis Done. Found {len(results)} pathways.")
                    else:
                        st.warning("No significant pathways found.")

            except Exception as e:
                st.error(f"Error: {e}")

    if 'raw_results' in st.session_state:
        df_raw = st.session_state['raw_results']
        
        st.divider()
        st.subheader("Step 2: Filter & Select")
        
        col_f1, col_f2 = st.columns([1, 2])
        with col_f1:
            st.markdown("##### 1. Reduce Redundancy")
            use_simplify = st.checkbox("Apply Similarity Filter", value=False)
            sim_threshold = st.slider("Similarity Threshold:", 0.1, 1.0, 0.7)
            
        if use_simplify:
            df_processed = simplify_results(df_raw, threshold=sim_threshold)
            st.info(f"Filtered: {len(df_raw)} -> {len(df_processed)} pathways.")
        else:
            df_processed = df_raw.copy()

        with col_f2:
            st.markdown("##### 2. Manual Selection")
        
        df_display = df_processed[['source', 'name', 'p_value', 'intersection_size', 'hit_genes']].copy()
        df_display.insert(0, "Select", False)
        
        edited_df = st.data_editor(
            df_display,
            column_config={
                "Select": st.column_config.CheckboxColumn("Plot?", default=False),
                "p_value": st.column_config.NumberColumn(format="%.2e")
            },
            disabled=["source", "name", "p_value", "intersection_size", "hit_genes"],
            hide_index=True,
            height=300
        )
        
        selected_indices = edited_df[edited_df["Select"]].index
        
        if len(selected_indices) > 0:
            df_final_plot = df_processed.iloc[selected_indices].copy()
            auto_top_n = False
        else:
            df_final_plot = df_processed.copy()
            auto_top_n = True

        st.divider()
        st.subheader("Step 3: Visualization Studio")
        
        viz_c1, viz_c2 = st.columns([1, 3])
        with viz_c1:
            with st.container(border=True):
                st.markdown("**1. Chart Settings**")
                plot_type = st.selectbox("Chart Type:", ["Dot Plot", "Bar Chart"])
                
                if auto_top_n:
                    top_n = st.slider("Top N Pathways:", 5, 50, 20)
                    plot_data = df_final_plot.head(top_n).sort_values('p_value', ascending=False)
                    plot_title = f"Top {top_n} Enriched Pathways"
                else:
                    plot_data = df_final_plot.sort_values('p_value', ascending=False)
                    plot_title = "Manually Selected Pathways"

                color_scale = st.selectbox("Color Theme:", ["Tealgrn", "Viridis", "Plasma", "Bluered", "Sunset"])

                # --- 新增：排版控制 ---
                st.divider()
                st.markdown("**2. Layout Adjustments (Fix PDF Cutoff)**")
                margin_left = st.slider("Left Margin (px):", 100, 800, 300, help="Increase this if pathway names are cut off in the PDF.")
                plot_width = st.slider("Total Width (px):", 600, 2000, 1000, help="Increase width to prevent chart squishing.")
                plot_height = st.slider("Total Height (px):", 400, 1500, max(600, len(plot_data)*25))

        with viz_c2:
            # 在绘图前，如果名字实在太长，我们可以稍微截断一下，或者完全保留由用户通过 Margin 控制
            # 这里选择保留完整名字，因为现在有了 Margin 控制
            # plot_data['short_name'] = plot_data['name'].apply(lambda x: x[:50] + '...' if len(x)>50 else x)
            plot_data['short_name'] = plot_data['name'] # 尽可能保留全名
            
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
            
            # --- 关键修改：应用用户自定义的 Margin 和 Width ---
            fig.update_layout(
                width=plot_width,
                height=plot_height,
                margin=dict(l=margin_left, r=20, t=50, b=50), # l=Left Margin
                plot_bgcolor='white', 
                font=dict(family="Arial", size=12)
            )
            
            st.plotly_chart(fig, use_container_width=True) # 网页上自适应，但导出时会用上面的 width/height

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

        # PDF Export with explicit error handling
        try:
            with st.spinner("Generating PDF..."):
                # 使用用户定义的 width/height 导出
                pdf_bytes = fig.to_image(format="pdf", engine="kaleido", scale=2, width=plot_width, height=plot_height)
                
            e3.download_button(
                label="📥 Plot (PDF)",
                data=pdf_bytes,
                file_name="enrichment_plot.pdf",
                mime="application/pdf"
            )
        except Exception as e:
            e3.error(f"PDF Error: {e}")
