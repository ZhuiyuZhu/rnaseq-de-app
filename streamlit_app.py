"""RNApipe Web — 两组 RNA-seq 差异表达分析 GUI (Streamlit).

本地运行:   streamlit run streamlit_app.py
云端部署:   https://share.streamlit.io 绑定本仓库, main file 填 streamlit_app.py
"""
from pathlib import Path
import logging
import tempfile

import pandas as pd
import streamlit as st

from rnaseq_de.pipeline import run as run_pipeline

st.set_page_config(page_title="RNA-seq DE 分析", page_icon="🧬", layout="wide")
st.title("🧬 RNA-seq 两组差异表达分析")
st.caption("fastq → k-mer 伪比对定量 → 差异分析 → 交互式结果 | rnaseq-de-app v0.1.0")


class _ListHandler(logging.Handler):
    """捕获流水线日志, 在页面上显示进度。"""
    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        self.records.append(self.format(record))


# ---------------- 侧边栏: 输入与参数 ----------------
with st.sidebar:
    st.header("① 输入数据")
    use_example = st.checkbox("使用内置示例数据（推荐先跑一遍）", value=True)
    up_fq = up_fa = up_ss = None
    cond_a = cond_b = None
    if not use_example:
        up_fq = st.file_uploader("fastq 文件（可多选，支持 .gz）",
                                 accept_multiple_files=True)
        up_fa = st.file_uploader("参考转录本 fasta（header: >tx_xxx gene=yyy）",
                                 type=["fa", "fasta", "gz"])
        up_ss = st.file_uploader("样本分组表 CSV", type=["csv"])
        st.download_button("下载分组表模板",
                           "sample,fastq_1,fastq_2,condition\n",
                           "samplesheet_template.csv")
        if up_ss is not None:
            _ss = pd.read_csv(up_ss)
            _conds = sorted(_ss["condition"].dropna().unique())
            up_ss.seek(0)
            if len(_conds) == 2:
                cond_a = st.selectbox("基准组 A（对照）", _conds, key="ca")
                cond_b = st.selectbox("处理组 B（case）",
                                      [c for c in _conds if c != cond_a],
                                      key="cb")
            else:
                st.warning(f"分组表需恰好 2 组，当前 {_conds}")
    else:
        cond_a, cond_b = "control", "case"

    st.header("② 参数")
    de_backend = st.selectbox("差异分析后端",
                              ["auto（有 pydeseq2 用 DESeq2）", "ttest", "pydeseq2"])
    de_backend = de_backend.split("（")[0]
    alpha = st.slider("显著性阈值 padj", 0.01, 0.20, 0.10, 0.01)
    lfc = st.slider("最小 |log2FC|（图中标注）", 0.0, 3.0, 1.0, 0.5)
    run_clicked = st.button("🚀 运行分析", type="primary",
                            use_container_width=True)


# ---------------- 运行 ----------------
def _save_uploads(files, dest: Path) -> set:
    dest.mkdir(parents=True, exist_ok=True)
    names = set()
    for f in files:
        (dest / f.name).write_bytes(f.getbuffer())
        names.add(f.name)
    return names


if run_clicked:
    handler = _ListHandler()
    handler.setFormatter(logging.Formatter("%%(message)s".replace("%%", "%")))
    logging.getLogger("rnaseq_de").addHandler(handler)
    try:
        with tempfile.TemporaryDirectory() as td:
            work = Path(td)
            if use_example:
                root = Path(__file__).resolve().parent
                fastq_dir = root / "data" / "fastq"
                ss_path = root / "data" / "samplesheet.csv"
                fa_path = root / "data" / "ref" / "transcriptome.fa"
            else:
                if not (up_fq and up_fa and up_ss):
                    st.error("请先在侧边栏上传 fastq、转录本 fasta 和分组表。")
                    st.stop()
                ss = pd.read_csv(up_ss)
                up_ss.seek(0)
                names = _save_uploads(up_fq, work / "fastq")
                ref = list(ss["fastq_1"]) + list(ss["fastq_2"])
                missing = [Path(str(x)).name for x in ref
                           if Path(str(x)).name not in names]
                if missing:
                    st.error("分组表引用但未上传的文件: " + ", ".join(missing[:8]))
                    st.stop()
                fa_path = work / "transcriptome.fa"
                fa_path.write_bytes(up_fa.getbuffer())
                ss_path = work / "samplesheet.csv"
                ss_path.write_bytes(up_ss.getbuffer())
                fastq_dir = work / "fastq"

            outdir = work / "results"
            with st.status("正在运行分析…", expanded=True) as status:
                res_df, norm_df = run_pipeline(
                    str(fastq_dir), str(ss_path), str(fa_path), str(outdir),
                    quant_backend="python", de_backend=de_backend,
                    cond_a=cond_a, cond_b=cond_b, alpha=alpha, lfc=lfc)
                status.update(label="✅ 分析完成", state="complete")
                st.code("\n".join(handler.records), language=None)

            st.session_state["de"] = res_df.reset_index(names="gene")
            st.session_state["figs"] = {p.stem: p.read_bytes()
                                        for p in sorted((outdir / "figures").glob("*.png"))}
            st.session_state["qc"] = pd.read_csv(outdir / "qc_mapping.csv")
            st.session_state["raw_csv"] = (outdir / "counts" / "raw_counts.csv").read_bytes()
            st.session_state["norm_csv"] = (outdir / "counts" / "normalized_counts.csv").read_bytes()
            st.session_state["de_csv"] = (outdir / "de" / "de_results.csv").read_bytes()
    except Exception as e:
        st.exception(e)
    finally:
        logging.getLogger("rnaseq_de").removeHandler(handler)


# ---------------- 结果展示 ----------------
if "de" in st.session_state:
    de = st.session_state["de"]
    sig = de[de["padj"] < alpha]
    m1, m2, m3 = st.columns(3)
    m1.metric(f"DE 基因（padj<{alpha}）", f"{len(sig)} / {len(de)}")
    m2.metric("上调", int((sig["log2FoldChange"] > 0).sum()))
    m3.metric("下调", int((sig["log2FoldChange"] < 0).sum()))

    tab1, tab2, tab3 = st.tabs(["📊 差异基因表", "🖼️ 图表", "💾 下载"])
    with tab1:
        st.dataframe(sig.sort_values("padj"), use_container_width=True,
                     height=400)
    with tab2:
        figs = st.session_state["figs"]
        cols = st.columns(2)
        order = ["volcano", "MA_plot", "PCA", "heatmap_top30"]
        for i, name in enumerate(order):
            if name in figs:
                cols[i % 2].image(figs[name], caption=name, use_container_width=True)
        with st.expander("比对率统计"):
            st.dataframe(st.session_state["qc"], use_container_width=True)
    with tab3:
        c1, c2, c3 = st.columns(3)
        c1.download_button("de_results.csv", st.session_state["de_csv"],
                           "de_results.csv", "text/csv")
        c2.download_button("raw_counts.csv", st.session_state["raw_csv"],
                           "raw_counts.csv", "text/csv")
        c3.download_button("normalized_counts.csv",
                           st.session_state["norm_csv"],
                           "normalized_counts.csv", "text/csv")
else:
    st.info("👈 在侧边栏点击「🚀 运行分析」开始；建议先勾选内置示例数据试跑。")
