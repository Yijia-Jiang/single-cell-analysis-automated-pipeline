configfile: "config.yaml"


# Check metasheet setup
from scripts.metasheet_setup import updateMeta
import pandas as pd
import yaml
import itertools

config = updateMeta(config)
print(f"ASSAY {config['assay']}")
print(f"BATCH {config['batch']}")
print(f"mergeRNA {config['mergeRNA']}")
print(f"mergeATAC {config['mergeATAC']}")
print(f"multiome {config['multiome']}")
metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#')

# Get comparsion columns for DEseq2
def _getColumn(comparison):
    return metadata["comp_{}".format(comparison)]

def _getComparison(name, group):
    comp = _getColumn(name)
    return metadata[comp == group].index

def _getSamples(wildcards):
    comp = _getColumn(wildcards.comparison)
    return comp.dropna().index


from string import Template

def getRuns(config):
    """parse metasheet for Run groupings"""
    ret = {}

    #LEN: Weird, but using pandas to handle the comments in the file
    #KEY: need skipinitialspace to make it fault tolerant to spaces!
    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
    f = metadata.to_csv().split() #make it resemble an actual file with lines
    #SKIP the hdr
    for l in f[1:]:
        tmp = l.strip().split(",")
        #print(tmp)
        ret[tmp[0]] = tmp[1:]

    #print(ret)
    config['runs'] = ret
    return config


config = getRuns(config)

#print(config["comparisons"])
# print(config)


#### scRNA #####

rule scrna_preprocess:
    input:
        lambda wildcards: config["rna_countsamples"][wildcards.sample]
    output:
        qc_out="analysis/scRNA/single_sample/{sample}/{sample}_QC.rds"
    shell:
        "Rscript scripts/1.1scRNA_preprocessing_quality_control_per_sample.R {input} {output.qc_out} {config[nCount_RNA_min]} {config[nCount_RNA_max]} {config[mito_rate_max]}"

rule scrna_doubletfinder:
    input:
        "analysis/scRNA/single_sample/{sample}/{sample}_QC.rds"
    output:
        doublet_filt="analysis/scRNA/single_sample/{sample}/{sample}_QC_doublet_removal.rds",
        doublet_df="analysis/scRNA/single_sample/{sample}/{sample}_QC_doublet.csv"
    shell:
        "Rscript scripts/1.2scRNA_preprocessing_doublet_detection_per_sample.R {input} {output.doublet_filt} {output.doublet_df}"


rule scrna_annotate:
    input:
        "analysis/scRNA/single_sample/{sample}/{sample}_QC_doublet_removal.rds"
    output:
        "analysis/scRNA/single_sample/{sample}/{sample}_DEG_annotated.rds"
    shell:
        "Rscript scripts/3.0DEG_annotation_per_sample.R {input} {output}"


rule scrna_singler:
    input:
        "analysis/scRNA/single_sample/{sample}/{sample}_QC_doublet_removal.rds"
    output:
        anno1="analysis/scRNA/single_sample/{sample}/{sample}_singler_annotation_hpca.csv",
        anno2="analysis/scRNA/single_sample/{sample}/{sample}_singler_annotation_encode.csv"
    shell:
        "Rscript scripts/3.1SingleR_annotation_per_sample.R {input} {output.anno1} {output.anno2}"

rule scrna_output_meta:
    input:
        sample_annot="analysis/scRNA/single_sample/{sample}/{sample}_DEG_annotated.rds",
        anno1="analysis/scRNA/single_sample/{sample}/{sample}_singler_annotation_hpca.csv",
        anno2="analysis/scRNA/single_sample/{sample}/{sample}_singler_annotation_encode.csv",
    output:
        metadata="analysis/scRNA/single_sample/{sample}/{sample}_final_metadata.csv",
        finaldata="analysis/scRNA/single_sample/{sample}/{sample}_final.rds"
    shell:
        "Rscript scripts/4.0scRNA_metadata.R {input.sample_annot} {input.anno1} {input.anno2} {output.metadata} {output.finaldata}"



rule scrna_deg_and_plot:
    input:
        finaldata="analysis/scRNA/single_sample/{sample}/{sample}_final.rds",
    output:
        p1="analysis/scRNA/single_sample/{sample}/{sample}_DEG_diffgene_per_cluster.csv",
        p2="analysis/scRNA/single_sample/{sample}/{sample}_DEG_diffgene_per_celltype_marker_annot.csv",
        p3="analysis/scRNA/single_sample/{sample}/{sample}_DEG_diffgene_per_celltype_singler_hpca.csv",
        p4="analysis/scRNA/single_sample/{sample}/{sample}_DEG_diffgene_per_celltype_singler_encode.csv"
    shell:
        "Rscript scripts/5.0scRNA_DEG_plotting_per_sample.R {input.finaldata} {output.p1} {output.p2} {output.p3} {output.p4}"

rule scrna_deg_and_GSEA:
    input:
        p1="analysis/scRNA/single_sample/{sample}/{sample}_DEG_diffgene_per_cluster.csv",
        p2="analysis/scRNA/single_sample/{sample}/{sample}_DEG_diffgene_per_celltype_marker_annot.csv",
        p3="analysis/scRNA/single_sample/{sample}/{sample}_DEG_diffgene_per_celltype_singler_hpca.csv",
        p4="analysis/scRNA/single_sample/{sample}/{sample}_DEG_diffgene_per_celltype_singler_encode.csv",
    output:
        out_dir=directory("analysis/scRNA/single_sample/{sample}/GSEA/"),
    shell:
        "Rscript scripts/3.1GSEA_hallmark.R {input.p1} {input.p2} {input.p3} {input.p4} {output.out_dir}"

rule scrna_barplot:
    input:
        finalmeta="analysis/scRNA/single_sample/{sample}/{sample}_final_metadata.csv",
    output:
        p1="analysis/scRNA/single_sample/{sample}/{sample}_prop_barplot_celltype_by_sample.pdf",
        # p2="analysis/scRNA/single_sample/{sample}/{sample}_prop_barplot_celltype_singler_hpca_by_sample.pdf",
        # p3="analysis/scRNA/single_sample/{sample}/{sample}_prop_barplot_celltype_singler_encode_by_sample.pdf",
        p2="analysis/scRNA/single_sample/{sample}/{sample}_prop_barplot_cellcycle_by_cluster.pdf",
        p3="analysis/scRNA/single_sample/{sample}/{sample}_prop_barplot_cellcycle_by_celltype.pdf",
    shell:
        "Rscript scripts/5.1scRNA_barplot.R {input.finalmeta} {output.p1} {output.p2} {output.p3}"

rule scrna_merge:
    input:
        samples=lambda wildcards: expand("analysis/scRNA/single_sample/{sample}/{sample}_QC_doublet_removal.rds", sample=config['mergeRNA'][wildcards.comparison]),
    output:
        # tmp="analysis/scRNA/integration/merge_{comparison}/",
        merge_samples="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_combine.rds",
        umap_p1="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_UMAP_by_cluster_before_integration.pdf",
        umap_p2="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_UMAP_by_sample_before_integration.pdf"
    params:
        #s1=lambda wildcards: ",".join(_getComparison(wildcards.comparison, 1)),
        #s1 = lambda wildcards: ",".join(config['mergeRNA'][wildcards.comparison][int(wildcards.group)])
        #LEN: returns ALL of the samples in the comparison, e.g.
        #for col mergeRNA_2, s1 = ['BI54_rna', 'MGH23039R_rna', 'BA-1_209_NormR_rna', 'BA-8_374_NormL_rna']
        s1 = lambda wildcards: ",".join(config['mergeRNA'][wildcards.comparison])
    shell:
        """Rscript scripts/2.0scRNA_merge.R "{input.samples}" "{params.s1}" "{output.merge_samples}" "{output.umap_p1}" "{output.umap_p2}"
        """


rule scrna_integrate_batch_correct:
    input:
        merge_samples="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_combine.rds"
    output:
        integrate_samples="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_combine_integrated.rds"
    params:
        annotFile=config['metasheet'],
    shell:
        "Rscript scripts/2.1scRNA_integrate.R \"{input.merge_samples}\" {params.annotFile} {output.integrate_samples}"


rule scrna_merge_annotation:
    input:
        merge_sample="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_combine_integrated.rds",
    output:
        merge_sample_annot="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_combine_integrated_annotated.rds"
    shell:
        "Rscript scripts/3.0DEG_annotation_integrate.R {input.merge_sample} {output.merge_sample_annot}"

rule scrna_merge_annotation_singler:
    input:
        merge_sample="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_combine_integrated.rds",
    output:
        anno1="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_singler_annotation_hpca.csv",
        anno2="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_singler_annotation_encode.csv"
    shell:
        "Rscript scripts/3.1SingleR_annotation_integrate.R {input.merge_sample} {output.anno1} {output.anno2}"

rule scrna_output_meta_merge:
    input:
        merge_sample_annot="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_combine_integrated_annotated.rds",
        anno1="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_singler_annotation_hpca.csv",
        anno2="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_singler_annotation_encode.csv",
    output:
        merge_metadata="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_final_metadata.csv",
        merge_finaldata="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_final.rds",
    shell:
        "Rscript scripts/4.0scRNA_metadata.R {input.merge_sample_annot} {input.anno1} {input.anno2} {output.merge_metadata} {output.merge_finaldata}"


rule scrna_custom_marker_annotate:
    input:
        sample_path="analysis/scRNA/single_sample/{sample}/{sample}_final.rds"
    output:
        "analysis/scRNA/single_sample/{sample}/{sample}_custom_marker_annotated.rds"
    params:
        annotFile=config['custom_marker_list'],
    shell:
        "Rscript scripts/6.2scRNA_annotation_input_markerlist.R {input.sample_path} {params.annotFile} {output}"


rule scrna_deg_plot_merge:
    input:
        merge_finaldata="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_final.rds",
    output:
        merge_deg_p1="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_DEG_diffgene_per_cluster.csv",
        merge_deg_p2="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_DEG_diffgene_per_celltype_marker_annot.csv",
        merge_deg_p3="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_DEG_diffgene_per_celltype_singler_hpca.csv",
        merge_deg_p4="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_DEG_diffgene_per_celltype_singler_encode.csv"
    shell:
        "Rscript scripts/5.0scRNA_DEG_plotting_integrate.R {input.merge_finaldata} {output.merge_deg_p1} {output.merge_deg_p2} {output.merge_deg_p3} {output.merge_deg_p4}"

rule scrna_barplot_merge:
    input:
        finalmeta="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_final_metadata.csv",
    output:
        p1="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_prop_barplot_celltype_by_sample.pdf",
        p2="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_prop_barplot_cellcycle_by_cluster.pdf",
        p3="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_prop_barplot_cellcycle_by_celltype.pdf",
        # p4="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_barplot_by_cell_cycle.pdf"
    shell:
        "Rscript scripts/5.1scRNA_barplot.R {input.finalmeta} {output.p1} {output.p2} {output.p3}"

rule scrna_deg_and_GSEA_merge:
    input:
        merge_deg_p1="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_DEG_diffgene_per_cluster.csv",
        merge_deg_p2="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_DEG_diffgene_per_celltype_marker_annot.csv",
        merge_deg_p3="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_DEG_diffgene_per_celltype_singler_hpca.csv",
        merge_deg_p4="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_DEG_diffgene_per_celltype_singler_encode.csv",
    output:
        merge_out_dir=directory("analysis/scRNA/integration/merge_{comparison}/GSEA/")
    shell:
        "Rscript scripts/3.1GSEA_hallmark.R {input.merge_deg_p1} {input.merge_deg_p2} {input.merge_deg_p3} {input.merge_deg_p4} {output.merge_out_dir}"


rule scrna_custom_marker_annotate_merge:
    input:
        sample_path="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_final.rds",
    output:
        "analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_custom_marker_annotated.rds"
    params:
        annotFile=config['custom_marker_list'],
    shell:
        "Rscript scripts/6.2scRNA_annotation_input_markerlist.R {input.sample_path} {params.annotFile} {output}"



rule scrna_merge_convert_mtx:
    input:
        merge_finaldata="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_final.rds",
    output:
        # h5ad_dir=directory("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad")
        h5ad_mtx="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/matrix.mtx",
        h5ad_barcode="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/barcodes.tsv",
        h5ad_feature="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/features.tsv",
        h5ad_metadata="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/metadata_final.csv",
    shell:
        "Rscript scripts/6.0scRNA_write_h5ad.R {input.merge_finaldata} {output.h5ad_mtx} {output.h5ad_barcode} {output.h5ad_feature} {output.h5ad_metadata}"


rule scrna_merge_convert_h5ad:
    input:
        h5ad_mtx="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/matrix.mtx",
        h5ad_barcode="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/barcodes.tsv",
        h5ad_feature="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/features.tsv",
        h5ad_metadata="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/metadata_final.csv",
    output:
        h5ad_data="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/merge_{comparison}.h5ad",
    shell:
        "python scripts/6.1write_h5ad_from_mtx.py --mtx {input.h5ad_mtx} --barcodes {input.h5ad_barcode} --features {input.h5ad_feature} --metadata {input.h5ad_metadata} --out {output.h5ad_data}"


rule scrna_stats:
    input:
        finalmeta="analysis/scRNA/single_sample/{sample}/{sample}_final_metadata.csv",
    output:
        df_stat="analysis/scRNA/single_sample/{sample}/{sample}_stats.csv",
    shell:
        "Rscript scripts/6.3scRNA_summary_per_sample.R {input.finalmeta} {output.df_stat}"


rule scrna_stats_cat:
    input:
        samples=lambda wildcards: expand("analysis/scRNA/single_sample/{sample}/{sample}_stats.csv", sample=config["rna_countsamples"])
    output:
        df_stat_cat="analysis/scRNA/all_sample_stats.csv"
    shell:
        """Rscript scripts/6.3scRNA_summary_all.R "{input.samples}" {output.df_stat_cat}
        """

rule rna_stats_merge:
    input:
        merge_metadata="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_final_metadata.csv",
    output:
        df_stat="analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_stats.csv",
    shell:
        "Rscript scripts/6.3scRNA_summary_merge.R {input.merge_metadata} {output.df_stat}"

# ####scATAC ##################

rule scatac_preprocess:
    input:
        samplepath=lambda wildcards: config["atac_countsamples"][wildcards.sample],
        fragpath=lambda wildcards: config["atac_fragments"][wildcards.sample],
        metapath=lambda wildcards: config["atac_metadata"][wildcards.sample],
    output:
        qc_out="analysis/scATAC/single_sample/{sample}/{sample}_QC.rds",
        df_output_dir="analysis/scATAC/single_sample/{sample}/{sample}_cell_filtered.csv",
        df_output_dir_barcodes="analysis/scATAC/single_sample/{sample}/{sample}_filtered_barcodes.csv"
    shell:
        "Rscript scripts/1.1scATAC_processing_quality_control_vSnakemake.R {input.samplepath} {input.fragpath} {input.metapath} {output.qc_out} {output.df_output_dir} {output.df_output_dir_barcodes}"


rule scatac_dim_reduction:
    input:
        qc_out="analysis/scATAC/single_sample/{sample}/{sample}_QC.rds"
    output:
        pro_out="analysis/scATAC/single_sample/{sample}/{sample}_QC_processed.rds"
    shell:
        "Rscript scripts/1.2scATAC_processing_dimension_reduction_v2.R {input.qc_out} {output.pro_out}"


def scatac_merge_inputFn(wildcards):
    samples = config['mergeATAC'][wildcards.comparison]
    #bedfiles = [f"analysis/scATAC/single_sample/{sample}/{sample}_peaks.bed" for sample in samples]
    bedfiles = [config['atac_peakbed'][sample] for sample in samples]
    fragfiles = [config['atac_fragments'][sample] for sample in samples]

    #fragfiles = [f"analysis/scATAC/single_sample/{sample}/{sample}_fragments.csv.gz" for sample in samples]
    metafiles = [f"analysis/scATAC/single_sample/{sample}/{sample}_cell_filtered.csv" for sample in samples]
    return {'bedfiles': bedfiles, 'fragfiles': fragfiles, 'metafiles':metafiles}

rule scatac_merge:
    input:
        # bedfiles=lambda wildcards: expand(config["atac_peakbed"][{sample}], sample=config['mergeATAC'][wildcards.comparison]),
        #fragfiles=lambda wildcards: expand(config["atac_fragments"][{sample}], sample=config['mergeATAC'][wildcards.comparison]),
        #metafiles=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/{sample}_cell_filtered.csv", sample=config["mergeATAC"][wildcards.comparison]),
        unpack(scatac_merge_inputFn)
    output:
	#merge_samples_path=directory("analysis/scATAC/integration/merge_{comparison}/"),
        merge_samples="analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_combine.rds",
        merge_peaks="analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_combined_peaks.rds",
        count_lists="analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_concat_count_lists.rds",
        umap_p1="analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_UMAP_by_sample_before_integration.pdf",
        # umap_p2="analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_UMAP_by_sample_before_integration.pdf",
    params:
        # s1 = lambda wildcards: ",".join(itertools.chain.from_iterable(config['mergeATAC'][wildcards.comparison]))
        s1 = lambda wildcards: ",".join(config['mergeATAC'][wildcards.comparison])
    shell:
        """Rscript scripts/2.0scATAC_merge.R "{input.bedfiles}" "{input.metafiles}" "{input.fragfiles}" "{params.s1}" "{output.merge_samples}" "{output.merge_peaks}" "{output.count_lists}" "{output.umap_p1}"
        """


rule scatac_integrate_batch_correct:
    input:
        merge_samples="analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_combine.rds"
    output:
        integrate_samples="analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_combine_integrated.rds",
        umap_p1="analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_UMAP_by_sample_after_integration.pdf",
    params:
        annotFile=config['metasheet'],
    shell:
        "Rscript scripts/2.1scATAC_integrate.R \"{input.merge_samples}\" {params.annotFile} {output.integrate_samples} {output.umap_p1}"


rule scatac_annotate_quickatac_step1_filterbarcode:
    input:
        fragpath=lambda wildcards: config["atac_fragments"][wildcards.sample],
        df_output_dir="analysis/scATAC/single_sample/{sample}/{sample}_filtered_barcodes.csv"
    output:
        frag_tmp1="analysis/scATAC/single_sample/{sample}/{sample}_cell_filtered_fragments.tsv"
    shell:
        "quick filter-barcodes {input.fragpath} --barcodes {input.df_output_dir} > {output.frag_tmp1}"

rule scatac_annotate_quickatac_step2_prepfiles:
    input:
        frag_tmp1="analysis/scATAC/single_sample/{sample}/{sample}_cell_filtered_fragments.tsv",
    params:
        spikefile=lambda wildcards: config["scatac_annotation_references"][wildcards.atlas]["spikein"],
        peak_bed=lambda wildcards: config["scatac_annotation_references"][wildcards.atlas]["peak_bed"],
    output:
        frag_tmp2="analysis/scATAC/single_sample/{sample}/{sample}_{atlas}_concat_fake_true_fragments.tmp.tsv",
        frag_tmp3="analysis/scATAC/single_sample/{sample}/{sample}_{atlas}_concat_fake_true_fragments.sorted.tsv"
    shell:
        "cat {params.spikefile} {input.frag_tmp1} > {output.frag_tmp2} && sort -k 1,1 -k2,2n {output.frag_tmp2} > {output.frag_tmp3}"


rule scatac_annotate_quickatac_step3_featurematrix:
    input:
        frag_tmp3="analysis/scATAC/single_sample/{sample}/{sample}_{atlas}_concat_fake_true_fragments.sorted.tsv",
    params:
        peak_bed=lambda wildcards: config["scatac_annotation_references"][wildcards.atlas]["peak_bed"],
        chromsize=config["chromsize_file"],
        count_dir= lambda wildcards: f"analysis/scATAC/single_sample/{wildcards.sample}/count_{wildcards.atlas}/{wildcards.sample}",
        # count_dir= lambda wildcards: f"analysis/scATAC/single_sample/{wildcards.sample}/count_{wildcards.atlas}/count_{wildcards.atlas}",
    output:
        tmp=directory("analysis/scATAC/single_sample/{sample}/count_{atlas, [A-Z]+}/"),
    shell:
        "mkdir -p {output.tmp} && quick agg-countmatrix {input.frag_tmp3} -g {params.chromsize} -p {params.peak_bed} --max-fragsize 999999999 -o {params.count_dir} "

rule scatac_annotate:
    input:
        count_dir="analysis/scATAC/single_sample/{sample}/count_{atlas}/",
    params:
        ref=lambda wildcards: config["scatac_annotation_references"][wildcards.atlas]["h5ad"],
        #chromsize=config["chromsize_file"]
        samplename=lambda wildcards: f"{wildcards.sample}",
        atlas=lambda wildcards: f"{wildcards.atlas}",
    output:
        #output_dir=directory("analysis/scATAC/single_sample/{sample}/count_{atlas}/scATAnno_{atlas}"),
        sample_h5ad="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}.h5ad",
        sample_metadata="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_query_annotated.csv",
    shell:
        "export OPENBLAS_NUM_THREADS=1 && python scripts/3.0scATAC_scATAnno_annotation.py --samplename {params.samplename} --input {input.count_dir} --atlas {params.atlas} --ref {params.ref} --out {output.sample_h5ad} --out_df {output.sample_metadata}"

rule scatac_output_meta:
    input:
        sample_annot="analysis/scATAC/single_sample/{sample}/{sample}_QC_processed.rds",
        anno1="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_query_annotated.csv",
    output:
        metadata="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_final_metadata.csv",
        finaldata="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_final.rds"
    shell:
        "Rscript scripts/4.0scATAC_metadata.R {input.sample_annot} {input.anno1} {output.metadata} {output.finaldata}"



rule scatac_annotate_multiple_samples:
    input:
        # count_dir=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/count_{wildcards.atlas}", sample=config["mergeATAC"][wildcards.comparison])
        count_dir=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/count_{atlas}", sample=config["mergeATAC"][wildcards.comparison], atlas=[wildcards.atlas])
    params:
        ref=lambda wildcards: config["scatac_annotation_references"][wildcards.atlas]["h5ad"],
        atlas=lambda wildcards: f"{wildcards.atlas}",
        samplename = lambda wildcards: " ".join(config['mergeATAC'][wildcards.comparison])
    output:
        merge_samples_h5ad="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}.h5ad",
        anno1="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/query_annotated.csv",
    shell:
        "export OPENBLAS_NUM_THREADS=1 && python scripts/3.0scATAC_scATAnno_annotation_multiple_samples.py --samplename {params.samplename} --input {input.count_dir} --atlas {params.atlas} --ref {params.ref} --out {output.merge_samples_h5ad} --out_df {output.anno1}"

rule scatac_output_meta_merge:
    input:
        sample_annot="analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_combine_integrated.rds",
        # merge_samples_h5ad="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}.h5ad",
        anno1="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/query_annotated.csv",
    output:
        metadata="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_final_metadata.csv",
        finaldata="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_final.rds"
    shell:
        "Rscript scripts/4.0scATAC_metadata_merge.R {input.sample_annot} {input.anno1} {output.metadata} {output.finaldata}"


rule scatac_barplot_merge:
    input:
        finalmeta="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_final_metadata.csv",
    output:
        p1="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_prop_barplot_celltype_by_sample.pdf",
    shell:
        "Rscript scripts/5.1scATAC_barplot.R {input.finalmeta} {output.p1}"


rule scatac_tf_analysis:
    input:
        finaldata="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_final.rds"
    output:
        da_peaks_dir="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/tf_analysis/{sample}_da_peaks_by_celltypes.rds",
        #motif_dir="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_enriched_motifs_by_celltypes.rds",
        #plot_motif="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/tf_analysis/{sample}_enriched_motifs.pdf",
    params:
        samplename=lambda wildcards: f"{wildcards.sample}",
    shell:
        "Rscript scripts/5.1scATAC_tf_analysis.R {input.finaldata} {params.samplename} {output.da_peaks_dir}"



rule scatac_barplot:
    input:
        finalmeta="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_final_metadata.csv",
    output:
        p1="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_prop_barplot_celltype_by_sample.pdf",
    shell:
        "Rscript scripts/5.1scATAC_barplot.R {input.finalmeta} {output.p1}"


rule scatac_stats:
    input:
        finalmeta="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_final_metadata.csv",
    output:
        df_stat="analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_stats.csv",
    params:
        atlas=lambda wildcards: f"{wildcards.atlas}",
    shell:
        "Rscript scripts/6.3scATAC_summary_per_sample.R {input.finalmeta} {output.df_stat}"


rule scatac_stats_cat:
    input:
        samples=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_stats.csv", sample=config["atac_countsamples"], atlas=config['scatac_annotation_references'])
    output:
        df_stat_cat="analysis/scATAC/all_sample_stats.csv"
    #params:
        #atlas=lambda wildcards: f"{wildcards.atlas}",
    shell:
        """Rscript scripts/6.3scATAC_summary_all.R "{input.samples}" {output.df_stat_cat}
        """

rule atac_stats_merge:
    input:
        # merge_metadata=lambda wildcards: expand("analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_final_metadata.csv", comparison=config["mergeATAC"], atlas=config['scatac_annotation_references'])
        merge_metadata="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_final_metadata.csv",
    output:
        df_stat="analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_stats.csv",
    shell:
        "Rscript scripts/6.3scATAC_summary_merge.R {input.merge_metadata} {output.df_stat}"


# #### multiomic intergration ##########
def getMultiomeRNA(comparison):
    """Use wildcard.comparisons to return samples in the multiome_{comparison}
    column and whose ASSAY == RNA"""
    #list of samples whose assay == RNA, which is first elm in config['assay']
    assay_RNA = config['assay'][0]
    #Get samples in the multiome_{comparison} that are also assay_RNA
    samples = [s for s in config['multiome'][comparison] if s in assay_RNA]
    #samples = [s for s in itertools.chain.from_iterable(config['multiome'][comparison]) if s in assay_RNA]
    #print(samples)
    return samples

# ## merge all samples on columns begin with "multiome_" and assay is RNA
rule scrna_multiomic_merge:
    input:
        #samples=lambda wildcards: expand("analysis/scRNA/single_sample/{sample}/{sample}_QC_doublet_removal.rds", sample=config["rna_countsamples"]),
        samples=lambda wildcards: expand("analysis/scRNA/single_sample/{sample}/{sample}_QC_doublet_removal.rds", sample=getMultiomeRNA(wildcards.comparison)),
    output:
        merge_samples="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine.rds",
        umap_p1="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_UMAP_by_cluster_before_integration.pdf",
        umap_p2="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_UMAP_by_sample_before_integration.pdf"
    params:
        #s1=lambda wildcards: ",".join(_getComparison(wildcards.comparison, 1)),
        s1=lambda wildcards: ",".join(getMultiomeRNA(wildcards.comparison))
    shell:
        """Rscript scripts/7.0multiome_scRNA_merge.R "{input.samples}" "{params.s1}" "{output.merge_samples}" "{output.umap_p1}" "{output.umap_p2}"
        """

rule scrna_multiomic_batch_correct:
    input:
        merge_samples="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine.rds",
    output:
        integrate_samples="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine_integrated.rds",
        umap_p1="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_UMAP_by_sample_after_integration.pdf",
        #merge_samples="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine.rds",
        #umap_p1="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_UMAP_by_cluster_before_integration.pdf",
        #umap_p2="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_UMAP_by_sample_before_integration.pdf"
    params:
        annotFile=config['metasheet'],
    shell:
        """Rscript scripts/7.1multiome_scRNA_integrate.R {input.merge_samples} {params.annotFile} {output.integrate_samples} {output.umap_p1}
        """

rule scrna_multiomic_annotation:
    input:
        integrate_samples="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine_integrated.rds",
    output:
        merge_sample_annot="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine_integrated_annotated.rds"
    shell:
        """Rscript scripts/7.3multiome_scRNA_DEG_annotation_integrate.R {input.integrate_samples} {output.merge_sample_annot}
        """


def getMultiomeATAC(comparison):
    """Use wildcard.comparisons to return samples in the multiome_{comparison}
    column and whose ASSAY == ATAC"""
    #list of samples whose assay == ATAC, which is second elm config['assay']
    assay_ATAC = config['assay'][1]
    #Get samples in the multiome_{comparison} that are also assay_ATAC
    samples = [s for s in config['multiome'][comparison] if s in assay_ATAC]
    #samples = [s for s in itertools.chain.from_iterable(config['multiome'][comparison]) if s in assay_ATAC]
    #print(samples)
    return samples

def multiome_scatac_merge_inputFn(wildcards):
    samples=getMultiomeATAC(wildcards.comparison)
    bedfiles = [config['atac_peakbed'][sample] for sample in samples]
    fragfiles = [config['atac_fragments'][sample] for sample in samples]
    metafiles = [f"analysis/scATAC/single_sample/{sample}/{sample}_cell_filtered.csv" for sample in samples]
    rfiles = [f"analysis/scATAC/single_sample/{sample}/{sample}_QC_processed.rds" for sample in samples]
    return {'bedfiles': bedfiles, 'fragfiles': fragfiles, 'metafiles':metafiles, 'rfiles':rfiles }


rule scatac_multiomic_merge:
    input:
        unpack(multiome_scatac_merge_inputFn)
    output:
        merge_samples="analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_combine.rds",
        #merge_peaks="analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_combined_peaks.rds",
        #count_lists="analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_concat_count_lists.rds",
        umap_p1="analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_UMAP_by_sample_before_integration.pdf",
    params:
        s1=lambda wildcards: ",".join(getMultiomeATAC(wildcards.comparison))
    shell:
        """Rscript scripts/7.0multiome_scATAC_merge.R "{input.bedfiles}" "{input.metafiles}" "{input.fragfiles}" "{params.s1}" "{output.merge_samples}" "{output.umap_p1}" "{input.rfiles}"
        """

rule scatac_multiomic_integrate_batch_correct:
    input:
        merge_samples="analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_combine.rds",
    output:
        merge_samples_int="analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_combine_integrated.rds",
        umap_p1="analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_UMAP_by_sample_after_integration.pdf",
    params:
        annotFile=config['metasheet'],
    shell:
        """Rscript scripts/7.1multiome_scATAC_integrate.R {input.merge_samples} {params.annotFile} {output.merge_samples_int} {output.umap_p1}
        """

## merge from output of above two rules based on on columns begin with "multiome_"

rule multiomic_merge:
    input:
        scRNA_samples="analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine_integrated_annotated.rds",
        scATAC_samples="analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_combine_integrated.rds",
    output:
        integrate_samples="analysis/scMultiome/integration/merge_{comparison}/merge_{comparison}_scrna_scatac_integrated.rds"
    #params:
        #annotFile=config['metasheet'],
    shell:
        """Rscript scripts/8.0multiome_scRNA_scATAC_integration.R "{input.scRNA_samples}" "{input.scATAC_samples}" {output.integrate_samples}"""

# TODO: perform this after stats rule for single and merged
rule summary_report:
    input:
        rna_df_stat_cat="analysis/scRNA/all_sample_stats.csv",
        rna_merge_df_stat=lambda wildcards: expand("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_stats.csv",comparison=config["mergeRNA"]),
        atac_df_stat_cat="analysis/scATAC/all_sample_stats.csv",
        atac_merge_df_stat=lambda wildcards: expand("analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_stats.csv",comparison=config["mergeATAC"], atlas=config['scatac_annotation_references']),
    params:
        #atlas=lambda wildcards: f"{wildcards.atlas}",
        single_path="analysis/scRNA/single_sample/",
        integrate_path="analysis/scRNA/integration/",
        single_path_atac="analysis/scATAC/single_sample/",
        integrate_path_atac="analysis/scATAC/integration/",
        integrate_path_scrna_scatac="analysis/scMultiome/integration",
        atlas=lambda wildcards: f"{wildcards.atlas}",
    output:
        report="analysis/summary_report_{atlas}.Rmd"
    shell:
        """Rscript scripts/9.0single-cell_report.R "{params.single_path}" "{params.integrate_path}" {params.atlas} "{params.single_path_atac}" "{params.integrate_path_atac}" "{params.integrate_path_scrna_scatac}" {output.report}"""

###### collect all output ############

rule all:
    input:
        scRNA_final_sample=lambda wildcards: expand("analysis/scRNA/single_sample/{sample}/{sample}_final.rds",sample=config["rna_countsamples"]),
        scRNA_metadata=lambda wildcards: expand("analysis/scRNA/single_sample/{sample}/{sample}_final_metadata.csv",sample=config["rna_countsamples"]),
        scRNA_deg_p1=lambda wildcards: expand("analysis/scRNA/single_sample/{sample}/{sample}_DEG_diffgene_per_celltype_marker_annot.csv",sample=config["rna_countsamples"]),
        scRNA_barplot_p1=lambda wildcards: expand("analysis/scRNA/single_sample/{sample}/{sample}_prop_barplot_celltype_by_sample.pdf",sample=config["rna_countsamples"]),
        scRNA_custom_anno=lambda wildcards: directory(expand("analysis/scRNA/single_sample/{sample}/{sample}_custom_marker_annotated.rds",sample=config["rna_countsamples"])),
        scRNA_gsea_out_dir=lambda wildcards: directory(expand("analysis/scRNA/single_sample/{sample}/GSEA/",sample=config["rna_countsamples"])),
        scRNA_merge_sample_final=lambda wildcards: expand("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_final.rds",comparison=config["mergeRNA"]),
        scRNA_merge_metadata=lambda wildcards: expand("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_final_metadata.csv",comparison=config["mergeRNA"]),
        scRNA_merge_deg_p1=lambda wildcards: expand("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_DEG_diffgene_per_cluster.csv",comparison=config["mergeRNA"]),
        scRNA_merge_barplot_p1=lambda wildcards: expand("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_prop_barplot_celltype_by_sample.pdf",comparison=config["mergeRNA"]),
        scRNA_merge_gsea_out_dir=lambda wildcards: directory(expand("analysis/scRNA/integration/merge_{comparison}/GSEA/",comparison=config["mergeRNA"])),
        #scRNA_merge_h5ad_mtx= lambda wildcards: expand("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/matrix.mtx", comparison=config["mergeRNA"]),
        #scRNA_merge_h5ad_h5ad= lambda wildcards: expand("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_h5ad/merge_{comparison}.h5ad", comparison=config["mergeRNA"]),
        scRNA_stats="analysis/scRNA/all_sample_stats.csv",
        scRNA_stats_merge=lambda wildcards: expand("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_stats.csv", comparison=config["mergeRNA"]),
        scRNA_custom_anno_merge=lambda wildcards: expand("analysis/scRNA/integration/merge_{comparison}/merge_{comparison}_custom_marker_annotated.rds",comparison=config["mergeRNA"]),
        
        ### scATAC output ##### merge output should based on columns begin with "mergeATAC_"
        scATAC_final_sample=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/{sample}_QC_processed.rds",sample=config["atac_countsamples"]),
        scATAC_quickatac_tmp=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/{sample}_cell_filtered_fragments.tsv", sample=config["atac_countsamples"]),
        scATAC_quickatac_tmp2=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/{sample}_{atlas}_concat_fake_true_fragments.sorted.tsv", sample=config["atac_countsamples"], atlas=config['scatac_annotation_references']),
        scATAC_annotation=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}.h5ad", sample=config["atac_countsamples"], atlas=config['scatac_annotation_references']),
        scATAC_annotation_df=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_query_annotated.csv", sample=config["atac_countsamples"], atlas=config['scatac_annotation_references']),
        scATAC_final=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_final.rds", sample=config["atac_countsamples"], atlas=config['scatac_annotation_references']),
        scATAC_motif=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/tf_analysis/{sample}_da_peaks_by_celltypes.rds", sample=config["atac_countsamples"], atlas=config['scatac_annotation_references']),
        scATAC_barplot_p1=lambda wildcards: expand("analysis/scATAC/single_sample/{sample}/scATAnno_{atlas}/{sample}_prop_barplot_celltype_by_sample.pdf",sample=config["atac_countsamples"], atlas=config['scatac_annotation_references']),
        merge_h5ad_final=lambda wildcards: expand("analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}.h5ad",comparison=config["mergeATAC"], atlas=config['scatac_annotation_references']),
        scATAC_merge_data=lambda wildcards: expand("analysis/scATAC/integration/merge_{comparison}/merge_{comparison}_combine_integrated.rds",comparison=config["mergeATAC"]),
        scATAC_merge_metadata=lambda wildcards: expand("analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_final_metadata.csv",comparison=config["mergeATAC"],atlas=config['scatac_annotation_references']),
        scATAC_merge_barplot_p1=lambda wildcards: expand("analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_prop_barplot_celltype_by_sample.pdf",comparison=config["mergeATAC"],atlas=config['scatac_annotation_references']),
        atac_df_stat_cat="analysis/scATAC/all_sample_stats.csv",
        atac_stats_merge=lambda wildcards: expand("analysis/scATAC/integration/merge_{comparison}/scATAnno_{atlas}/merge_{comparison}_stats.csv", comparison=config["mergeATAC"], atlas=config['scatac_annotation_references']),


        ### scMultiome output ### merge output should based on columns begin with "multiome_"

        scMultiome_merge_rna=lambda wildcards: expand("analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine.rds",comparison=config['multiome']),
        scMultiome_merge_rna_batch=lambda wildcards: expand("analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine_integrated.rds",comparison=config['multiome']),
        scMultiome_merge_rna_batch_annot=lambda wildcards: expand("analysis/scMultiome/scRNA/merge_{comparison}/merge_{comparison}_combine_integrated_annotated.rds",comparison=config['multiome']),
        scMultiome_merge_atac=lambda wildcards: expand("analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_combine.rds",comparison=config['multiome']),
        scMultiome_merge_atac_int=lambda wildcards: expand("analysis/scMultiome/scATAC/merge_{comparison}/merge_{comparison}_combine_integrated.rds",comparison=config['multiome']),
        scMultiome_merge_sample_final=lambda wildcards: expand("analysis/scMultiome/integration/merge_{comparison}/merge_{comparison}_scrna_scatac_integrated.rds",comparison=config['multiome']),
        report=lambda wildcards: expand("analysis/summary_report_{atlas}.Rmd",atlas=config['scatac_annotation_references']),
