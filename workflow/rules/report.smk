rule report_by_alg:
    message:
        "Take the results and provide a report"
    input:
        pred_rds = os.path.join(odir, "{algorithm}/prs.rds")
    output:
        report_rmd = os.path.join(odir, "{algorithm}/report.html")
    conda:
        "../envs/rmd.yaml"
    script:
        "../scripts/collect_results.Rmd"

rule prs_plot:
    input:
        pred_rds = os.path.join(odir, "{algorithm}/prs.rds"),
    output:
        prs_dist = report(os.path.join(odir, "prs_dist_{algorithm}.png"), 
            caption = "../report/distribution.rst", 
            category = lambda w: "{pheno}",
            subcategory = lambda w: "{algorithm}",
            labels={"trait": "{pheno}", "model": "{algorithm}"}
            )
    script:
        "../scripts/plot_prs_dist.R"
 