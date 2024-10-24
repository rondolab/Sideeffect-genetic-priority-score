# Side effect genetic priority score (SE-GPS)

<h2>Summary</h2>
<p>We created an <i>in-silico</i> side effect genetic priority score that can inform the likelihood of an adverse side effect across 13,114 genes and 466 phenotypes. This score is constructed as a weighted sum using a multivariable mixed-effect regression model of five constructed genetic features (described below) with drug side effects to obtain the effect sizes from the association of each feature. This was applied to 19,422 protein-coding genes and 470 phenotypes. We further incorporated the direction of genetic effect in a complementary score, the SE-GPS-DOE, to mimic the mechanism of the drug using predictions of loss-of-function (LOF) and gain-of-function (GOF) from LoGoFunc<sup>1</sup> for clinical variants and QTL estimates for GWAS phenotypes. Positive directional scores reflect the likelihood of an adverse side effect from target inhibition, and a negative SE-GPS-DOE reflects target activation.</p>

<h2>Genetic Evidence</h2>
    
<p>Using publicly available sources, we collected nine data sources from three types of genetic evidence (clinical variants, coding variants and GWAS phenotypes) and combined these into five features. We mapped these features to phecodeX terms and restricted these to phecode integer terms. Phecodes that mapped to phecode categories: Neoplasms and Neonatal and Pregnancy categories were excluded from the analysis. Each of these genetic features is described below.</p>

<h3>1. Clinical variants </h3>
<p>We collected clinical variant genetic evidence from three sources:</p>
                           
<p><b>OMIM:</b> Mendelian genes from the Online Mendelian Inheritance in Man (OMIM)<sup>2</sup> database.</p>
    
<p><b>HGMD:</b> Disease-causing and likely disease-causing genes for human inherited diseases from the human mutation database (HGMD)<sup>3</sup>.</p>

<p><b>EVA-ClinVar:</b> Clinically relevant genetic variants and diseases from ClinVar<sup>4</sup>. This evidence is obtained from the Open Target platform and was curated by European Variation Archive (EVA)<sup>5</sup>. ClinVar evidence was filtered on clinical significance: 'likely pathogenic' and 'pathogenic' and the confidence of the submission was filtered on: 'criteria provided, multiple submitters, no conflicts', 'reviewed by expert panel', and 'practice guideline'.</p>

<h3>2. Gene burden</h3>
<p>We used two datasets to combine gene burden genetic evidence from rare variant collapsing analyses:</p>
              
<p><b>a) Open Target GB:</b> Gene burden data from whole exome and whole genome sequencing data sourced from Open Targets<sup>6</sup>. We extracted traits labeled ‘ICD first occurrence’ and restricted to missense and pLOF variants and <i>P</i> < 4.3 x 10<sup>-7</sup>. </p>

<p><b>b) RAVAR GB:</b> Gene level disease associations from the Rare Variant Association Repository (RAVAR)<sup>7</sup> an open database that compiles rare variant associations obtained via a literature search using a minor allele frequency (MAF) of less than 0.02 and <i>P</i> < 1.0 x 10<sup>-4</sup>. </p>

 <h3>3. Single-Variant</h3>
<p>We used two different datasets to combine single variant evidence from rare variant analyses:</p>
              
<p><b>a) Genebass SV:</b> Exome-based single variant association statistics from the UK Biobank hosted by Genebass<sup>8</sup> which uses <i>P</i> < 1.0 x 10<sup>-7</sup>. </p>

<p><b>b) RAVAR SV:</b> Rare variant associations from the Rare Variant Association Repository (RAVAR)<sup>7</sup> using a minor allele frequency (MAF) of less than 0.02 and <i>P</i> < 1.0 x 10<sup>-6</sup>. </p>
 

<h3>4. eQTL phenotype</h3>
<p>Genes with a GWA phenotype driven by gene expression regulation through shared variants. These phenotypes were obtained from the Pan-UK Biobank<sup>9</sup> and eQTL summary statistics were obtained from the Genotype-Tissue Expression (GTEx, v8) Portal across 49 tissues<sup>10</sup>.</p>

<h3>5. Locus2gene:</h3>
<p>Genes with a GWA phenotype identified by the Locus2Gene machine learning model from the Open Targets Genetics Portal which prioritizes likely causal genes at each GWAS locus by integrating fine-mapped associations with functional genomic features<sup>11</sup>.</p>
 
<h2>Study</h2>

<p>For further details about the methods and analysis behind this study, please see our paper: Duffy, A et al. Development of a Genetic Priority Score to Predict Drug Side Effects Using Human Genetic Evidence. Submitted. </p>
   
<p>The Side effect Genetic Priority Score R shiny app is accessible at: https://rstudio-connect.hpc.mssm.edu/sideeffect-geneticpriorityscore/ </p>


<p></p>


<h2>References</h2>

<p>1. David, S. et al. Genome-wide prediction of pathogenic gain- and loss-of-function variants from ensemble learning of a diverse feature set. <i>Genome Med</i> <b>15</b> 103 (2023).</p>
<p>2.	Hamosh, A. et al. Online Mendelian Inheritance in Man (OMIM), a knowledgebase of human genes and genetic disorders. <i>Nucleic Acids Res</i> <b>30</b>, 52-5 (2002).</p>
<p>3.	Stenson, P.D. et al. The Human Gene Mutation Database (HGMD(®)): optimizing its use in a clinical diagnostic or research setting. <i>Hum Genet</i> <b>139</b>, 1197-1207 (2020).</p>
<p>4.	Landrum, M.J. et al. ClinVar: improving access to variant interpretations and supporting evidence. <i>Nucleic Acids Research</i> <b>46</b>, D1062-D1067 (2017).</p>
<p>5.	Cook, C.E. et al. The European Bioinformatics Institute in 2016: Data growth and integration. <i>Nucleic Acids Research</i> <b>44</b>, D20-D26 (2015).</p>
<p>6.	Ochoa, D. et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. <i>Nucleic Acids Res</i> <b>49</b>, D1302-D1310 (2021).</p>
<p>7. Cao, C. et al. RAVAR: a curated repository for rare variant-trait associations. <i>Nucleic Acids Res</i> <b>52</b>, D990-D997 (2024).</p>
<p>8.	Karczewski, K.J. et al. Systematic single-variant and gene-based association testing of thousands of phenotypes in 394,841 UK Biobank exomes. <i>Cell Genomics</i> <b>2</b>, 100168 (2022).</p>
<p>9. Pan-UKB team https://pan.ukbb.broadinstitute.org. (2020).</p>
<p>10. Aguet, F. et al. The GTEx Consortium atlas of genetic regulatory effects across human tissues. <i>Science</i> <b>369</b>, 1318-1330 (2020).</p>
<p>11. Mountjoy, E. et al. An open approach to systematically prioritize causal variants and genes at all published human GWAS trait-associated loci. <i>Nat Genet</i> <b>53</b>, 1527-1533 (2021).</p>
 
