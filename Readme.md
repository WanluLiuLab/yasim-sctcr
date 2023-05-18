# `yasim` -- Yet Another SIMulator

**Markdown compatibility guide** This file is written in [Myst-flavored Markdown](https://myst-parser.readthedocs.io/), and may show errors on the default landing page of PYPI or Git Hostings. You can correctly preview it on generated Sphinx documentation or [Visual Studio Code](https://code.visualstudio.com) with [ExecutableBookProject.myst-highlight](https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight) plugin.

---

With the development of Third-Generation Sequencing (TGS) and related technologies, accurate quantification of transcripts in isoform level with precise detection of novel isoforms from Alternative Splicing (AS) events or relocation of Transposable Elements (TEs) had become possible. YASIM is the tool that simulates Next- or Third-Generation bulk RNA-Sequencing raw FASTQ reads with ground truth genome annotation and realistic gene expression profile (GEP). It can be used to benchmark tools that are claimed to be able to detect isoforms (e.g., [StringTie](https://ccb.jhu.edu/software/stringtie/)) or quantify reads on an isoform level (e.g., [featureCounts](https://subread.sourceforge.net/featureCounts.html)).

YASIM serves for different simulation purposes. For example, it can be used to simulate count matrix from reference genome annotation, or to simulate raw FASTQ reads from user-provided count matrix. When combined with other tools, user can also simulate reads from genome with Single Nucleotide Polymorphism (SNP), Insertions \& Deletions (InDels), Structural Variations (SVs) and other genomic variations.

YASIM does not have the ability to generate machine noises for each sequencer, and third-party DNA- or RNA-Seq simulators (Referred to as Low-Level Read generators, LLRGs) are needed to convert cDNA sequences to reads with machine errors and quality information [^qual]. This gives YASIM extreme flexibility over sequencer models. Till now, YASIM can simulate most Illumina NGS sequencers and most pacBio/ONT TGS sequencing platforms.

[^qual]: Except PacBio Sequel model.

YASIM is designed to be modularized, as some of the modules are general-purpose and can be used in other simulation tasks. Implemented in Python 3, YASIM follows Object-Oriented Programming (OOP) styles and can be easily extended. Theoretically, YASIM can run on any platform that supports Python3. However, most LLRG are POSIX-only (i.e., work on GNU/Linux, MacOS and friends). So it is recommended to deploy this tool inside major GNU/Linux distributions like Ubuntu, Debian, CentOS, Fedora, etc. Using YASIM on Microsoft Windows Subsystem of Linux (WSL), version 1 or 2, is **NOT** recommended -- It would lead to impaired performance and may cause other problems due to LLRG incompatibilities. Using YASIM on other platforms (e.g., Oracle Solaris) is neither tested nor recommended.

## Installation

### Using pre-built Library from PYPI

You need Python interpreter (CPython implementation) >= 3.8 (recommended 3.8) and latest [`pip`](https://pip.pypa.io/) to install this software from [PYPI](https://pypi.org). Command:

```shell
pip install yasim==1.0.0
```

You are recommended to use this application inside a virtual environment like [`venv`](https://docs.python.org/3/library/venv.html), [`virtualenv`](https://virtualenv.pypa.io), [`pipenv`](https://pipenv.pypa.io), [`conda`](https://conda.io) or [`poetry`](https://python-poetry.org).

### Build from Source

You need Python interpreter (CPython implementation) >= 3.7, latest PYPA [`build`](https://pypa-build.readthedocs.io) and latest [`setuptools`](https://setuptools.pypa.io/) to build this software. You are recommended to build the software in a virtual environment provided by [`virtualenv`](https://virtualenv.pypa.io), etc.

Build the simulator using:

```shell
python3 -m build
pip install dist/yasim-1.0.0.tar.gz
```

### Installation of Third-Party Programs

For TGS LLRGs:

- PBSIM, which simulates PacBio RS C1 and C2 chemistry, with CCS support.
  - Official site: [GitHub](https://github.com/yukiteruono/pbsim)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/pbsim) [Debian](https://packages.debian.org/stable/pbsim)
  - Y. Ono, K. Asai, and M. Hamada, "Pbsim: Pacbio reads simulator–toward accurate genome assembly.," _Bioinformatics (Oxford, England)_, vol. 29, pp. 119–121, 1 Jan. 2013, ISSN : 1367-4811. DOI: [10.1093/bioinformatics/bts649](https://doi.org/10.1093/bioinformatics/bts649)
- PBSIM2, which simulates PacBio RS II P4C2, P5C3 and P6C4 chemistry, CLR only; ONT R94, R95 and R103 chemistry.
  - Official site: [GitHub](https://github.com/yukiteruono/pbsim2)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/pbsim2)
  - Y. Ono, K. Asai, and M. Hamada, "PBSIM2: A simulator for long-read sequencers with a novel generative model of quality scores," _Bioinformatics (Oxford, England)_, vol. 37, no. 5, pp. 589–595, May 5, 2021, Number: 5, ISSN: 1367-4811. DOI: [10.1093/bioinformatics/btaa835](https://doi.org/10.1093/bioinformatics/btaa835)
- PBSIM3, which simulates PacBio RS II and Sequel model, with CCS support.
  - Official site: [GitHub](https://github.com/yukiteruono/pbsim3)
  - Y. Ono, M. Hamada, and K. Asai, "Pbsim3: A simulator for all types of pacbio and ont long reads.," _NAR genomics and bioinformatics_, vol. 4, lqac092, 4 Dec. 2022, ISSN: 2631-9268. DOI: [10.1093/nargab/lqac092](https://doi.org/10.1093/nargab/lqac092)
- BadRead, which simulates arbitrary PacBio and ONT models.
  - Official site: [GitHub](https://github.com/rrwick/Badread)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/badread)
  - R. Wick, "Badread: Simulation of error-prone long reads," Journal of Open Source Software, vol. 4, no. 36, p. 1316, Apr. 2019. DOI: [10.21105/joss.01316](https://doi.org/10.21105/joss.01316)

For NGS LLRGs:

- ART, which simulates Illumina GenomeAnalyzer I, GenomeAnalyzer II, HiSeq 1000, HiSeq 2000, HiSeq 2500, HiSeqX PCR free, HiSeqX TruSeq, MiniSeq TruSeq, MiSeq v3, NextSeq500 v2.
  - Official site: [NIEHS](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/art) [Debian](https://packages.debian.org/stable/art-nextgen-simulation-tools)
  - W. Huang, L. Li, J. R. Myers, and G. T. Marth, "Art: A next-generation sequencing read simulator.," _Bioinformatics (Oxford, England)_, vol. 28, pp. 593–594, 4 Feb. 2012, ISSN: 1367-4811. DOI: [10.1093/bioinformatics/btr708](https://doi.org/10.1093/bioinformatics/btr708)
- DWGSIM, which simulates arbitrary Illumina models.
  - Official site: [GitHub](https://github.com/nh13/DWGSIM)
  - Other installation sources: [Conda](https://anaconda.org/bioconda/dwgsim) [Debian](https://packages.debian.org/stable/dwgsim)
  - NO PUB

You may refer to LLRG tutorial for detailed guidance on utilization of these software.

## NEWS

- 2022-11-09: The `yasim.commonutils` and `yasim.bioutils` package would be separated into a new package (still under internal development) and `yasim` of later versions would depend on it.
