Modified from https://github.com/guogenglin/Bac_fetch

1. Download NCBI datasets at here: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
For macOS: 
Download datasets: curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
Download dataformat: curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat'
Make them executable: chmod +x datasets dataformat

2. Customize the output summary table. Parameters see here: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-assembly/?utm_source=chatgpt.com
  
2. Run python Bac_fetch.py
