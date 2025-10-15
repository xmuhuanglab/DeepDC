![logo](imag/DeepDC logo_20250319.png)

DeepDC, developed by Jialiang Huang's Lab at Xiamen University (XMU) and Teng Fei's Lab at Northeastern University (NEU), is a deep learning-based computational model designed to predict dual-editing efficiency in non-coding regions. It combines a HyenaDNA block and XGBoost to extract epigenomic features and sequence information from each paired sgRNA and predicts editing efficiency.

![workflow](imag/workflow.png)

We also provide a [webserver](https://deepdc.huanglabxmu.com/) for user to design pgRNA and scoring.

# Run the Source Code:

## Dependencies

## Install


## Run the demo datasets
Here, we provide training and test demo input to showcase the DeepDC workflow:
| Chromosome | Start | End |
|:-----------:|:------:|:---:|
| chr1        | 10000  | 10500 |
| chr7        | 56000  | 56500 |
| chr12       | 91000  | 91500 |

### Preparation of input (Feature processing)
This Jupyter notebook is used to collect information related to the target fragment.: 
1. pgRNA: Includes paired gRNA sequences and PAM information.
2. Epigenetic Features: Examples include DNase, ATAC, H3K4me3, and H3K27ac.
3. Deleted fragment length and GC Content: Distance between paired sgRNAs is limited to 50–200 bp.
4. gRNA Efficiency: Predicted by DeepCpf1, DeepCRISPR, CRISPRedit, and Ruleset2.
5. Fragment Score: Calculated using the HyenaDNA scoring system.

### For training and scoring:
```bash
Training.ipynb Prediction.ipynb
```

# Contact us
If you have any further questions or encounter any issues, please feel free to contact us:
- [Liquan Lin: 21620241153548@stu.xmu.edu.cn](mailto:21620241153548@stu.xmu.edu.cn)
- [Shijie Luo: 24520230157443@xmu.edu.cn](mailto:24520230157443@xmu.edu.cn)
