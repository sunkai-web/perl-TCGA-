# perl-TCGA-
perl脚本整理TCGA转录组
#使用方法
1.文件配置
project_dir/
│
├── metadata.json          # GDC 下载时附带的元数据文件
├── file/                  # 存放所有样本的TSV表达文件
│    ├── a1f3...STAR.tsv
│    ├── b7e4...STAR.tsv
│    └── ...
└── make_matrix.pl         # 就是这个Perl标本
