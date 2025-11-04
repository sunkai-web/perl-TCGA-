#!/usr/bin/env perl

use strict;
use warnings;
use utf8;
use feature 'say';

# --- 核心模块 ---
# 1. JSON: 用于解析 metadata.json 文件
# 2. File::Spec: 用于跨平台安全地拼接文件路径
# 3. Cwd: 用于获取当前目录
use JSON;
use File::Spec;
use Cwd 'cwd';

# --- 配置 ---
my $metadata_filename = 'metadata.json'; # 你的 meta 文件名
my $data_dir_name     = 'file';          # 包含所有TSV文件的文件夹 (根据你上个问题)
my $output_filename   = 'fpkm_matrix.tsv'; # 最终生成的矩阵文件名

# --- 列配置 (基于GDC STAR-Counts文件格式) ---
# Perl 数组索引从 0 开始
my $GENE_SYMBOL_COL = 1; # 第2列 (gene_name)
my $FPKM_VALUE_COL  = 7; # 第8列 (fpkm_unstranded)
# --------------------

say "--- 开始处理 FPKM 矩阵生成 ---";

# === 步骤 1: 解析 metadata.json 文件 ===
say "1. 正在解析 $metadata_filename ...";

my %file_to_tcga; # 哈希: $file_to_tcga{文件名} = TCGA ID

eval {
    open(my $fh, '<:utf8', $metadata_filename)
        or die "无法打开 $metadata_filename: $!";
    
    # 一次性读入所有内容
    my $json_text = do { local $/; <$fh> };
    close $fh;

    # 解码 JSON
    my $json_data = decode_json($json_text);

    # 遍历JSON数组，构建映射
    foreach my $entry (@$json_data) {
        my $file_name = $entry->{'file_name'};
        my $tcga_id   = $entry->{'associated_entities'}->[0]->{'entity_submitter_id'};

        if ($file_name && $tcga_id) {
            $file_to_tcga{$file_name} = $tcga_id;
        } else {
            warn "警告: 在 metadata 中找到不完整的条目，已跳过。";
        }
    }
};
if ($@) {
    # 检查是否是 JSON 模块未安装
    if ($@ =~ /Can't locate JSON\.pm/) {
        die "致命错误: 'JSON' 模块未安装。\n请运行: cpanm JSON\n然后重试。\n";
    }
    # 其他错误
    die "解析 $metadata_filename 失败: $@";
}

my $total_files_in_meta = scalar(keys %file_to_tcga);
if ($total_files_in_meta == 0) {
    die "致命错误: 未能在 $metadata_filename 中找到任何 文件->TCGA ID 映射。";
}
say "   ...完成。在 metadata 中找到 $total_files_in_meta 个文件映射。";

# === 步骤 2: 遍历数据文件并构建矩阵 ===
say "2. 正在处理 'file' 文件夹中的 TSV 文件...";

my $base_dir      = cwd();
my $data_dir_path = File::Spec->catdir($base_dir, $data_dir_name);

my %data_matrix;        # 存储所有数据: $data_matrix{基因}{TCGA_ID} = FPKM
my %all_genes;          # 存储所有见过的基因名 (用于最后排序)
my @sample_order;       # 存储所有样本的TCGA ID (用于按顺序输出)

# 打开数据目录
opendir(my $dh, $data_dir_path)
    or die "无法打开数据文件夹 '$data_dir_path': $!";

# 遍历目录中的所有文件
while (my $file_name = readdir($dh)) {
    
    # 检查这个文件是否是我们在 metadata 中关心的文件
    unless (exists $file_to_tcga{$file_name}) {
        # say "   ...跳过文件 (不在 metadata 中): $file_name";
        next;
    }
    
    # 是我们关心的文件，获取其 TCGA ID
    my $tcga_id = $file_to_tcga{$file_name};
    push @sample_order, $tcga_id; # 记录样本顺序
    
    my $file_path = File::Spec->catfile($data_dir_path, $file_name);
    # say "   ...正在处理: $file_name (ID: $tcga_id)";

    open(my $fh, '<:utf8', $file_path)
        or do { warn "警告: 无法打开 $file_path: $!。已跳过。"; next; };

    # 逐行读取 TSV 文件
    while (my $line = <$fh>) {
        chomp $line;
        
        # 跳过注释行和统计行
        next if $line =~ /^(#|N_unmapped|N_multimapping|N_noFeature|N_ambiguous|gene_id)/;

        # 按 Tab 键分割
        my @cols = split(/\t/, $line);
        
        # 确保数据完整
        unless (@cols > $GENE_SYMBOL_COL && @cols > $FPKM_VALUE_COL) {
            # warn "警告: $file_name 中有格式不正确的行，已跳过。";
            next;
        }

        my $gene_symbol = $cols[$GENE_SYMBOL_COL];
        my $fpkm_value  = $cols[$FPKM_VALUE_COL];
        
        # 存储数据
        $data_matrix{$gene_symbol}{$tcga_id} = $fpkm_value;
        $all_genes{$gene_symbol} = 1; # 记录这个基因
    }
    close $fh;
}
closedir $dh;

say "   ...完成。处理了 " . scalar(@sample_order) . " 个样本文件。";

# === 步骤 3: 写入最终的矩阵文件 ===
say "3. 正在写入最终矩阵到 $output_filename ...";

# 对基因名和样本ID进行排序，确保输出文件顺序一致
my @sorted_genes   = sort keys %all_genes;
my @sorted_samples = sort @sample_order;

open(my $out_fh, '>:utf8', $output_filename)
    or die "无法创建输出文件 $output_filename: $!";

# 3.1 打印表头 (TCGA IDs)
# 先打印左上角的空单元格 (或基因列的标题)
print $out_fh "GeneSymbol";
# 循环打印所有 TCGA ID
foreach my $tcga_id (@sorted_samples) {
    print $out_fh "\t$tcga_id";
}
print $out_fh "\n"; # 表头行结束

# 3.2 打印数据行
foreach my $gene_symbol (@sorted_genes) {
    # 打印行头 (基因名)
    print $out_fh $gene_symbol;
    
    # 按照表头的顺序，依次打印该基因在每个样本中的值
    foreach my $tcga_id (@sorted_samples) {
        
        # 检查是否存在这个数据
        if (defined $data_matrix{$gene_symbol}{$tcga_id}) {
            my $value = $data_matrix{$gene_symbol}{$tcga_id};
            print $out_fh "\t$value";
        } else {
            # 如果某个样本中没有这个基因的数据 (理论上不应该，但作为安全措施)
            # 打印 "NA" (Not Available)
            print $out_fh "\tNA";
        }
    }
    print $out_fh "\n"; # 数据行结束
}

close $out_fh;

say "--- 全部完成！矩阵文件已保存为 $output_filename ---";