# bacteria-phylogeny-procedure
The procedure of construct a phylogenetic tree



### 1.介绍

### 2.软件需求
基因组功能注释软件: prokka
基因组比对软件: clustalo
基因组模糊位点过滤: Gblocks
系统发育树构建软件: iqtree



### 3. 程序一站式直接运行

直接打包下载这些程序，然后直接运行`tree_build_orthofinder.py`程序即可。
```powershell
python tree_build_orthofinder.py -i genomes -a faa_dir -t 10
```
> -i: 基因组所在的文件夹
> -a: 注释的蛋白所在的文件夹
> -t: 线程数

### 4. 系统发育树构建流程(自己分步做)
##### 1.选择基因组
选择自己的基因组（`fasta`格式的）。
##### 2. 序列注释(得到蛋白序列) 
`prokka`是原核生物功能注释的很强大的工具，构建系统发育树之前需要对所有的基因组进行注释。以`GCA_000069185.1`为例

```powershell
prokka --outdir GCA_000069185.1 --prefix GCA_000069185.1 --noanno --cpus 8 --locustag 'GCA_000069185.1|ORF' GCA_000069185.1_genomic.fna
```

> --cpus：使用的线程数
>--prefix：文件结果前缀
>--locuatag： ORF的命名
>--outdir： 输出文件夹

##### 3. orthofinder 聚类(mcl)，寻找单拷贝基因
`orthofinder`采用mcl聚类的方法，对基因进行聚类，从而得到直系同源基因，以及用于构建物种树的单拷贝基因。

```powershell
orthofinder -og -f faa_dir -t 60
```

> 将上一步prokka注释的蛋白文件“`.faa`”结尾的文件放在`faa_dir`文件夹下，然后开始`orthofinder`。
> 最后在如下的目录中会生成很多单拷贝基因，用于构建系统发育树，文件夹对应的日期跟现实日期保持一致。
![在这里插入图片描述](https://img-blog.csdnimg.cn/656fd3253855476ba060f10723a3f976.png)
##### 4. 单拷贝基因的比对和对齐
在构建进化树之前，需要对单拷贝基因比对`clustalo`和对齐`alignment`，以及`Gblock`来过滤模糊位点（保留保守位点）。

```powershell
# 单拷贝基因的序列比对
# 首先需要新建文件夹`final_clustalo_fasta`
# 需要循环对所有的单拷贝基因做clustalo, 以该目录下的OG000001.fa为例子
clustalo -i faa_dir/OrthoFinder/Result_/Single_Copy_Orthologue_Sequences/OG000001.fa -o final_clustalo_fasta

#单拷贝基因组的序列对齐: 生成文件 `final.aln`
python AlignConcat.py -i final_clustalo_fasta -o final.aln

# 对齐序列的位点过滤: 生成文件`final.aln-gb`
Gblocks final.aln -t=p -e=-gb.fst -b4=5

```
##### 5. 使用`iqtree`对齐序列的构建系统进化树
使用对齐文件`final.aln-gb`进行构建系统发育树，模型参数按照需求更改。
```powershell
# 
iqtree -s final.aln-gb -m LG+F+R4 -bb 1000 -nt 45
```
