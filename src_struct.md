# 代码结构说明

* `writesac`, `readsac` 读取sac文件
* `wavePlot`, sac文件画图，效果未可知
* `suf`，返回文件的后缀名
* `siteInfo`，查找BHZ结尾的文件，获取经纬度和海拔。
* `seisShift`，不是移动时间，而是移动整体的数组
* `preEval`，对台站数据进行评估，基于`readPar`函数读取配置文件。

## RTW3D 调用思路

### 读取配置信息后的检查

norm0, norm1，用于归一化的时间范围。

至144行，参考`preEval`

meshgrid 网格，生成XI，YI

## 计算流程 1-3D_RBF

`Seisintep1D_RBF`, `Seisintep2D_RBF`, `Seisintep3D_RBF`

* `nodesLoc`：插值节点的位置
* `nodesValue`：插值节点的值
* `inteX`，X方向的网格， Xi 向量
* `rbf`， 进项接收函数的名称，如，gaussian, multiquadric,cubic,linear and thinplate
* `const`，高斯和多元二次径向接收函数的参数，通常是距离的均方差，默认为2
* `smooth`，平滑参数，默认为0
* `ZI`，插值得到的结果，1*N向量

处理方法：

1. lower， deblank，
2. 使用meshgrid，调整TI, XI, YI?
3. rbfcreate 创建rbf 插值
4. rbfIntep 重新插值

## cfg文件说明

1. norm，最小与最大的经纬度，插值step
2. 最小与最大的index of normalization window？需要与提供的数据量作比较。
3. DATAF -- 数据文件
   数据文件夹，lst文件的数量，名称，输出的目录
4. DIMtype
   插值的类型，'1D-RBF', '2D-RBF', '3D-RBF', '2D-CUBIC'、
    multiquadric, gaussian, cubic, linear, thinplate

## 台站数据评估`preEval`运算

`preEval`，函数

首先，接受cfg文件，cfg文件包含有一定的信息；

进行一些函数的预处理：

1. 归一化处理，根据所有记录中，振幅的最大值统一计算
2. 低通滤波？ 0.1 hz
3. 一些评估手段和补救措施

几个内部的参数：

srd, 搜索半径，0.5
cor, 互相关阈值，0.6
thN， thN ， 0.33 比例的 阈值

1. 搜索处于半径里的台站，并计算互相关
2. isgood判断好坏，半径里要有足够台站，互相关系数超过阈值的要占到一定比例
3. 如果较坏，需要简单做个互相关的到时校正（有方案一二，系统默认为方案二，即选一个互相关程度最好的台站。

## RBF计算


主要是`rbfcreate`和`rbfinterp`

### rbf_interp
 
几个参数：

* `x`， 插值点，维度与nodes相同？
* `phi`
* `rbfconst`
* `nodes`
* `rbfcoeff`

MATLAB size？ dimPoints， nPoints？

函数的输入参数包括插值点`x`和一个包含插值选项的结构体`options`。`options`结构体中包含了RBF插值所需的各种参数，如RBF函数`rbfphi`、RBF常数`RBFConstant`、插值节点`x`和RBF系数`rbfcoeff`。

函数首先检查插值点`x`的维度是否与插值节点`nodes`的维度相同。如果不同，函数会抛出一个错误。

然后，函数初始化了插值结果`f`和距离向量`r`。`f`是一个1行nPoints列的零矩阵，用于存储插值结果。`r`是一个1行n列的零矩阵，用于存储插值点到插值节点的距离。

接下来，函数进入一个循环，对每一个插值点进行插值。在每次循环中，函数首先计算插值点到插值节点的距离`r`，然后计算RBF插值的值`s`。RBF插值的值是由RBF函数、RBF系数和距离向量共同决定的。

在计算了RBF插值的值后，函数还会计算线性部分的值，并将其加到RBF插值的值上。最后，函数将计算得到的插值值存储到插值结果`f`中。

在所有插值点都进行了插值后，如果`options`结构体中的`Stats`字段为`'on'`，函数会打印出插值的统计信息，包括插值点的数量和插值的计算时间。

总的来说，这段代码实现了RBF插值的主要步骤，包括计算距离、计算RBF插值的值、计算线性部分的值和存储插值结果等。

## 辅助功能

* `cprintf`：C风格的printf
* `defval`：将string里的下划线进行替换，95为`_`的ASCII码
* `osdep`：检查文件系统类型，linux or solaris
