<!--
 * @Author: Darth_Eternalfaith darth_ef@hotmail.com
 * @Date: 2022-09-10 20:58:17
 * @LastEditors: Darth_Eternalfaith darth_ef@hotmail.com
 * @LastEditTime: 2022-09-11 00:08:36
 * @FilePath: \PrimitivesTGT-2D_Editor\js\import\PrimitivesTGT\readme.md
 * @Description: 
 * 
 * Copyright (c) 2022 by Darth_Eternalfaith darth_ef@hotmail.com, All Rights Reserved. 
-->
# PrimitivesTGT 图元对象
*正在从 PrimitvesTGT-2D 重构中*   
用于图形渲染和其他计算的一个库，使用面向对象的手法编写   
使用 jsDoc 格式注释

---
## 这个库依赖于 Darth_Eternalfaith 的另一个库: <a href="https://gitee.com/darth_ef/basics">basics</a> 并且需要在同一级目录下,目录结构如下:   
...   
 * import   
    * basics
    * PrimitvesTGT


----
## 所有暴露可用的类、函数 ↓↓↓
###  一点点数学 Nittle Math Library: <a href="./NML.js.md"> NML.js </a>

---
---
---
# 正则表达式查找替换: 
* 静态成员函数的函数注释 --> md文档
    *  ` */\*\*(.*)([\s\S\n]*?)\*\/[\s\S\n]*?(\w.*)\{[\s\S\n]*?\/\*`
    *  `\n### $3  $1$2---\n/*`