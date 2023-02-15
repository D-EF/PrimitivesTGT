/*!
 * @Author: Darth_Eternalfaith darth_ef@hotmail.com
 * @LastEditors: Darth_Eternalfaith darth_ef@hotmail.com
 * @LastEditTime: 2022-11-04 21:43:26
 * @FilePath: \site\js\import\PrimitivesTGT\Geometry_2D.js
 * @Description: 2D 几何图形
 * 
 * Copyright (c) 2022 by Darth_Eternalfaith darth_ef@hotmail.com, All Rights Reserved. 
 */

import { NML_CONFIG as CONFIG, Polygon, Vector } from "../NML/NML.js";

/** @typedef {import("../NML/NML").List_Value} Vec */


// open * 多边形 * open
    /** 多边形 */
    class Polygon{
        /** 
         * @param {List_Value} data 数据来源
         * @param {int} dimension   维 默认使用data的维度 或 二维(2D)
         */
        constructor (data,dimension){
            this.data=Array.from(data);
            var d=dimension||data.dimension||2;
            /** @type {int} 表示这个多边形在哪个维度 */
            this.dimension=d;
        }

        /** 使用向量数组创建
         * @param {List_Value[]} point_list 
         */
        static create_ByVectorList(point_list){
            var rtn=new Polygon([],point_list[0].length);
            
            return rtn;
        }

        /** 追加一个点
         * @param {List_Value}  p 点的数据
         * @param {Polygon}     p 点的数据
         * @return {}
         */
        static addPrint
    }
// end  * 多边形 * end 

// open * 矩形 * open
    /** 矩形 */
    class Rect extends CONFIG.VALUE_TYPE{
        /** 构造函数
         * @param {[x:number,y:number,width:number,height:number]} data 
         */
        constructor(data){
            super([
                data[0]||0,
                data[1]||0,
                data[2]||0,
                data[3]||0
            ])
        }

        /** 求最小坐标
         * @param {Rect}    rect 矩阵数据
         * @return {Vector} 返回一个2D向量
         */
        static get_Min(rect){
            return new Vector([
                rect[2]>=0?rect[0]:rect[0]+rect[2],
                rect[3]>=0?rect[1]:rect[1]+rect[3]
            ]);
        }

        /** 求最大坐标
         * @param {Rect}    rect 矩阵数据
         * @return {Vector} 返回一个2D向量
         */
        static get_Min(rect){
            return new Vector([
                rect[2]<=0?rect[0]:rect[0]+rect[2],
                rect[3]<=0?rect[1]:rect[1]+rect[3]
            ]);
        }

        /** 判断点是否在内部
         * @param {Rect} rect 矩阵数据
         * @param {Vec}  vec     点 2D向量
         * @return {Boolean} 返回 v 是否在矩形内部
         */
        static is_Inside(rect,vec){
            var x_min,x_max,y_min,y_max;
            if(rect[2]>=0){
                x_min = rect[0];
                x_max = rect[0]+rect[2];
            }else{
                x_max = rect[0];
                x_min = rect[0]+rect[2];
            }
            if(rect[3]>=0){
                y_min = rect[1];
                y_max = rect[1]+rect[3];
            }else{
                y_min = rect[1]+rect[3];
                y_max = rect[1];
            }

            return (
                vec[0]>x_min &&
                vec[0]<x_max &&
                vec[1]<y_max &&
                vec[1]>y_min
            );
        }

        /**
         * @return {Polygon}
         */
        static create_Polygon(){
            // todo
        }
    }
// end  * 矩形 * end 


export {
    Rect
}