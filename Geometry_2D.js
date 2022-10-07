/*
 * @Author: Darth_Eternalfaith darth_ef@hotmail.com
 * @LastEditors: Darth_Eternalfaith darth_ef@hotmail.com
 * @LastEditTime: 2022-10-07 16:57:52
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
         * @param {List_Value[]} list 
         */
        static create_ByVectorList(list){
            var rtn=new Polygon([],list[0].length);
            
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
         * @param {[x:Number,y:Number,width:Number,height:Number]} data 
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
         * @param {Rect}    r 矩阵数据
         * @return {Vector} 返回一个2D向量
         */
        static get_Min(r){
            return new Vector([
                r[2]>=0?r[0]:r[0]+r[2],
                r[3]>=0?r[1]:r[1]+r[3]
            ]);
        }

        /** 求最大坐标
         * @param {Rect}    r 矩阵数据
         * @return {Vector} 返回一个2D向量
         */
        static get_Min(r){
            return new Vector([
                r[2]<=0?r[0]:r[0]+r[2],
                r[3]<=0?r[1]:r[1]+r[3]
            ]);
        }

        /** 判断点是否在内部
         * @param {Rect} r 矩阵数据
         * @param {Vec}  v     点 2D向量
         * @returns {Boolean} 返回 v 是否在矩形内部
         */
        static is_Inside(r,v){
            var x_min,x_max,y_min,y_max;
            if(r[2]>=0){
                x_min = r[0];
                x_max = r[0]+r[2];
            }else{
                x_max = r[0];
                x_min = r[0]+r[2];
            }
            if(r[3]>=0){
                y_min = r[1];
                y_max = r[1]+r[3];
            }else{
                y_min = r[1]+r[3];
                y_max = r[1];
            }

            return (
                v[0]>x_min &&
                v[0]<x_max &&
                v[1]<y_max &&
                v[1]>y_min
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