/*
 * @Author: Darth_Eternalfaith darth_ef@hotmail.com
 * @Date: 2022-09-10 19:51:58
 * @LastEditors: Darth_Eternalfaith darth_ef@hotmail.com
 * @LastEditTime: 2022-09-10 19:54:28
 * @FilePath: \PrimitivesTGT-2D_Editor\js\import\PrimitivesTGT\test\test_matrix.js
 * @Description: test Matrix
 * 
 * Copyright (c) 2022 by Darth_Eternalfaith darth_ef@hotmail.com, All Rights Reserved. 
 */
import { Matrix, Matrix_2, Matrix_3, Vector } from "../NML";

// console.log("create_Print [1,2,3,4,5,6,7,8,9] ",Matrix.create_Print([1,2,3,4,5,6,7,8,9]));
// console.log("create_Identity",Matrix.create_Print(Matrix.create_Identity(3)));
// console.log("[1,2,3,4]*[5,6,7,8]",Matrix.create_Print(Matrix.multiplication([1,2,3,4],[5,6,7,8])));
// console.log("[1,2,3,4,5,6,7,8,9]*[9,8,7,6,5,4,3,2,1]",Matrix.create_Print(Matrix.multiplication([1,2,3,4,5,6,7,8,9],[9,8,7,6,5,4,3,2,1])));
// console.log("transposed [1,4,7,2,5,8,3,6,9]  ->  [1,2,3,4,5,6,7,8,9] ",Matrix.create_Print(Matrix.transpose([1,4,7,2,5,8,3,6,9])));
// console.log(Matrix.create_Print([4,1,1,1,1,4,1,1,1,1,4,1,1,1,1,4]));
// console.log("det [4,1,1,1,1,4,1,1,1,1,4,1,1,1,1,4]",Matrix.calc_Det([4,1,1,1,1,4,1,1,1,1,4,1,1,1,1,4]));
// console.log("det--transform [4,1,1,1,1,4,1,1,1,1,4,1,1,1,1,4]",Matrix.calc_Det__Transform([4,1,1,1,1,4,1,1,1,1,4,1,1,1,1,4]));
// console.log(Matrix.create_Print([4,1,1,1,0,4,1,1,1,1,0,1,1,1,1,0]));
// console.log("det [4,1,1,1,0,4,1,1,1,1,0,1,1,1,1,0]",Matrix.calc_Det([4,1,1,1,0,4,1,1,1,1,0,1,1,1,1,0]));
// console.log("det--transform [4,1,1,1,0,4,1,1,1,1,0,1,1,1,1,0]",Matrix.calc_Det__Transform([4,1,1,1,0,4,1,1,1,1,0,1,1,1,1,0]));
// console.log("create_Inverse [4,1,1,1,0,4,1,1,1,1,0,1,1,1,1,0] ",Matrix.create_Print(Matrix.create_Inverse([4,1,1,1,0,4,1,1,1,1,0,1,1,1,1,0])));
// console.log("create_Inverse [5,1,1,3,6,2,2,4,8] ",Matrix.create_Print(Matrix.create_Inverse([5,1,1,3,6,2,2,4,8])));


// console.log("create_Inverse [0, 7, 2, 2, 9, 9, 9, 8, 1] ",Matrix.create_Print(Matrix.create_Inverse([0, 7, 2, 2, 9, 9, 9, 8, 1])));
console.log("change size 2x2 Matrix [1,2,3,4] --> 3x3 Matrix ",Matrix.create_Print(Matrix.create_NewSize([1,2,3,4],2,4)));


function test_(){
    var test_matrix=[],
        temp,
        n;
    var i,j;
    console.log("创建测试矩阵数据 open");
    console.time();
    for(i=100;i>0;--i){
        // n=parseInt(Math.random()*5+3);
        n=3
        temp=[];
        for(j=n*n;j>0;--j){
            temp.push(parseInt(Math.random()*10));
        }
        test_matrix.push({
            data:temp,
            n:n
        })
    }
    // test_matrix=[{
    //     // data:[4,1,1,1,1,4,1,1,1,1,4,1,1,1,1,4],
    //     data:[4,1,1,1,0,4,1,1,1,1,0,1,1,1,1,0],
    //     n:4
    // }]
    console.timeEnd();
    console.log("创建测试矩阵数据 End",test_matrix);
    
    console.log("计算行列式 open");
    console.time();
        for(i=test_matrix.length-1;i>=0;--i){
            // console.log("det_c>"
            // +test_matrix[i].data.join(','),
            Matrix.calc_Det(test_matrix[i].data)
            // );
        }
    console.timeEnd();
    console.log("计算行列式 end");
    
    console.log("计算逆矩阵 open");
    console.time();
        for(i=test_matrix.length-1;i>=0;--i){
            // console.log("Inverse>"
            // +test_matrix[i].data.join(','),
            Matrix.create_Inverse(test_matrix[i].data)
            // );
        }
    console.timeEnd();
    console.log("计算逆矩阵 end");

    console.log("抽取样本 open");
        for(i=3;i>0;--i){
            // j=parseInt(Math.random()*100);
            j=0;
            temp=Matrix.create_Inverse(test_matrix[j].data);
            console.log(
                "原矩阵->",
                // Matrix.create_Print(
                    test_matrix[j].data,
                // ),
                "行列式->",
                Matrix.calc_Det(test_matrix[j].data),
                "逆矩阵->",
                // Matrix.create_Print(
                    temp,
                // ),
                "原矩阵右乘逆矩阵",
                // Matrix.create_Print(
                    Matrix.multiplication(test_matrix[j].data,temp),
                // )
            )
        }
    console.log("抽取样本 end");

    return test_matrix;
}

test_();

Object.assign(window,{
    Vector,
    Matrix,
    Matrix_2,
    Matrix_3,
});