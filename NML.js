/*
 * @Author: Darth_Eternalfaith darth_ef@hotmail.com
 * @Date: 2022-08-26 01:29:45
 * @LastEditors: Darth_Eternalfaith darth_ef@hotmail.com
 * @LastEditTime: 2022-09-16 22:06:55
 * @FilePath: \PrimitivesTGT-2D_Editor\js\import\PrimitivesTGT\NML.js
 * @Description: Nittle Math Library
 * 
 * Copyright (c) 2022 by Darth_Eternalfaith darth_ef@hotmail.com, All Rights Reserved. 
 */

import { approximately, deg } from "../basics/math_ex.js";

/** @typedef {Float32Array} globalThis.NML_VALUE_TYPE */
globalThis.NML_VALUE_TYPE=globalThis.NML_VALUE_TYPE||Float32Array;

/** 配置  */
const CONGFIG={
    /** 向量使用的数据类型; 可选值为 Float32Array, Float64Array */
    VALUE_TYPE:globalThis.NML_VALUE_TYPE,
    /** @type {float} 计算容差 */
    APPROXIMATELY_TOLERANCE:1e-6
};

const {sin,cos,asin,acos,abs,sqrt}=Math;

// 类型注释 open
    /** @typedef {Float32Array} CONGFIG.VALUE_TYPE 矩阵计算时缓存下标的类型; 决定了计算时矩阵的n的大小 可选值为 Uint_N_Array, Int_N_Array */
    /** @typedef {Number} int      整形数字 */
    /** @typedef {Number} double   双浮点数字 */
    /** @typedef {Number} float    单浮点数字 */
    /** @typedef {Number[]|Float32Array|Float64Array|Vector} Vec 2D向量的各种存储形式 */
    /** @typedef {Number[]|Float32Array|Float64Array|Matrix} Mat 矩阵的各种存储形式 */
// 类型注释 end


/** 向量 */
class Vector extends CONGFIG.VALUE_TYPE{
    // 继承使用 CONGFIG.VALUE_TYPE 的构造函数

    /** 判断2d向量在哪个象限上, 规定 0 视作正
     * @param  {Vec} v 向量
     * @return {int} 返回在哪个象限
     */
    static v2__get_Quadrant(v){
        var f1=v[0]>=0,f2=v[1]>=0;
        if(f1){
            if(f2)  return 1;
            else    return 4;
        }else{
            if(f2)  return 2; 
            else    return 3; 
        }
    }
    
    /** 求模长
     * @param  {Vec} v 向量
     * @return {Number} 返回模长
	 */
	static mag(v) {
        var temp=0;
        for(var i =v.length-1;i>=0;--i){
            temp+=v[i]*v[i];
        }
		return sqrt(temp);
    }

    /** 判断某个向量是否为单位向量
     * @param {Vec} v 向量
     * @param {boolean} 返回是否为单位向量
     */
    is_Unit(){
        return abs(Vector.dot(v,v))<CONGFIG.APPROXIMATELY_TOLERANCE;
    }

    /** 创建标准化向量
     * @param  {Vec} v 向量
     * @return {Vector} 返回新的向量
	 */
    static create_Normalization(v){
        return Vector.normalize(Vector.copy(v));
    }

    /** 标准化向量
     * @param  {Vec} v 向量
     * @return {Vec} 修改并返回 v
	 */
	static normalize(v) {
        if(!Vector.is_Zero(v))throw new Error("This is a zero Vector.");
		var magSq = v.get_Mag(),oneOverMag=0;
		if (magSq>0) {
			oneOverMag = 1.0/magSq;
            for(var i =v.length-1;i>=0;--i){
                v[i] *= oneOverMag;
            }
		}
        return v;
    }

    /** 判断向量是不是零向量
     * @param  {Vec} v 向量
     * @return {Number} 返回0或非0
     */
    static is_Zero(v){
        var i=v.length;
        do{
            --i;
        }while((!v[i])&&i>0)
        return v[i];
    }
    
    /** 判断向量是否相等
     * @param  {Vec} v1 向量1
     * @param  {Vec} v2 向量2
     * @return {Boolean}
     */
    static is_Equal(v1,v2){
        var i=v1.length;
        if(i!==v2.length)return false;
        do{
            --i;
        }while((v1[i]===v2[i])&&i>0)
        return v1[i]===v2[i];
    }
    
    /** 取反
     * @param  {Vec} v 向量
     * @return {Vector} 返回新的向量
     */
    static instead(v){
        return Vector.instead_b(new Vector(v));
    }
    
    /** 取反
     * @param  {Vec} v 向量
     * @return {Vec} 修改并返回v
     */
    static instead_b(v){
        for(var i=v.length-1;i>=0;--i){
            v[i]*=-1;
        }
        return v;
    }

    /** 求向量和
     * @param  {Vec} v1 向量1
     * @param  {Vec} v2 向量2
     * @return {Vector} 返回新的向量
     */
    static sum(v1,v2){
        return Vector.translate(new Vector(v1),v2);
    }
    
    /** 再平移
     * @param {Vec} v1  原向量
     * @param {Vec} v2  偏移量向量
     * @return {Vec} 修改并返回 v1
     */
    static translate(v1,v2){
        if(v1.length!==v2.length) throw new Error("They vectors have different length!")
        for(var i=v1.length-1;i>=0;--i){
            v1[i]+=v2[i];
        }
        return v1;
    }
    
    /** 求向量差 1-2
     * @param {Vec} v1 向量1
     * @param {Vec} v2 向量2
     * @return {Vector} 返回一个新向量
     */
    static dif(v1,v2){
        // return Vector.translate(Vector.instead(v2),v1);
        
        if(v1.length!==v2.length) throw new Error("They vectors have different length!")
        var rtn=new Vector(v1);
        for(var i=rtn.length-1;i>=0;--i){
            rtn[i]-=v2[i];
        }
        return rtn;
    }
    
    /** 数字乘向量 
     * @param {Vec} v    向量
     * @param {Number} n 标量
     * @return {Vector} 返回新的向量
     */
     static np(v,n){
        return Vector.np_b(new Vector(v),n);
    }

    /** 数字乘向量 
     * @param {Vec} v    向量
     * @param {Number} n 标量
     * @return {Vector} 修改并返回 v
     */
    static np_b(v,n){
        for(var i=v.length-1;i>=0;--i){
            v[i]*=n;
        }
        return v;
    }

    /** 向量内积
     * @param {Vec} v1 向量1
     * @param {Vec} v2 向量2
     * @return {Number} 返回 v1 * v2
     */
    static dot(v1,v2){
        if(v1.length!==v2.length) throw new Error("They vectors have different length!")
        var rtn=0;
        for(var i=v1.length-1;i>=0;--i){
            rtn+=v1[i]*v2[i];
        }
        return rtn;
    }

    /** 向量外积 仅支持 3D 和 2D 向量
     * @param {Vec} v1 向量1
     * @param {Vec} v2 向量2
     * @return {Number|Vec} 返回 v1 x v2
     */
    static cross(v1,v2){
        if(v1.length===2&&v2.length===2)return v1[0]*v2[1]-v1[1]*v2[0];
        else if(v1.length===3&&v2.length===3) return  new Vector([
            v1[1]*v2[2]-v1[2]*v2[1],    // x : y1z2-z1y2
            v1[2]*v2[0]-v1[0]*v2[2],    // y : z1x2-x1z2
            v1[0]*v2[1]-v1[1]*v2[0]     // z : x1y2-y1x2
        ]);
        else throw new Error("This function only run in 2D or 3D Vector! ");
    }

    /** 计算向量夹角 ∠AOB 的 cos
     * @param {Vec} v1 表示角的一边的射线上 的 向量A
     * @param {Vec} v2 表示角的一边的射线上 的 向量B
     * @return {Number} 返回夹角的cos值
     */
    static cos_VV(v1,v2){
        return Vector.dot(v1,v2)/(Vector.mag(v1)*Vector.mag(v2));
    }
}


/* ! 变换矩阵    ----------------------------------------------------------------------------------------------------------------------------------------------------------    变换矩阵 */


/**
 * 矩阵的数据类型为1维线性表:
 * [1, 2, 3]
 * [4, 5, 6]  >>  [1,2,3,4,5,6,7,8,9]
 * [7, 8, 9]
 */

/** 矩阵  */
class Matrix extends CONGFIG.VALUE_TYPE{
    // 继承使用 CONGFIG.VALUE_TYPE 的构造函数

    /** 创建打印用的二维数组
     * @param {Mat} m 矩阵
     * @param {int} w 矩阵有多少列(宽度)
     * @return {Number[][]} 
     */
    static create_Print(m,w){
        var l=m.length,i,
            n=parseInt(w||parseInt(sqrt(l)));
        var rtn=[];
        i=0;
        do{
            rtn.push(m.slice(i,i+=n));
        }while(i<l);
        return rtn;
    }

    /** 校验矩阵是否为方阵
     * @param {Mat} m    矩阵
     * @param {int} [_n]    n阶矩阵
     * @return {int} 返回n
     * @throws {Error} 当 _n 和 m的长度 无法形成方阵时 将会报错
     */
    static check_Square(m,_n){
        var n=parseInt(_n||sqrt(m.length));
        if(n*n!==m.length) throw new Error("This is not a square matrix! It should be a (n*n)!");
        return n;
    }

    /** 矩阵升级
     * @param {Matrix} m    原矩阵
     * @param {int} low_w           原矩阵宽度
     * @param {int} new_w           新矩阵宽度
     * @param {int} [_low_h]        原矩阵高度 无输入时将使用 low_w
     * @param {int} [_new_h]        新矩阵高度 无输入时将使用 new_w
     * @param {int} [_shift_top]    旧矩阵拷贝到新矩阵时的上方偏移 默认为
     * @param {int} [_shift_left]   旧矩阵拷贝到新矩阵时的左侧偏移 默认使用 _shift_top
     * @return {Matrix} 返回一个新矩阵
     */
     static create_NewSize(m,low_w,new_w,_low_h,_new_h,_shift_top,_shift_left){
        var low_h=_low_h||low_w,new_h=_new_h||new_w,
            shift_top  = (_shift_top&&((new_w+_shift_top)%new_w))||0,
            shift_left = _shift_left===undefined?shift_top:((new_h+_shift_left)%new_h),
            l=new_w*new_h,
            temp_u,temp_v,
            i,u,v;
        var rtn=new Matrix(l);
        u=new_w-1;
        v=new_h-1;
        temp_u=u-shift_left;
        temp_v=v-shift_top;
        for(i=l-1;i>=0;--i){
            if(temp_u>=low_w||temp_v>=low_h){
                rtn[i]=(u===v)?1:0;
            }else{
                rtn[i]=m[temp_v*low_w+temp_u];
            }
            --u;
            --temp_u;
            if(temp_u<0)temp_u=new_w-1;
            if(u<0){
                u=new_w-1;
                --v;
                --temp_v;
                if(temp_v<0)temp_v=new_h-1;
            }
        }
        return rtn;
    }

    /** 计算张量积
     * @param {Mat} m1 矩阵1
     * @param {Mat} m2 矩阵2
     * @param {int} [_w1] 矩阵1的宽度 默认认为 m1 是列向量(w1=1)
     * @param {int} [_h1] 矩阵1的高度 默认认为 m1 是列向量(h1=m1.length)
     * @param {int} [_w2] 矩阵2的宽度 默认认为 m2 是行向量(w2=m2.length)
     * @param {int} [_h2] 矩阵2的高度 默认认为 m2 是行向量(h2=1)
     * @return {Matrix} 返回一个新的矩阵
     */
    static create_TensorProduct(m1,m2,_w1,_h1,_w2,_h2){
        var w1=_w1||1,
            h1=_h1||m1.length,
            w2=_w2||m2.length,
            h2=_h2||1,
            i=w1*h1;
        var rtn=new Array(i);
        for(--i;i>=0;--i){
            rtn[i]=Matrix.np(m2,m1[i]||0);
        }
        return Matrix.concat(rtn,w1,w2,h1,h2);
    }

    /** 合并矩阵
     * @param  {Mat[]} m_list 传入多个矩阵,矩阵应该拥有相同大小
     * @param  {int} w_l      m_list中一行放几个矩阵
     * @param  {int} w_m      m_list[i]的宽度
     * @param  {int} [_h_l]   m_list中一列放几个矩阵
     * @param  {int} [_h_m]   m_list[i]的高度
     * @return {Matrix} 返回一个新的矩阵
     ```javascript
        Matrix.create_Concat([[1,2,3,4], [5,6,7,8]], 2, 2);
        // [1,2]   [5,6] >> [1,2,5,6] >> [1,2,5,6,3,4,7,8]
        // [3,4] , [7,8]    [3,4,7,8]
     ```
     */
    static concat(m_list,w_l,w_m,_h_l,_h_m){
        var h_l=_h_l||Math.ceil(m_list.length/w_l),
            h_m=_h_m||Math.ceil(m_list[0].length/w_m),
            l_l=w_l*h_l,
            l_m=w_m*h_m,
            l=l_l*l_m,
            w=w_l*w_m,
            u_l,v_l,u,v,i,j,k;
        var rtn=new Matrix(l);
        k=l_l;
        for(v_l=h_l-1;v_l>=0;--v_l){
            for(u_l=w_l-1;u_l>=0;--u_l){
                --k;
                j=l_m;
                for(v=h_m-1;v>=0;--v){
                    i=(v_l*h_m+v)*w+w_m*(u_l+1);
                    if(m_list[k])
                    for(u=w_m-1;u>=0;--u){
                        --i;
                        --j;
                        rtn[i]=m_list[k][j];
                    }
                }        
            }
        }
        return rtn;
    }

    /** 矩阵乘标量
     * @param {Mat}     m   矩阵
     * @param {Number}  k   标量
     * @return {Matrix} 返回一个新的矩阵
     */
    static np(m,k){
        return Matrix.np_b(new Matrix(m),k);
    }

    /** 矩阵乘标量
     * @param {Mat}     m   矩阵
     * @param {Number}  k   标量
     * @return {Mat} 修改m并返回
     */
    static np_b(m,k){
        var i;
        for(i=m.length-1;i>=0;--i){
            m[i]*=k;
        }
        return m;
    }

    /** 使用 uv 获取 index 
     * @param {int} n 矩阵宽度 (列数)
     * @param {int} u 元素的 u 坐标 (第u列)
     * @param {int} v 元素的 v 坐标 (第v行)
     */
    static get_Index(n,u,v){
        return v*n+u;
    }

    /** 创建单位矩阵
     * @param {int}  n   n阶矩阵
     * @return {Matrix} 
     */
    static create_Identity(n){
        var l=n*n, sp=n+1, i=l-1;
        var rtn=new Matrix(l);
        do{
            rtn[i]=1.0;
            i-=sp;
        }while(i>=0);
        return rtn;
    }
    
    /** 初等变换 换行操作
     * @param {Mat|Mat[]} m 一个或多个矩阵
     * @param {int} n       n阶矩阵 用来表示一行的长度
     * @param {int} v1      矩阵v坐标1 (要对调的行下标1)
     * @param {int} v2      矩阵v坐标2 (要对调的行下标2)
     * @return {m} 修改并返回m
     */
    static transform_Exchange(m,n,v1,v2){
        var i,j,k,l,t;
        var f=ArrayBuffer.isView(m[0])||Array.isArray(m[0]);
        // 换行
        for(i=v1*n,j=v2*n,k=n; k>0; --k,++i,++j){
            if(f){
                for(l=m.length-1;l>=0;--l){
                    t=m[l][i];
                    m[l][i]=m[l][j];
                    m[l][j]=t;  
                }
            }else{
                t=m[i];
                m[i]=m[j];
                m[j]=t;
            }
        }
        return m;
    }

    /** 初等变换 某行乘标量
     * @param {Mat|Mat[]} m     矩阵
     * @param {int} n           n阶矩阵
     * @param {int} v           矩阵的v坐标(行下标)
     * @param {Number} value    乘法中的标量部分
     */
    static transform_multiplication(m,n,v,value){
        var i,j,l;
        var f=ArrayBuffer.isView(m[0])||Array.isArray(m[0]);
        // 换行
        if(f){
            for(i=v*n,k=n; k>0; --k,++i){
                for(l=m.length-1;l>=0;--l) m[l][j]*=value;  
            }
        }else{
            for(i=v*n,k=n; k>0; --k,++i) m[j]*=value;
        }
        return m;
    }
    
    /** 将矩阵某个为0的项 通过初等变换的换行操作, 变成非0
     * @param {Mat|Mat[]} m     一个或多个矩阵
     * @param {int} index       当前下标
     * @param {int} v           当前v坐标(行下标)
     * @param {int} spl         寻址步长,应为 ±n
     * @param {int} [_index_m]   传入多个矩阵时使用哪个矩阵的值 默认0
     * @return {m} 
     */
    static exchange_zero(m,index,v,spl,_index_m){
        if(!spl) return m;
        var _v=v,i;
        var f=spl>0?1:-1;
        var tm=(ArrayBuffer.isView(m[0])||Array.isArray(m[0]))?m[_index_m||0]:m;
        for(i=index;tm[i]!==undefined;i+=spl,_v+=f){
            if(tm[i]){
                if(_v===v)  return m;
                else        return Matrix.transform_Exchange(m,Math.abs(spl),_v,v);
            }
        }
        // 找不到可以替换的
        return m;
    }

    /** 矩阵乘法    如果不传入矩阵宽高信息将视为方阵
     * @param {Mat} m1 矩阵1
     * @param {Mat} m2 矩阵2
     * @param {int} [_h1]   左矩阵的行数(高度)
     * @param {int} [_w1h2] 左矩阵的列数(宽度) 与 右矩阵的行数(高度)
     * @param {int} [_w2]   右矩阵的列数(宽度)
     * @return {Matrix} 返回一个新的矩阵
     */
    static multiplication(m1,m2,_h1,_w1h2,_w2){
        var n=_h1||_w1h2||_w2;
        if(!n){
            n=Matrix.check_Square(m1);
            Matrix.check_Square(m2,n);
        }
        var h1=_h1||n,
            w1h2=_w1h2||n,
            w2=_w2||n;
        if(h1<=0&&w1h2<=0&&w2<=0) throw new Error ("This is a null matrix!");
        var l=w2*h1, _u, _v, u, v, index, i;
        var rtn=new Matrix(l);
        
        for(v=h1-1;v>=0;--v){
            for(u=w2-1;u>=0;--u){
                _u=v*w1h2;
                _v=u;
                index=v*w2+u;
                for(i=w1h2;i>0;--i,++_u,_v+=w2){
                    rtn[index]+=m1[_u]*m2[_v]
                }
            }
        }
        return rtn;
    }
    
    /** 矩阵转置
     * @param {Mat} m 矩阵
     * @param {int} [_w] 矩阵宽度(列数)
     * @param {int} [_h] 矩阵高度(行数)
     * @return {m} 修改m并返回
     */
    static transpose(m,_w,_h){
        var u, v, tu, tv, temp;
        var n=_w||Matrix.check_Square(m,_n);
        if(n===(_h||n)){  //方阵
            for(v=n-1; v>0; --v){
                for(u=v-1; u>=0; --u){
                    tu=v*n+u;
                    tv=u*n+v;
                    temp=m[tv];
                    m[tv]=m[tu];
                    m[tu]=temp;
                }
            }
        }else{
            tu=m.length
            temp=new Array(m.length);
            for(--tu,v=_w-1,u=_h-1;tu>=0;--tu,--u){
                if(u===-1){
                    u=_h-1;
                    --v;
                }
                temp[tu]=m[u*_w+v]
            }
            u=m.length-1;
            do{
                m[u]=temp[u];
                --u;
            }while(u>=0)
        }
        return m;
    }

    /** 创建矩阵的转置
     * @param {Mat} m 矩阵
     * @param {int} [_n] 矩阵为n阶矩阵
     * @return {m} 返回一个新的矩阵
     */
    static create_Transposed(m,_n){
        return Matrix.transpose(new Matrix(m),_n);
    }
    
    /** 计算矩阵行列式
     * @param {Mat} temp_m 矩阵
     * @param {int} [_n] 矩阵为n阶矩阵
     * @return {Number} 返回矩阵的行列式
     */
    static calc_Det(m,_n){
        if(m.length===1)return  m[0];

        if(m.length===4)return  m[0]*m[3]-m[1]*m[2];

        if(m.length===9)return  m[0] * (m[4]*m[8] - m[5]*m[7])+
                                m[1] * (m[5]*m[6] - m[3]*m[8])+
                                m[2] * (m[3]*m[7] - m[5]*m[6]);
        // if(m.length===16)return m[0];
        return Matrix.calc_Det__Transform(m,_n);
    }

    /** 计算矩阵行列式 --使用初等变换
     * @param {Mat} temp_m 矩阵
     * @param {int} [_n] 矩阵为n阶矩阵
     * @return {Number} 返回矩阵的行列式
     */
    static calc_Det__Transform(m,_n){
        var n, uv, uv_i, i, j, flag=1;
        var k, sp;
        var m__transform, temp_row;
        n=Matrix.check_Square(m,_n);
        m__transform=new Matrix(m)
        for(uv=n-1; uv>0; --uv){
            uv_i=uv*n+uv;
            if(!m__transform[uv_i]){
                // 换行
                Matrix.exchange_zero(m__transform,uv_i,uv,-n);
                if(!m__transform[uv_i])return 0;
                else flag*=-1;
            }
            temp_row=m__transform.slice(uv_i-uv,uv_i);
            // 单位化
            if(!(m__transform[uv_i]===1)){
                sp=1/(m__transform[uv_i]);
                i = uv;
                for(--i; i>=0; --i){
                    temp_row[i]*=sp;
                }
            }
            // 消元
            for(uv_i-=n;uv_i>=0;uv_i-=n){
                k=m__transform[uv_i];
                if(k===0)continue;
                for(i=uv,j=uv_i; i>=0; --i,--j){
                    m__transform[j]-=k*[temp_row[i]];
                }
            }
        }
        for(sp=1,j=m__transform.length-1; j>=0; j-=n+1){
            sp*=m__transform[j];
        }
        return sp*flag;
    }

    /** 变换得到矩阵逆
     * @param {Mat} m       矩阵
     * @param {int} [_n]    矩阵为n阶矩阵
     * @return {Matrix}     修改 m 并返回
     */
    static inverse(m,_n){
        var n,uv,uv_i,i,j,v,temp;
        var k,sp;
        var _m=[];

        n=Matrix.check_Square(m,_n);
        _m[0]=m;
        _m[1]=Matrix.create_Identity(n);
        
        for(uv=n-1; uv>=0; --uv){
            uv_i=uv*n+uv;
            if(!_m[0][uv_i]){
                // 换行
                Matrix.exchange_zero(_m,uv_i,uv,-n);
                if(!_m[0][uv_i]){
                    console.warn("This is a singular matrix!");
                    return m;
                }
            }
            k=uv_i-uv;
            // 单位化
            if(!(_m[0][uv_i]===1)){
                sp=1/(_m[0][uv_i]);
                for(i = n-1; i>=0; --i){
                    _m[0][k+i]*=sp;
                    _m[1][k+i]*=sp;
                }
            }
            // 消元
            for(v=n-1,uv_i=v*n+uv;uv_i>=0;--v,uv_i=v*n+uv){
                k=_m[0][uv_i];
                if((k===0)||(uv===v))continue;
                temp=n*uv;
                for(i=n-1,j=(v+1)*n-1; i>=0; --i,--j){
                    if(i<=uv)_m[0][j]-=k*_m[0][temp+i];
                    _m[1][j]-=k*_m[1][temp+i];
                }
            }
        }
        return _m[1];
    }

    /** 求矩阵的逆 (创建逆矩阵)
     * @param {Mat} m       矩阵
     * @param {int} [_n]    矩阵为n阶矩阵
     * @return {Matrix}     返回一个新的矩阵
     */
    static create_Inverse(m,_n){
        return Matrix.inverse(new Matrix(m),_n);
    }
}

/** 用于创建2d变换矩阵的静态类 */
class Matrix_2 extends Matrix{
    /** 创建旋转矩阵
     * @param {Number} theta 顺时针 旋转角弧度
     * @return {Matrix_2}
     */
     static create_Rotate(theta){
        var s=Math.sin(theta),
            c=Math.cos(theta);
        return new Matrix_2([c,s,-s,c]);
    }

    /** 创建旋转矩阵 使用向量
     * @param {Vector} _v 2d向量
     * @return {Matrix_2}
     */
    static create_Rotate__v(_v){
        var v=Vector.is_Unit(_v)?_v:Vector.create_Normalization(_v);
        return new Matrix_2([v[0],v[1],-v[1],v[0]]);
    }

    /** 创建缩放矩阵
     * @param {Number} x x 轴方向上的缩放系数
     * @param {Number} y y 轴方向上的缩放系数
     * @return {Matrix_2}
     */
    static create_Scale(x,y){
        return new Matrix_2([x,0,0,y]);
    }

    /** 创建镜像矩阵(对称)
     * @param {Number} x 对称轴的法向 x 坐标
     * @param {Number} y 对称轴的法向 y 坐标
     * @return {Matrix_2}
     */
    static create_Horizontal (x,y){
        return new Matrix_2(
            1-2*x*x ,   -2*x*y,
            -2*x*y  ,   1-2*y*y
        )
    }

    /** 创建切变矩阵
     * @param {Number} axis 方向轴 0:x 非零:y
     * @param {Number} k 切变系数
     * @return {Matrix_2}
     */
    static create_Shear(axis,k){
        if(axis){
            // y轴
            return new Matrix_2([1,0,k,1]);
        }
        else{
            // x轴
            return new Matrix_2([1,k,0,1]);
        }
    }

    /** 创建单位矩阵
     * @return {Matrix_2}
     */
    static create_Identity(){
        return new Matrix_2([1,0,0,1]);
    }

    /** 创建等比缩放&旋转矩阵 根据向量生成矩阵
     * @param {Vec} v2 2d向量
     * @return {Matrix_2} 返回一个矩阵
     */
    static create_ByVector(v2){
        return new Matrix_2([v2[0],v2[1],-1*v2[1],v2[0]]);
    }
}

Matrix_2.ROTATE_90=new Matrix_2([0,1,-1,0]);
Matrix_2.ROTATE_90_I=new Matrix_2([0,-1,1,0]);
Matrix_2.FLIP_HORIZONTAL=new Matrix_2([-1,0,0,1]);


/** 用于创建3D变换矩阵的静态类
 *  规定统一使用左手坐标系
 *            ^  +y
 *            |     7 +z
 *            |  /  
 * -----------+-----------> +x
 *         /  |   
 *      /     |   
 *            |   
 */
class Matrix_3 extends Matrix{
    /** 创建缩放矩阵
     * @param {flot} x x坐标中的缩放系数
     * @param {flot} y y坐标中的缩放系数
     * @param {flot} z z坐标中的缩放系数
     * @returns {Matrix_3} 返回一个矩阵
     */
    static create_Scale(x,y,z){
        return new Matrix_3([
            x,0,0,
            0,y,0,
            0,0,z
        ]);
    }
    /** 创建旋转矩阵
     * @param {Number} theta 旋转弧度
     * @param {int} axis 旋转中心轴  [z,y,x] 默认为 0
     * @return {Matrix_3} 返回一个矩阵
     */
    static create_Rotate(theta,axis){
        var s=sin(theta),
            c=cos(theta);
        Matrix.create_NewSize([ c, s,-s, c,],2,2,2,2,axis,axis);
    }

    /** 创建正交投影矩阵
     * @param {Vec} v 投影面的法线
     * @return {Matrix_3} 返回一个矩阵
     */
    static create_Projection__Orthographic(v){
        var v=Vector.is_Unit(_v)?_v:Vector.create_Normalization(_v);
        // todo
    }

    /** 创建旋转矩阵, 使用旋转向量
     * @param {Vec} _v  3d向量
     * @return {Matrix_3} 返回一个矩阵
     */
    static create_Rotate__v(_v){
        var v=Vector.is_Unit(_v)?_v:Vector.create_Normalization(_v);
        // todo
    }

    /** 创建旋转矩阵, 使用欧拉角
     * @param {Vec} value 欧拉角参数
     * @param {Vec} axis 欧拉角的轴向顺序 [x,y,z]
     */
    static create_Rotate__EulerAngles(value,axis){
        // todo
    }
}

Matrix_3.ROTATE_X_CCW_90DEG = new Matrix_3([1, 0, 0, 0, 0, 1, 0, -1, 0 ]);
Matrix_3.ROTATE_X_CW_90DEG  = new Matrix_3([1, 0, 0, 0, 0, 1, 0, -1, 0 ]);
Matrix_3.ROTATE_X_180DEG    = new Matrix_3([1, 0, 0, 0, -1, 0, 0, -0, -1 ]);
Matrix_3.ROTATE_Y_CCW_90DEG = new Matrix_3([0, 0, -1, 0, 1, 0, 1, 0, 0 ]);
Matrix_3.ROTATE_Y_CW_90DEG  = new Matrix_3([0, 0, -1, 0, 1, 0, 1, 0, 0 ]);
Matrix_3.ROTATE_Y_180DEG    = new Matrix_3([-1, 0, -0, 0, 1, 0, 0, 0, -1 ]);
Matrix_3.ROTATE_Z_CCW_90DEG = new Matrix_3([0, 1, 0, -1, 0, 0, 0, 0, 1 ]);
Matrix_3.ROTATE_Z_CW_90DEG  = new Matrix_3([0, 1, 0, -1, 0, 0, 0, 0, 1 ]);
Matrix_3.ROTATE_Z_180DEG    = new Matrix_3([-1, 0, 0, -0, -1, 0, 0, 0, 1 ]);


export{
    Vector,
    Matrix,
    Matrix_2,
    Matrix_3,
}