//
//  main.cpp
//  AETG3
//
//  Created by 李刚 on 2018/12/20.
//  Copyright © 2018 李刚. All rights reserved.
//

#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
using namespace std;

const int M = 50; //迭代次数
int max_row = 0;
//每行都存放着在dim维下该行所表示的参数组合
vector<vector<int> > p;
vector<vector<int> > p_lt;
int num;//参数个数
int t;//要生成的维度，设置成可选的
int nrow; //总行数不会变，定义为全局

int lt_combs[10][10000] = {0}; //类似于uncover_target_combs的存在，只是行数为t
int lt_ncol[10]; //支持的最大维度为10 可以更改 lt_ncol[i]:代表lt_combs[i]的size大小

vector<vector<bool> > uncover_target_combs; //未被覆盖的组合 false代表覆盖了 true表没有

//初始化p[][]
void init_p(const int r,const int c,vector<vector<int> >&p){
    for (int i = 0; i < r; ++i) {
        p.resize(i+1);
        for (int j = 0; j < c; ++j) {
            p[i].push_back(j);
        }
    }
}
//找最大值 也即找哪一行覆盖的组合数最多
int p_max(vector<int> p_combs){
    int k = 0;
    int max = p_combs[0];
    for (int i = 1; i < p_combs.size(); ++i) {
        if (p_combs[i] > max) {
            k = i;
            max = p_combs[i];
        }
    }
    return k;
}

//组合数 m表示参数的个数，n表示需覆盖的维度
int cmn(int m,int n){
    if (m <= n) {
        return 1;
    }
    int k = 1,s = 1;
    for (int i = m - n + 1; i <= m; ++i) {
        k *= i;
    }
    for (int i = 2; i <= n; ++i) {
        s *= i;
    }
    return k / s;
}

//给我一个行数，我就知道哪两个参数
//rnum:行数 dim:需要的是几维的表 p[]:存储dim个参数的数组 num：参数个数
void row_num(int rnum,vector<int> &p,int num,int dim,int j){
    if(rnum == 0 || num <= dim || num <= 0 || dim <= 0)return;
    int k = cmn(num - 1, dim - 1);
    if(rnum >= k) {
        for (int i = j;i < t; ++i) {
            ++p[i];
        }
        row_num(rnum - k, p, num - 1, dim, j);
    }
    else {
        row_num(rnum, p, num - 1, dim - 1,j + 1);
    }
}

//init all_combs
void init_all_combs(vector<vector<bool> > &v,const vector<int> &vp){
    for (int i = 0; i < nrow; ++i) { //i代表第几行
        v.resize(i+1);
        int ncol = 1; //代表此行应该的总列数
        int r = t;
        while (r > 0) {
            int m = p[i][--r]; //m表示参数序号
            ncol *= vp[m]; //vp[m] 表示第m号参数有多少个参数
        }
        for (int k = 0; k < ncol; ++k) { //给第i行赋值
            v[i].push_back(true); //未被覆盖的参数用true表示
        }
    }
}

void invert(int i,vector<int> vp,int *b,int m_row){
    int m = i;
    for (int j = t - 1; j >= 0; --j) {
        int pi = p[m_row][j];
        b[j] = m % vp[pi];
        m /= vp[pi];
    }
}

void init_lt_combs(const vector<vector<bool>> &v,const vector<int> &vp,int m_row){
    //先对lt_combs初始化
    for (int i = 0; i < t; ++i) { //从lt_combs的每一行开始赋值
        lt_ncol[i] = 1; //ncol表示lt_combs这一行的总列数
        for (int j = 0; j < t - 1; ++j) {
            int c = p_lt[i][j];
            int pi = p[m_row][c];
            lt_ncol[i] *= vp[pi];
        }
        for (int j = 0; j < lt_ncol[i]; ++j) {
            lt_combs[i][j] = 0; //每个元素初始化为0 表示此组合出现的次数为0
        }
    }
    //再对m赋值
    for (int i = 0; i < v[m_row].size(); ++i) { //对max_row一行遍历
        if (v[m_row][i]) { //未被覆盖时才遍历
            int b[10]; //把v[m_row][i]转换成数组放在b[]里 支持最大的维度为10
            invert(i, vp, b,m_row);
            for (int j = 0; j < t; ++j) { //v的每个元素相当于t条t-1维 遍历m的t行往里填
                int a = 1,col = 0; //a是每个参数的量级
                for (int k = t - 2; k >= 0; --k) {
                    int pj = p_lt[j][k]; // pj代表着一个参数号
                    col += b[pj]*a; //用代数的简单方法将p[]中的参数号 转换成列号
                    a *= vp[p[m_row][pj]]; //量级在变化 主要c根据vp[]里参数取值个数的大小
                }
                ++lt_combs[j][col];
            }
        }
    }
}

//判断 未被覆盖的组合集 是否为空
bool uncover_isEmpty(const vector<vector<bool>> &v){
    for (int i = 0; i < v.size(); ++i) {
        for (int j = 0; j < v[i].size(); ++j) {
            if (v[i][j]) {
                return false; //存在一个组合未被覆盖 返回false
            }
        }
    }
    return true;
}

//求最大整个表里最大的未被覆盖数是多少
int freq_occur(){
    //首先认为出现最多的是第一个参数 de 第一个值
    int k = lt_combs[0][0]; //k表示当前最大的出现次数
    //筛选
    for (int i = 0; i < t; ++i) {
        for (int j = 0; j < lt_ncol[i]; ++j) {
            if (lt_combs[i][j] > k) {
                k = lt_combs[i][j];
            }
        }
    }
    return k;
}

//求哪行哪列在未被覆盖的组合中出现的次数最多
void freq_occur2(int &r,int &c){
    //首先认为出现最多的是第一个参数 de 第一个值
    r = 0;
    c = 0;
    int k = lt_combs[0][0]; //k表示当前最大的出现次数
    //筛选
    for (int i = 0; i < t; ++i) {
        for (int j = 0; j < lt_ncol[i]; ++j) {
            if (lt_combs[i][j] > k) {
                r = i ;
                c = j;
                k = lt_combs[i][j];
            }
        }
    }
}

//随机函数
int myrandom(int i){
    return std::rand() % i;
}

//前t - 1个参数已经复制，其他参数重排 a是t - 1维参数
void permute(int *a,int *r){
    std::srand(unsigned(std::time(0)));
    for (int i = 0,k = 0; k < num - t + 1; ++i) {
        int j = 0;
        for (; j < t - 1; ++j) {
            if (i == a[j]) { //说明第i个参数被赋值
                break;
            }
        }
        if (j == t - 1) {    //说明第i个参数未被赋值
            r[k] = i;
            ++k;
        }
    }
    // using built-in random generator:
    std::random_shuffle ( r, r + num - t + 1);
    
    // using myrandom:
    std::random_shuffle ( r, r + num - t + 1, myrandom);
    
}

//返回覆盖集合的数目 dim表示覆盖的维度
int uncover_num(const vector<vector<bool> > &v,const vector<int> &candidate,const vector<int> &vp){
    int m = 0;
    int pj = 0;
    for (int i = 0; i < nrow; ++i) { //遍历覆盖表v的每一行
        int a = 1,col = 0; //a是每个参数的量级
        for (int j = t - 1; j >= 0; --j) {
            pj = p[i][j]; // pj代表着一个参数号
            if (candidate[pj] == -1) { //如果此组参数中存在-1值 跳出本层循环
                break;
            }
            col += candidate[pj]*a; //用代数的简单方法将p[]中的参数号 转换成列号
            a *= vp[pj]; //量级在变化 主要c根据vp[]里参数取值个数的大小
        }
        if (v[i][col] && candidate[pj] != -1) { //如果那个元素未被覆盖 m加1
            ++m;
        }
    }
    return m;
}

//choose第i个参数的值
void c_pi_val(const vector<vector<bool> > &v,vector<vector<int>> &candidate,const vector<int> &vp){
    srand(unsigned(time(NULL)));
    int r[num - t + 1]; //剩下的num - t + 1个参数未被赋值
    int a[10]; //记录已被赋值的t-1个参数 最大也不会超过10
    int pi; //pi是参数
    int k,cl;
    
    //把下标已经赋值过的参数
    for (int i = 0,j = 0; i < num; ++i) {
        if (candidate[0][i] != -1) { //不等于-1说明此参数已被赋值
            a[j++] = i; //记录下来 然后j自增1
        }
    }
    
    //给剩下的参数重排
    permute(a, r);
    //再对剩下的赋值
    for (int m = 0; m < M; ++m) {
        vector<int> c_candidate(candidate[m]);
        for (int i = 0; i < num - t + 1; ++i) {
            pi = r[i];
            k = rand() % vp[pi];
            //if (m[pi][k]) {}
            candidate[m][pi] = k; //每次初始化为随机数
            cl = uncover_num(v, candidate[m], vp); //cl为赋值过后的未覆盖组合的数量
            for (int j = 0; j < vp[pi]; ++j) {
                c_candidate[pi] = j;
                int k1 = uncover_num(v, c_candidate, vp);
                if (k1 > cl) {
                    candidate[m][pi] = j;
                    cl = k1;
                }
            }
            c_candidate[pi] = candidate[m][pi];
        }
    }
}

//给测试集赋值
void assign(vector<vector<int> > &v,const vector<int> &p,const int num){
    int nrow = int(v.size());
    v.resize(nrow + 1);
    for (int i = 0; i < num; ++i) {
        v[nrow].push_back(p[i]);
    }
}

//移去candidate被覆盖的组合
void mov_combs(vector<vector<bool> > &v,const vector<int> &candidate,const vector<int > &vp){
    for (int i = 0; i < nrow; ++i) { //遍历覆盖表v的每一行
        int a = 1,col = 0; //a是每个参数的量级
        for (int j = t - 1; j >= 0; --j) {
            int pj = p[i][j]; // pj代表着一个参数号
            col += candidate[pj]*a; //用代数的简单方法将p[]中的参数号 转换成列号
            a *= vp[pj]; //量级在变化 主要c根据vp[]里参数取值个数的大小
        }
        if (v[i][col]) { //如果那个元素未被覆盖 将其变成被覆盖
            v[i][col] = false;
        }
    }
}

int main(int argc, const char * argv[]) {
    //重定向输入流
    freopen("/Users/ligang/code_study/test2.txt", "r", stdin);
    freopen("/Users/ligang/code_study/testcase1.txt", "w", stdout);
    
    num = 0; //参数个数初始化为0
    vector<int> vp; //参数个数取值列表
    string str; //保存文件的每一行
    while(cin.peek() != EOF){ //文件不为空时遍历
        getline(cin, str); //获得第一行 理应为[parameter]
        
        if(str.find("parameter") != str.npos){ //如果第一行是[parameter]
            cin >> t; //先读入维度
            getline(cin, str); //把这一行剩下的内容忽略掉
            while (cin.peek() == 'p') { //第一个字符为‘p’ 说明这是一个参数行
                ++num; //一行代表一组参数 参数个数加1
                vp.push_back(0) ; //开始时参数取值个数为0
                getline(cin, str); //把这一行内容存入str里面
                std::size_t pos = 0; //记录查询位置
                
                while((pos = str.find(',',pos+1)) != str.npos){
                    ++vp[num - 1] ; //找有多少个逗号 参数能取值的个数为逗号数加1
                }
                ++vp[num - 1] ; //最后还要加一次
            }
        }
        else if(str.find("constriant") != str.npos){
            getline(cin, str);
            break;
        }
    }
    fclose(stdin);
    //参数个数即相应行数 不支持只有一个参数 不支持生成1维覆盖表
    if (num < t || t <= 1) {
        cerr << "文件内容不符合要求！";
        fclose(stdout);
        return 0;
    }
    //t维覆盖表总行数
    nrow = cmn(num,t);
    
    //对nrow行dim列的数组赋初值 每行都是该行所表示的参数组合
    init_p(nrow,t,p);
    for(int i = 0; i < nrow; ++i){
        row_num(i, p[i], num, t, 0);
    }
    
    //init all_combs
    init_all_combs(uncover_target_combs, vp) ;
    
    init_p(t, t-1, p_lt);
    for(int i = 0; i < t; ++i){
        row_num(i, p_lt[i], t, t-1, 0);
    }
    
    //候选集为nrow条
    vector<vector<int> > candidate(M);
    int case_num = 1;
    
    //主要筛选和打印过程-----------------------------------------------
    while (!uncover_isEmpty(uncover_target_combs)) {
        int max_lt_combs = 0;
        for (int m = 0; m < nrow; ++m) {
            //uncover_combs的每行都对应一个lt_combs 对其初始化
            init_lt_combs(uncover_target_combs, vp, m);
            int k = freq_occur(); //给r c赋值
            if (k > max_lt_combs) { //此组未覆盖的数目与当前未覆盖最大数比较
                max_lt_combs = k;
                max_row = m; //记录行号
            }
        }
        //用最大的行号重新初始化lt_combs
        init_lt_combs(uncover_target_combs, vp, max_row);
        int r,c; //记录lt_combs中值最大的行号列号
        freq_occur2(r, c);
        int h = c;
        vector<int> a(num,-1); //a给赋t-1的值
        for (int i = t - 2; i >= 0; --i) {
            int j = p_lt[r][i];
            int pi = p[max_row][j];
            a[pi] = h % vp[pi];
            h /= vp[pi];
        }
        for (int m = 0; m < M; ++m) { //每行生成一个测试用例
            //赋值第i条candidate
            candidate[m].assign(a.begin(), a.end());
        }
        c_pi_val(uncover_target_combs, candidate, vp) ;
        //choose candidate[M]并给testcase赋值--------------------------------------
        int k = 0;
        int max_uncover = 0;
        for (int i = 0; i < M; ++i) {
            int m = uncover_num(uncover_target_combs, candidate[i], vp) ;
            if (m > max_uncover) {
                k = i;
                max_uncover = m;
            }
        }
        //输出测试用例-----------------------------------
        std::cout << "测试用例" << case_num << ": " ;
        for (int j = 0; j < num; ++j) {
            std::cout << candidate[k][j];
        }
        std::cout << "\t\t" <<"此测试用例覆盖了" << max_uncover << "组" << std::endl;
        ++case_num; //测试用例数目加1
        
        //移除被覆盖的组合--------------------------------
        mov_combs(uncover_target_combs, candidate[k], vp) ;
    }
    //关闭输出流-----------------------------------------
    fclose(stdout);
    return 0;
}
