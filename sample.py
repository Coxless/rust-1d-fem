#1次元Poisson方程式を、有限要素法で解く
#d/dx[p(x)du(x)/dx]=f(x) (x_min<x<x_max)
#u(x_min)=alpha, du(x_max)/dx=beta
import time  #時刻を扱うライブラリ
import numpy as np  #NumPyライブラリ
import scipy.linalg  #SciPyの線形計算ライブラリ
import matplotlib.pyplot as plt  #データ可視化ライブラリ



x_min = -1.0 #計算領域のXの最小値
x_max = 1.0 #計算領域のXの最大値

node_total = 4 #節点数(>=2)
ele_total = node_total-1 #要素数

func_f = 1.0 #定数関数f

#境界条件u(x_min)=alpha,du(x_max)/dx=beta。
#alphaやbetaが"inf"の時は、何も処理しないようにする。
alpha = 1.0 #左端のディリクレ境界条件
beta = -1.0 #右端のノイマン境界条件



#任意の節点に境界条件を実装する
def boundary(node_num_glo, Dirichlet, Neumann):
    #ディリクレ境界条件
    if (Dirichlet!="inf"):  #Dirichlet=無限大の時は処理しない
        vec_b_glo[:] -= Dirichlet*mat_A_glo[node_num_glo,:]  #定数ベクトルに行の値を移項
        vec_b_glo[node_num_glo] = Dirichlet  #関数を任意の値で固定
        mat_A_glo[node_num_glo,:] = 0.0  #行を全て0にする
        mat_A_glo[:,node_num_glo] = 0.0  #列を全て0にする
        mat_A_glo[node_num_glo,node_num_glo] = 1.0  #対角成分は1にする

    #ノイマン境界条件
    if (Neumann!="inf"):  #Neumann=無限大の時は処理しない
        vec_b_glo[node_num_glo] += Neumann #関数を任意の傾きで固定



############### 入力データの用意 ###############
#計算の開始時刻を記録
compute_time = time.time()

print("node_total = ",node_total, ",  ele_total = ",ele_total)

#Global節点のx座標を定義(x_min〜x_max)
print("Global node X")
node_x_glo = np.empty(node_total, np.float64) #Global節点のx座標
node_x_glo = np.linspace(x_min,x_max, node_total) #計算領域を等分割
print(node_x_glo)

#各要素のGlobal節点番号
print("Node number of each elements")
node_num_glo_in_seg_ele = np.empty((ele_total,2), np.int64) #各要素のGlobal節点番号
for e in range(ele_total):
    node_num_glo_in_seg_ele[e,0] = e
    node_num_glo_in_seg_ele[e,1] = e+1
print(node_num_glo_in_seg_ele)

#各要素のLocal節点座標
print("Local node X")
node_x_ele = np.empty((ele_total,2), np.float64) #各要素のLocal節点のx座標
for e in range(ele_total):
    for i in range(2):
        node_x_ele[e,i] = node_x_glo[ node_num_glo_in_seg_ele[e,i] ]
print(node_x_ele)



########## 要素行列を求める ##########
#各線分要素の長さを計算
print("Element length")
length = np.empty(ele_total, np.float64) #各要素の長さ
for e in range(ele_total):
    length[e] = np.absolute( node_x_ele[e,1] -node_x_ele[e,0] )
print(length)

#要素行列の初期化
mat_A_ele = np.zeros((ele_total,3,3), np.float64) #要素係数行列(ゼロで初期化)
vec_b_ele = np.zeros((ele_total,3), np.float64) #要素係数ベクトル(ゼロで初期化)

#要素行列の各成分を計算
print("Local matrix")
for e in range(ele_total):
    for i in range(2):
        for j in range(2):
            mat_A_ele[e,i,j] = ( (-1)**(i+1) *(-1)**(j+1) ) / length[e]
        vec_b_ele[e,i] = -func_f *length[e]/2.0



########## 全体行列を組み立てる ##########
#全体行列の初期化
mat_A_glo = np.zeros((node_total,node_total), np.float64) #全体係数行列(ゼロで初期化)
vec_b_glo = np.zeros(node_total, np.float64) #全体係数ベクトル(ゼロで初期化)

#要素行列から全体行列を組み立てる
print("Global matrix (constructed)")
for e in range(ele_total):
    for i in range(2):
        for j in range(2):
            mat_A_glo[ node_num_glo_in_seg_ele[e,i], node_num_glo_in_seg_ele[e,j] ] += mat_A_ele[e,i,j]
        vec_b_glo[ node_num_glo_in_seg_ele[e,i] ] += vec_b_ele[e,i]

if (node_total<20): #節点数が20未満の時、全体行列を確認
    for i in range(node_total):
        for j in range(node_total):
            print("{:6.2f}".format(mat_A_glo[i,j]), end='')
        print(";{:6.2f}".format(vec_b_glo[i]))



############### 境界条件処理 ###############
print("Boundary conditions")
boundary(0, alpha, 0.0) #左端はディリクレ境界
boundary(node_total-1, "inf", beta) #右端はノイマン境界

print("Post global matrix") #節点数が20未満の時、全体行列を確認
if (node_total<20):
    for i in range(node_total):
        for j in range(node_total):
            print("{:6.2f}".format(mat_A_glo[i,j]), end='')
        print(";{:6.2f}".format(vec_b_glo[i]))



############### 連立方程式を解く ###############
print('node_total = ',node_total, ',  ele_total = ',ele_total)
#print('detA = ', scipy.linalg.det(mat_A_glo))  #Aの行列式
#print('Rank A = ', np.linalg.matrix_rank(mat_A_glo))  #AのRank(階数)
#print('Inverse A = ', scipy.linalg.inv(mat_A_glo))  #Aの逆行列

print('Unkown vector u = ')  #未知数ベクトル
unknown_vec_u = scipy.linalg.solve(mat_A_glo,vec_b_glo)  #Au=bから、未知数ベクトルuを求める
print(unknown_vec_u)
print('Max u = ', max(unknown_vec_u), ',  Min u = ',min(unknown_vec_u))  #uの最大値、最小値

#計算時間の表示
compute_time = time.time() -compute_time
print ('compute_time: {:0.5f}[sec]'.format(compute_time))



############### 計算結果を可視化 ###############
#plt.rcParams['font.family'] = 'Times New Roman'  #全体のフォントを設定
#plt.rcParams['font.size'] = 10  #フォントサイズを設定
#plt.rcParams['lines.linewidth'] = 2  # 線の太さ設定
#plt.title("Finite element analysis of 1D Poisson's equation")  #グラフタイトル
plt.xlabel("$x$")  #x軸の名前
plt.ylabel("$u(x)$")  #y軸の名前
plt.grid()  #点線の目盛りを表示

#厳密解をプロット
exact_x = np.arange(x_min,x_max,0.01)
exact_y = (func_f/2)*exact_x**2 +(-func_f*x_max +beta)*exact_x \
          -(func_f/2)*x_min**2 -(-func_f*x_max +beta)*x_min +alpha
plt.plot(exact_x,exact_y, label="$u(x)$", color='#ff0000')  #折線グラフを作成

#近似解をプロット
plt.plot(node_x_glo,unknown_vec_u, label="$\hat{u}(x)$", color='#0000ff')  #折線グラフを作成
plt.scatter(node_x_glo,unknown_vec_u)  #点グラフを作成

#更に体裁を整える
plt.axis('tight') #余白を狭くする
plt.axhline(0, color='#000000')  #x軸(y=0)の線
plt.axvline(0, color='#000000')  #y軸(x=0)の線
plt.legend(loc='best')  #凡例(グラフラベル)を表示
for n in range(node_total):  #節点番号をグラフ内に表示
    plt.text(node_x_glo[n],unknown_vec_u[n], n, ha='center',va='bottom', color='#0000ff')

plt.show()  #グラフを表示
#plt.savefig('fem1d_poisson.pdf')
#plt.savefig('fem1d_poisson.png')