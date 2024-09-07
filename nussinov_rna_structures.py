import numpy as np



def main():
    ## RNA配列
    sequence = 'UCGUUCGGUAACA' # 課題
    #sequence = 'CGAUCGUACUAAGUGC' # デバッグ用

    ## 動的計画法
    W, W_map = nussinov(sequence) # W[i][j] = -(sequence[i:j+1] に含まれる塩基対の最大の数)
    print()
    print('W =')
    print(W)
    print()

    ## 塩基対の情報
    pairs_info = identfyPairs(W, W_map) # 動的計画法と逆の動きで計算していく
    
    ## 2次元構造の表示
    if len(pairs_info) == 0:
        print('Structure does not exist.')
    else:
        for i in range(len(pairs_info)):
            print(f'------ structure {i+1} ------')
            print()
            print(printStructure(sequence, pairs_info[i]))
            print()
    

    return





def printStructure(sequence, pairs_info):
    '''配列の二次構造を可視化するため，文字列として出力する

    Args:
        sequence (str): RNA配列 (N)
        pairs_info (numpy ndarray): RNA配列のペア情報を格納した配列 (N)
    
    Returns:
        result (str): 二次構造を可視化する文字列
    '''

    # pairs_info の各要素には，sequenceの各要素に対応するペアのインデックスが格納されている
    # ただし，対応するペアがない場合には-1が格納されている
    #
    # 入力例）
    # sequence = 'UCGUUCGGUAACA', pairs_info = np.array([12, -1, 11, 10, 9, 7, -1, 5, -1, 4, 3, 2, 0])
    #
    # 出力例)
    # U C G U U C G G U A A C A
    # |   | | | |   |   | | | |
    # |   | | | ¯¯¯¯¯   | | | |
    # |   | | ¯¯¯¯¯¯¯¯¯¯¯ | | |
    # |   | ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ | |
    # |   ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ |
    # ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯

    ## sequence と pairs_info は長さが同じ
    if len(sequence) != len(pairs_info):
        return 'エラー1 (pairs_info の長さが sequence の長さと一致しません)'
    
    ## pairs_info に -1 でない同じ数字は2つ以上あらわれない
    for p in pairs_info:
        if len(pairs_info[pairs_info==p]) >= 2 and p != -1:
            return 'エラー2 (pairs_info に -1 でない同じ数字が複数含まれています)'


    ## 見やすさのためにRNA配列 sequence の間に空白を入れる
    ## 例）'UCGU' => 'U C G U'
    sequence = ' '.join(sequence) 
    
    ## sequence に合わせて pairs_info も整形する
    tmp = np.full(len(pairs_info)*2-1, -1,int)
    tmp[::2] = pairs_info * 2
    tmp[tmp==-2] = -1
    pairs_info = tmp

    ## 設定
    result = sequence + '\n' # 出力
    n_pairs = len(pairs_info[pairs_info!=-1]) // 2 # (A,U), (G,C)ペアの数
    n_sequence = len(sequence) # sequence の長さ

    ## 出力の作成
    for i in range(n_pairs + 1):
        j = 0

        if np.all(pairs_info == -1): # 全ての塩基対を書き込めたら終了
            break

        while j < n_sequence:
            ind = pairs_info[j]

            if ind == -1:
                result += ' '
            
            else:
                if  i != 0 \
                    and j != n_sequence - 1 \
                    and ind > j \
                    and (all(element == -1 for element in pairs_info[j+1:ind]) \
                         or ind == j + 2):

                    result += '¯' * (ind + 1 - j)
                    pairs_info[j:ind+1] = -1
                    j += ind - j
                
                else:
                    result += '|'
            
            j += 1

        result += '\n'


    return result





def isBasePair(x, y):
    '''塩基対か否か判定

    Args:
        x (str): 塩基1
        y (str): 塩基2
    
    Returns:
        塩基対か(True)否か(False)
    '''
    
    if (x == 'A' and y == 'U') or (x == 'U' and y == 'A'):
        return True
    
    elif (x == 'G' and y == 'C') or (x == 'C' and y == 'G'):
        return True
    
    else:
        return False





def nussinov(sequence):
    '''RNA配列の(A,U)(G,C)ペアの数が最大となる疑似ノットなしの２次構造を動的計画法アルゴリズム(NussinovAlgorithm)で求める

    Args:
        sequence (str): RNA配列 (N)
    
    Returns:
        W (numpy ndarray): 動的計画法アルゴリズムで求めた行列 (N, N)
        W_map (numpy ndarray): 動的計画の動きを記録した行列 (N, N)
    '''

    ## 設定
    n_sequence = len(sequence) # 配列の長さ
    W = np.zeros((n_sequence, n_sequence), int)
    W_map = [[None for _ in range(n_sequence)] for _ in range(n_sequence)]

    ## 塩基対を持つRNA配列は最低でも3つの塩基が必要
    if n_sequence < 3:
        return W, W_map

    ## 動的計画法
    for i in range(n_sequence-3, -1, -1):
        for j in range(i+2, n_sequence):
            value = np.ones(4)
            value[0] = W[i+1, j].copy()
            value[1] = W[i, j-1].copy()
            if isBasePair(sequence[i], sequence[j]):
                value[2] = W[i+1, j-1] - 1

            if i+4 < j:
                W_tmp = np.array([W[i, k] + W[k+1, j] for k in range(i+2, j-2)])
                k_tmp = np.where(W_tmp == np.min(W_tmp))[0] + i + 2
                k = []

                for kk in k_tmp:
                    if W[i, kk] != 0 and W[kk+1, j] != 0:
                        k.append(kk)
                
                if len(k) > 0:
                    value[3] = np.min(W_tmp)
            
            min_ind = np.where(value == np.min(value))[0]
            if (0 in min_ind) and (value[0] == 0):
                min_ind = [0]
            
            W[i, j] = value[min_ind[0]]
            W_map[i][j] = min_ind


    return W, W_map





def identfyPairs(W, W_map):
    '''動的計画法アルゴリズム(NussinovAlgorithm)で求めた行列から塩基対を逆算する

    Args:
        W (numpy ndarray): 動的計画法アルゴリズムで求めた行列 (N, N)
        W_map (numpy ndarray): 動的計画の動きを記録した行列 (N, N)
    
    Returns:
        pairs_info (numpy ndarray): RNA配列のペア情報を格納した配列 (?, N)
    '''

    ## 設定
    n_sequence = len(W) # 配列の長さ
    pairs_info = [np.full(n_sequence, -1, int)]

    ## 逆算用の関数            
    def component(n=0, i=0, j=n_sequence-1):
        if i + 2 > j:
            return
        
        pairs_info_n = pairs_info[n].copy()

        for next in W_map[i][j]:
            
            if next == 0:
                component(n, i+1, j)

            elif next == 1:
                if next == W_map[i][j][0]:
                    component(n, i, j-1)

                elif W[i, j-1] != 0:
                    pairs_info.append(pairs_info_n.copy())
                    component(len(pairs_info)-1, i, j-1)

            elif next == 2:
                if next == W_map[i][j][0]:
                    pairs_info[n][i] = j
                    pairs_info[n][j] = i
                    component(n, i+1, j-1)

                else:
                    pairs_info.append(pairs_info_n.copy())
                    pairs_info[len(pairs_info)-1][i] = j
                    pairs_info[len(pairs_info)-1][j] = i
                    component(len(pairs_info)-1, i+1, j-1)

            elif next == 3:
                W_tmp = np.array([W[i, k] + W[k+1, j] for k in range(i+2, j-2)])
                k = np.where(W_tmp == np.min(W_tmp))[0] + i + 2

                for kk in k:
                    if next == W_map[i][j][0] and kk == k[0]:
                        len_before = len(pairs_info)
                        component(n, i, kk)
                        len_after = len(pairs_info)
                        component(n, kk+1, j)
                        
                        if len_before != len_after:
                            for nn in range(len_before, len_after):
                                component(nn, kk+1, j)

                    else:
                        pairs_info.append(pairs_info_n.copy())
                        len_before = len(pairs_info)
                        component(len_before-1, i, kk)
                        len_after = len(pairs_info)

                        for nn in range(len_before-1, len_after):
                            component(nn, kk+1, j)

        return
    
    ## 逆算による塩基対の計算
    component()
    pairs_info = np.unique(pairs_info, axis=0) # 被りの消去

    ## 塩基対が見つからなかったとき
    if np.all(pairs_info == -1):
        return []


    return pairs_info

    



if __name__ == "__main__":
    main()