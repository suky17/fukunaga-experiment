#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;



vector<vector<double>> CalcOzscore(){
    //１．コンセンサス表（頻度表）を作成
        ifstream ifs ("../MATa1");//jikken_2の外側のファイルだから..つける
        if(!ifs){
            cerr << "cant" << endl;
            exit(1);
        }

        string seq;
        
        string first;
        getline(ifs, first);//配列のlengthをとってきて、表の横を設定するのに使う。
        double length=first.size();

        vector<vector<double>> consensus(4, vector<double>(length, 1)); //初期化ですべての表要素を１にする（最後の＋１不要）
        
        //一行目の配列は、getlineしてfirstに格納してしまったので、一行目は別でやってあげる。
        
        for(int i=0; i<length; i++){
                if(first[i]=='A'){
                    consensus[0][i]++;
                }
                else if(first[i]=='C'){
                    consensus[1][i]++;
                }
                else if(first[i]=='G'){
                    consensus[2][i]++;
                }
                else if(first[i]=='T'){
                    consensus[3][i]++;
                }
            }

        //２行目以降はloopで。
        int seqcount = 1;//今後のために、配列の本数を数えておく。(一行目は上でやってしまっているので。)
        while(getline(ifs, seq)){ //ifsの2行目から読んでseqに格納。行がenterされるごとにloopが回る。
            for(int i=0; i<seq.size(); i++){
                if(seq[i]=='A'){
                    consensus[0][i]++;
                }
                else if(seq[i]=='C'){
                    consensus[1][i]++;
                }
                else if(seq[i]=='G'){
                    consensus[2][i]++;
                }
                else if(seq[i]=='T'){
                    consensus[3][i]++;
                }
            }
            seqcount++;
        }

        //コンセンサ表を表示
        cout << "consensus" << endl;
        for(int i=0; i<4; i++){
            for(int j=0; j<length; j++){
                cout<<consensus[i][j];
            }
            cout << endl;
        }
        

    //３．コンセンサス表から、i番目に塩基A、、が出現する確率をprob表に。
        vector<vector<double>> prob(4, vector<double>(length, 0));
        for(int i = 0; i<4; i++){
            for(int j=0; j<length; j++){
                prob[i][j]=consensus[i][j]/(seqcount+4);//+4は初期設定の１が各要素に入ってるから
            }
        }


        //prob表を表示。手順３
        cout << endl;
        cout << "prob" << endl;
        for(int i=0; i<4; i++){
            for(int j=0; j<length; j++){
                cout<<prob[i][j];
            }
            cout << endl;
            cout << endl; 
        }


    //４．塩基ｘのバックグラウンド出現確率q(x)を計算
        //出芽酵母の全ゲノムに含まれる塩基数
        double ATsum = 7519429;
        double CGsum = 4637676;
        double sum = ATsum*2 + CGsum*2;

        vector<double> background(4, 0);
        background[0]=ATsum/sum;
        background[1]=CGsum/sum;
        background[2]=CGsum/sum;
        background[3]=ATsum/sum;

        //background配列表示
        cout << "background" << endl;
        for(int i = 0; i<background.size(); i++){
            cout << background[i] << endl;
        }
        cout << endl;


    //５．対数オッズスコア表の作成
        vector<vector<double>> ozscore(4, vector<double>(length, 0));

        for(int i=0; i<4; i++){
            for(int j=0; j<length; j++){
                ozscore[i][j]=log(prob[i][j]/background[i]);
            }
        }

        //ozscore表の表示
        cout << "ozscore" << endl;
        for(int i=0; i<4; i++){
            for(int j=0; j<length; j++){
                cout<<ozscore[i][j];
            }
            cout << endl;
        }
        cout << endl;


    ifs.close();
    return ozscore();
}

//d.スコア行列を用いた結合部位の探索。別のfunctionを作る。
double CalcHit(vector<vector<double>>ozscore, string sequence){
    double hit=0;

    for(int j=0; j<sequence.size(); j++){
        if(sequence[j]=='A'){
            hit +=ozscore[0][j];
        }
        else if(sequnece[j]='C'){
            hit +=ozscore[1][j];
        }
        else if(sequnece[j]='G'){
            hit +=ozscore[2][j];
        }
        else if(sequnece[j]='T'){
            hit +=ozscore[3][j];
        }
    }

    return hit;
}

int main(){

    vector<vector<double>> ozscore;
    CalcOzscore;

    string sequence = "AAAAAAAA";
    CalcHit(ozscore, sequence);

//本来は転写因子ごとにozscore表を作成すので、fileも＆してretrieveするわ。

    return 0;
}


