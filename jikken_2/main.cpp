#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime> 
#include <random>
#include <algorithm> //sort用
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
        //cout<<consensus.size()<<endl;
        //cout<<consensus[0].size()<<endl;
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
                cout<<consensus[i][j]<<" ";
            }
            cout << endl;
        }
        

    //３．コンセンサス表から、i番目に塩基Xが出現する確率をprob表に。
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
        }
        cout << endl;

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
    return ozscore;
}

//d.スコア行列を用いた結合部位の探索。別のfunctionを作る。
double CalcHit(vector<vector<double>>ozscore, string sequence){
    double hit=0;
  
    for(int j=0; j<sequence.size(); j++){
        //cout << "a" << endl;
        if(sequence[j]=='A'){
            hit +=ozscore[0][j];
        }
        else if(sequence[j]=='C'){
            hit +=ozscore[1][j];
        }
        else if(sequence[j]=='G'){
            hit +=ozscore[2][j];
        }
        else if(sequence[j]=='T'){
            hit +=ozscore[3][j];
        }
    }

    return hit;
}

//position番号を保存しておけば、あとから配列はそのwindowの配列は取得できるので、一緒に保存しなくていい。
struct Scores {
    int position;
    double score;
};

//引数にozscoreを入れることで、その中の情報を使えるようになる。
void RunWindow(string promoter, vector<vector<double>> ozscore,vector<double> &scores){
    int length = ozscore[0].size();

    for(int i = 0; i < promoter.size() - length + 1; i++){
        //.substr(i番目から,n文字)　取得
        string window = promoter.substr(i, length);
        scores[i] = CalcHit(ozscore, window);
    }
}

//ランダムに塩基配列を作る
string RandomBase(double BaseLength, mt19937 mt, vector<vector<double>> ozscore){
    string rand_sequence;
    uniform_real_distribution<> dist(0,1);
    
    for(int i=0; i<BaseLength; i++){
        //毎回新しいランダム値を出してもらう
        double random = dist(mt);

        if(random<0.25){
            rand_sequence.push_back('A');
        }
        else if(random<0.50){
            rand_sequence.push_back('C');
        }
        else if(random<0.75){
            rand_sequence.push_back('G');
        }
        else if(random<1.0){
            rand_sequence.push_back('T');
        }
    }
    return rand_sequence;
}

void result_output(vector<string> promoters_set, vector<string> promoters_name_sets, vector<vector<double>> ozscore, double threshold){
    for(int i=0; i<promoters_set.size(); i++){//複数のプロモーターがある
        double score; 
        int length = ozscore[0].size();
        //一つのプロモーターについて。その文字列
        for(int j = 0; j < promoters_set[i].size() - length + 1; j++){
            //.substr(i番目から,n文字)　取得
            string window = promoters_set[i].substr(j, length);
            score = CalcHit(ozscore, window);
            if(score > threshold){
                cout << "Promoter name " << promoters_name_sets[i] << endl;
                cout << "Sequence [" << window << "]  Score: " << score << endl;
                cout << "Position " << j+1 << endl;
                cout << endl;
            }
        }
        cout << "/////////////////////" << endl;
    }
}



int main(){

    //各モチーフに対してここは手動で変える。
    cout << "Motif name: " << "MATa1" << endl;

    //プロモーターの名前、配列をvectorに書き込む。
    ifstream ifs("../promoters");
    if(!ifs){
            cerr << "cant" << endl;
            exit(1);
    }

    vector<string> promoters_set;
    vector<string> promoters_names_set;

    for(int i =0; i<8; i++){
        string promo_name;
        getline(ifs, promo_name);
        promoters_names_set.push_back(promo_name);

        string promoter;
        getline(ifs, promoter);
        promoters_set.push_back(promoter);
    }


    //タンパク質に特有のozscore、ファイルのインポートはmainではなく、ozscoreで作る。
    vector<vector<double>> ozscore;
    ozscore = CalcOzscore();


    //ランダムの塩基配列を作る。
    double BaseLength; 
    cout << "Enter Base Length" << endl; 
    cin >> BaseLength;


    random_device rnd;
    mt19937 mt(rnd());

    string random_sequence;
    random_sequence = RandomBase(BaseLength, mt, ozscore);
    vector<double> Rand_scores(random_sequence.size());
    /*random_sequnece配列にたいしてwindowを行って、score達だけを配列に保存。けどさっきの使うから、位置も保存されるわ。*/
    RunWindow(random_sequence, ozscore, Rand_scores);
    
    //閾値を計算する
    double pvalue; 
    cout << "Enter p-value" << endl; 
    cin >> pvalue;
    cout << endl;

    //random_sequnece配列についてwindowを回して得られたscoreをsortした、上からtop番目のscoreを閾値に設定。
    sort(Rand_scores.begin(), Rand_scores.end());
    reverse(Rand_scores.begin(), Rand_scores.end());
    
    int top;
    top = pvalue*BaseLength;
    
    //閾値以上のものをDetThreshold行列に入れる。
    double threshold = Rand_scores[top-1];
    result_output(promoters_set,promoters_names_set,ozscore,threshold);
    return 0;
}


