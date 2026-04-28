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

int NUM_FEATURES = 53;
int NUM_SEQS = 10000;

struct TreeNode{
    int feature_id;
    double threshold;
    int left_class_id;
    int right_class_id;
};


//&参照だと、メインの中で呼んだ時に、直接書き込んでくれる。
void LoadSolubilityFile(string filename, vector<string> &feature_name, vector<vector<double>> &dataset, vector<int> &labels){
    ifstream ifs (filename);//ここで"filename"のquotationは不要。stringで定義されているから。
    if(!ifs){
        cerr << "can't open file" << endl;
        exit(1);
    }

    string id;
    ifs >> id; //不要
    for(int i=0; i<53; i++){
        ifs >> feature_name[i];
    }

    string label;
    ifs >> label; //不要

    for(int i=0; i<10000; i++){
        //ここで勝手に改行
        string protein_number;
        ifs >>  protein_number; //不要

        for(int j=0; j<53; j++){
            ifs >> dataset[i][j];
        }
        ifs >> labels[i];
    }

    return;
    
}


//mainの中でrandomは使わないから、引数にmt入れなくてもいい。
//＆参照だとメインの中でも共通で使える
void DivideDataset(vector<vector<double>> dataset, vector<int> &labels, vector<vector<double>> &training_dataset, vector<int> &training_labels, vector<vector<double>> &test_dataset, vector<int> &test_labels, double &test_ratio){
    mt19937 mt(0);
    //シャッフルするための土台の0～9999をつくる。
    vector<int> temp;
    for(int i=0; i<10000; i++){
        temp.push_back(i);
    }
    //これでランダムに並び替える
    shuffle(temp.begin(), temp.end(), mt);

    //〇番までテストに入れる
    int test_number = test_ratio*10000;
    //大きさがメインの中で定義されていないので、定義してあげる。
    //定義ではないので、ここで大きさを確保するには .resize 
    training_dataset.resize(10000-test_number, vector<double>(53));
    training_labels.resize(10000-test_number);
    test_dataset.resize(test_number, vector<double>(53));
    test_labels.resize(test_number);


    for(int i =0; i<test_number; i++){
        for(int j =0; j<53; j++){
            test_dataset[i][j] = dataset[temp[i]][j];
        }
        test_labels[i]=labels[temp[i]];
    }

    for(int i=0; i<10000-test_number; i++){
        for(int j =0; j<53; j++){
            training_dataset[i][j]=dataset[temp[i+test_number]][j];
        }
        training_labels[i]=labels[temp[i+test_number]];
    }    
    return;
}

//テストデータを流し込んで、閾値をもととに左か右かで判定していくシステム。
void Evaluation(vector<TreeNode> &decision_tree, vector<vector<double>> test_dataset, vector<int> test_labels, double &test_ratio){
    int test_number = test_ratio*10000;
    double TP=0;
    double FP=0;
    double FN=0;
    double TN =0;

    vector<int> pred_label(test_labels.size());

    //branch 1,2,3で３回判定
    for(int b =0; b<3; b++){
        for(int i =0 ;i<test_labels.size(); i++){
            if(test_dataset[i][decision_tree[b].feature_id] <= decision_tree[b].threshold){
                //leftと判定された＝０ではなくtrainingデータの左にいるものの代表値をとる。
                pred_label[i]=decision_tree[b].left_class_id;
            }
            else if(test_dataset[i][decision_tree[b].feature_id] > decision_tree[b].threshold){
                //right 
                pred_label[i]=decision_tree[b].right_class_id;
            }
        }

        for(int i=0; i<test_number; i++){
            if(test_labels[i]==1 && pred_label[i]==1){
                TP++;
            }
            else if(test_labels[i]==0 && pred_label[i]==1){
                FP++;
            }
            else if(test_labels[i]==1 && pred_label[i]==0){
                FN++;
            }
            else if(test_labels[i]==0 && pred_label[i]==0){
                TN++;
            }
        }
    }                                               
    //intで定義してるから、計算はdoubleで行ってね。
    double accuracy=(TP+TN)/(TP+FP+FN+TN);
    double precision=TP/(TP+FP);
    double recall=TP/(TP+FN);
    double f_score=(2*precision*recall)/(recall+precision);
    
    cout << endl;
    //cout << TP << " " << FP << " " << FN << " " << TN << endl;
    cout << "Accuracy: " << accuracy << " Precision: " << precision << " Recall: " << recall << endl;
    cout << "f-score: " << f_score << endl;
}


void TrainDecisionNode(vector<vector<double>> training_dataset, vector<int> training_labels, TreeNode &decision_tree){
    //一つのfeatureについての値を大きい順に並べる。縦のボックスをtemp_featureとする。
    vector<double> temp_feature(training_labels.size());
    double min =10000000000;
    int left_class_id =0;
    int right_class_id=0;

    //横→にずれていく
    for(int j=0; j<53; j++){
        //cout << "j: " << j << " " << endl;
        //一つの縦列について
        for(int i=0; i<training_labels.size(); i++){
            //temporary格納する
            temp_feature[i]=training_dataset[i][j];
        }
        //大きい順に並べる
        sort(temp_feature.begin(), temp_feature.end());

        //9999の区切りで計算するのは難しいので、1～99パーセンタイルという方法を使う。すると、区切りが99回で済む。
        

        for(int c=0; c<99; c++){

            //○○パーセント番目のdata_setの値
            double threshold = temp_feature[c*0.01*training_labels.size()];
            double L0=0;
            double L1=0;
            double R0=0;
            double R1=0;

            for(int n=0; n<training_labels.size(); n++){
                if(training_dataset[n][j]<threshold){
                    if(training_labels[n]==0){
                        L0++;
                    }
                    else{
                        L1++;
                    }

                }
                else{
                    if(training_labels[n]==0){
                        R0++;
                    }
                    else{
                        R1++;
                    }
                }
            }

            double pL = 0.0;
            double pR=0.0; 
            pL=L1/(L0+L1);
            pR=R1/(R0+R1);
            
            double GL= 0.0;
            double GR = 0.0;

            GL = 2.0*pL*(1.0-pL);
            GR = 2.0*pR*(1.0-pR);

            double Jini = 0;

            Jini = (L0+L1)/(training_dataset.size())*GL + (R0+R1)/(training_dataset.size())*GR;
            
            if(Jini<min){

                min=Jini;
                decision_tree.feature_id = j;
                decision_tree.threshold = threshold;
                if(L0>L1){
                    left_class_id=0;
                }
                else if(L0<L1){
                    left_class_id=1;
                }
                if(R0>R1){
                    right_class_id=0;
                }
                else if(R0<R1){
                    right_class_id=1;
                }
                
                decision_tree.left_class_id = left_class_id;
                decision_tree.right_class_id = right_class_id; 
            }
        }
    }
    cout <<"ジニ不純度が最も小さくなる最適な特徴量：" << decision_tree.feature_id << " 閾値：" << decision_tree.threshold << endl;
    cout << "left_class_id: " << decision_tree.left_class_id << " right_class_id: "<< decision_tree.right_class_id << endl;
    cout << endl;
}



void TrainDecisionTree(vector<vector<double>> training_dataset, vector<int> training_labels, vector<TreeNode> &decision_tree){
    //この中でTraindecisionnodeを3回呼び出す                
    TrainDecisionNode(training_dataset, training_labels, decision_tree[0]);

    vector<vector<double>> training_dataset1;
    vector<vector<double>> training_dataset2;

    vector<int> training_labels1;
    vector<int> training_labels2;

    //上からｎ番目のプロテインについて。左か右に分類する。
    for(int i=0; i<training_dataset.size(); i++){
        if(training_dataset[i][decision_tree[0].feature_id] <= decision_tree[0].threshold){
            ///左に行く要素だけで、新たに作る。
            training_dataset1.push_back(training_dataset[i]);
            training_labels1.push_back(training_labels[i]);
        }
        else if(training_dataset[i][decision_tree[0].feature_id] > decision_tree[0].threshold){
            ///右に行く要素だけで、新たに作る。
            training_dataset2.push_back(training_dataset[i]);
            training_labels2.push_back(training_labels[i]);
        }
    }

    //左に分類されたdatasets1について、同じTrainDecisionNode判定を実装する
    TrainDecisionNode(training_dataset1, training_labels1, decision_tree[1]);
    TrainDecisionNode(training_dataset2, training_labels2, decision_tree[2]);

}



int main(void){
    //上のカテゴリー５３個
    vector<string> feature_name(NUM_FEATURES, "");
    vector<vector<double>> dataset(NUM_SEQS, vector<double>(NUM_FEATURES, 0.0));
    vector<int> labels(NUM_SEQS);
    
    //○○＝になってないから返り値が無いことが分かる。
    LoadSolubilityFile("protein_solubility_dataset.txt", feature_name, dataset, labels);

    vector<vector<double>> training_dataset;
    vector<int> training_labels;
    vector<vector<double>> test_dataset;
    vector<int> test_labels;
    double test_ratio = 0.2;

    DivideDataset(dataset, labels, training_dataset, training_labels, test_dataset, test_labels, test_ratio);

    vector<TreeNode> decision_tree(3);

    TrainDecisionTree(training_dataset, training_labels, decision_tree);

    Evaluation(decision_tree, test_dataset, test_labels, test_ratio);


    return 0;
}