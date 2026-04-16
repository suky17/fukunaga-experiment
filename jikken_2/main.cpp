#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

int main(){
    
    ifstream ifs ("MATa1");
    if(!ifs){
        cerr << "cant" << endl;
        exit(1);
    }
    char base;
    string file;
    cout << 2 << endl;
    while(getline(ifs, file)){
        cout << file.size() << endl;
        cout << 1 << "file" << endl;
    }
    cout << file.size() << endl;
    cout << 2 << endl; //ifsの１行目を読んでfileに格納
    ifs.close();
}

