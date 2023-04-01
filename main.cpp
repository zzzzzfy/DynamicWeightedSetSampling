#include <unordered_map>
#include <random>
#include <fstream>
#include <string>
#include <cstring>
#include "ScapegoatTree.hpp"
#include "ChunkScapegoatTree.hpp"
#include "value_bucket.hpp"
#include "BF.hpp"
#include "BSTSampling.hpp"
#include "utility.h"
#include "basic_alias.hpp"
#include "linear_search.hpp"
#include "size_block.hpp"
#include "basic_bst.hpp"
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
using namespace std;

//Element element[50000000];
Element *element;
int key_max = 1000000000;
int value_max = 1000;
int SAMPLE_TIMES = 10000;
int csize;
// U max weight value, dis_type weight distribution type

void WIRS_sgt_test(int n, Element* ele_vec, vector<Opt>& opt_vec) {
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();

    sort(ele_vec, ele_vec+n, cmp_element_key);
    ScapegoatTree<sizeBlockMethod> scapegoatTree;
    scapegoatTree.init(ele_vec, n);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "sgt tree build_time:" << elapsed_seconds.count() << "\n";
    start = chrono::system_clock::now();
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            scapegoatTree.ask(i.key, i.value, i.weight);
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            scapegoatTree.insert(tmp_ele);
        }
        if (i.opt_type == 2) {
            scapegoatTree.erase(i.key);
        }
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    puts("");
    double cnt = opt_vec.size();
    cout << "sgt Operation_time_average:" << elapsed_seconds.count()/cnt<<endl;
    cout << "sgt Operation_time:" << elapsed_seconds.count() << "\n";
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024 ;
    cout <<"sgt_memory " << gms << "mb\n";
}

void WIRS_sgt_bucket_test(int n, Element* ele_vec, vector<Opt>& opt_vec,double U) {
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();

    sort(ele_vec, ele_vec+n, cmp_element_key);
    ScapegoatTree<valueBucketMethod> scapegoatTree;
    scapegoatTree.init(ele_vec, n,U);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "sgt tree build_time:" << elapsed_seconds.count() << "\n";
    start = chrono::system_clock::now();
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            scapegoatTree.ask(i.key, i.value, i.weight);
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            scapegoatTree.insert(tmp_ele);
        }
        if (i.opt_type == 2) {
            scapegoatTree.erase(i.key);
        }
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    puts("");
    double cnt = opt_vec.size();
    cout << "bucket sgt Operation_time_average:" << elapsed_seconds.count()/cnt<<endl;
    cout << "bucket sgt Operation_time:" << elapsed_seconds.count() << "\n";
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024 ;
    cout <<"sgt_memory " << gms << "mb\n";
}

void WIRS_chunk_bucket_test(int n, Element* ele_vec, vector<Opt>& opt_vec,double U) {
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();

    sort(ele_vec, ele_vec + n, cmp_element_key);
    ChunkScapegoatTree<valueBucketMethod> scapegoatTree;
    scapegoatTree.setChunkSize(csize);
    scapegoatTree.init(ele_vec, n,U);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "sgt chunk build_time:" << elapsed_seconds.count() << "\n";

    int tot_query = 0;
    double tot_sum = 0;
    double cnt= 0;
    start = chrono::system_clock::now();
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            vector<int> ans;
            scapegoatTree.ask(i.key, i.value, i.weight);
            cnt+=1;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            scapegoatTree.insert(tmp_ele);
        }
        if (i.opt_type == 2) {
            scapegoatTree.erase(i.key);
        }if(i.opt_type == 0 && cnt == 50) {
            end = chrono::system_clock::now();
            elapsed_seconds = end - start;
            cout << "chunk bucket Operation_time_average:" << elapsed_seconds.count() / cnt << endl;
            cout << "chunk bucket Operation_time:" << elapsed_seconds.count() << "\n";
            start = chrono::system_clock::now();
            cnt = 0;
            puts("");
        }
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "chunk bucket Operation_time:" << elapsed_seconds.count() << "\n";
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cout << "chunk_memory " << gms << "mb\n";
    //cout << "sizeBlock:" << tot_sum / tot_query << "\n";
    puts("-----------------------------------");
}

void WIRS_chunk_test(int n, Element* ele_vec, vector<Opt>& opt_vec) {
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();

    sort(ele_vec, ele_vec + n, cmp_element_key);
    ChunkScapegoatTree<sizeBlockMethod> scapegoatTree;
    scapegoatTree.setChunkSize(csize);
    scapegoatTree.init(ele_vec, n);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cerr<<"buildOK"<<endl;
    cout << "sgt chunk build_time:" << elapsed_seconds.count() << "\n";

    int tot_query = 0;
    double tot_sum = 0;
    double cnt = 0;
    start = chrono::system_clock::now();
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            vector<int> ans;
            scapegoatTree.ask(i.key, i.value, i.weight);
            cnt+=1;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            scapegoatTree.insert(tmp_ele);
        }
        if (i.opt_type == 2) {
            scapegoatTree.erase(i.key);
        }
        if(i.opt_type == 0 && cnt == 50) {
            end = chrono::system_clock::now();
            elapsed_seconds = end - start;
            cout << "chunk block Operation_time_average:" << elapsed_seconds.count() / cnt << endl;
            cout << "chunk block Operation_time:" << elapsed_seconds.count() << "\n";
            start = chrono::system_clock::now();
            cnt = 0;
            puts("");
        }
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "chunk block Operation_time:" << elapsed_seconds.count() << "\n";
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cout << "chunk_memory " << gms << "mb\n";
    //cout << "sizeBlock:" << tot_sum / tot_query << "\n";
    puts("-----------------------------------");
}

void WIRS_bstTest(int n, Element* ele_vec, vector<Opt>& opt_vec) {
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();

    BSTSampling bst;

    sort(ele_vec, ele_vec + n, cmp_element_key);
    bst.init(ele_vec, n);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "bst tree build_time:" << elapsed_seconds.count() << "\n";

    int tot_query = 0;
    double tot_sum = 0;
    start = chrono::system_clock::now();
    double cnt = 0;
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            bst.ask(i.key, i.value, i.weight);
            cnt+=1;
            //tot_sum += tmp_method.random_sample_value();
            //tot_query++;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            bst.insert(tmp_ele);
            //tmp_method.insert_element(tmp_ele);
        }
        if (i.opt_type == 2) {
            bst.erase(i.key);
            //tmp_method.delete_element(i.key);
        }
        if(i.opt_type == 0 && cnt == 50) {
            end = chrono::system_clock::now();
            elapsed_seconds = end - start;
            cout << "bst Operation_time_average:" << elapsed_seconds.count() / cnt << endl;
            cout << "bst Operation_time:" << elapsed_seconds.count() << "\n";
            start = chrono::system_clock::now();
            cnt = 0;
            puts("");
        }
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "bst Operation_time:" << elapsed_seconds.count() << "\n";
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cout << "bst_memory " << gms << "mb\n";
    puts("-----------------------------------");
}

void WIRS_produce_data(int tot_size, int opt_size, int ins_del_ratio, int query_ratio, int coverage, int sample_times,string filename,string file_s) {
    ifstream inFile;
    inFile.open(filename.c_str());

    ofstream tmpf;
    tmpf.open(file_s.c_str());

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> key_random(0, tot_size);

    set<int> selected_key;
    if(query_ratio == 0)
    {
        for(int i=1;i<=opt_size;i++)
        {
            int now = key_random(gen);
            while(selected_key.find(now)!=selected_key.end())
                now = key_random(gen);
            selected_key.insert(now);
        }
    }
    if(query_ratio == 0&&ins_del_ratio==100) tmpf<<tot_size-opt_size<<'\n';
    else tmpf<<tot_size<<'\n';
    int mx = 0;
    vector<Element> opt;
    for(int i = 0;i<tot_size;i++)
    {
        int key,value,weight;
        inFile>>key>>value>>weight;
        if(query_ratio == 0&&selected_key.find(i)!=selected_key.end()){
            opt.push_back(Element(key,value,weight));
            if(ins_del_ratio==100) continue;
        }
        tmpf << key << " " << value << " " << weight << "\n";
        mx = max(weight,mx);
    }
    cout<<mx<<endl;
    uniform_int_distribution<int> opt_random(1, 100);

    if(query_ratio == 100){
        tmpf << opt_size * 8 << '\n';
        sample_times = 100000;
        for(int j = 1;j<=4;j++) {
            coverage = 20*j;
            int offset = tot_size * (coverage / 100.0);
            for (int i = 0; i < opt_size; i++) {
                double l = key_random(gen);
                double r = l + offset;
                while (r > tot_size) {
                    l = key_random(gen);
                    r = l + offset;
                }
                int t = sample_times;
                int L = l,R = r;
                tmpf << "0" << " " << L << " " << R << " " << t << "\n";
            }
        }
        sample_times = 100;
        coverage = 50;
        for(int j = 1;j<=4;j++) {
            sample_times *= 10;
            int offset = tot_size * (coverage / 100.0);
            for (int i = 0; i < opt_size; i++) {
                double l = key_random(gen);
                double r = l + offset;
                while (r > tot_size) {
                    l = key_random(gen);
                    r = l + offset;
                }
                int t = sample_times;
                int L = l,R = r;
                tmpf << "0" << " " << L << " " << R << " " << t << "\n";
            }
        }
    }else {
        tmpf << opt_size << '\n';
        for (int i = 0; i < opt_size; i++) {
            if (ins_del_ratio == 100) {
                tmpf << "1 " << opt[i].key << " " << opt[i].value << " " << opt[i].weight << '\n';
            } else tmpf << "2 " << opt[i].key << "\n";
        }
    }
//    tmpf<<opt_size<<'\n';
//
//    int offset = tot_size*(coverage/100.0);
//    for (int i = 0; i < opt_size; i++) {
//        if (query_ratio==100) {
//            int l = key_random(gen);
//            int r = l+offset;
//            while(r>tot_size){
//                l = key_random(gen);
//                r = l+offset;
//            }
//            int t = sample_times;
//            tmpf << "0" << " " << l << " " << r <<" "<<t << "\n";
//        }
//        else if(ins_del_ratio == 100)
//        {
//            tmpf<<"1 "<<opt[i].key<<" "<<opt[i].value<<" "<<opt[i].weight<<'\n';
//        }
//        else tmpf<<"2 "<<opt[i].key<<"\n";
//    }
}

void WIRS_produce(int n, int m, int ins_del_ratio, int query_ratio, int coverage,int sample_times,string filename) {

    string file_s = "/home/DynamicSetSampling/data/WIRS_n" +to_string(n) + "_m" + to_string(m) + "_r1_" +
                    to_string(ins_del_ratio) + "_r2_" + to_string(query_ratio) + "_coverage"+ to_string(coverage)+"_times"+
                                                                                                                  to_string(sample_times);
    WIRS_produce_data(n, m, ins_del_ratio, query_ratio, coverage,sample_times,filename,file_s);
}

void WIRS_test_instance(string file_s,double U,int op) {
    ifstream tmpf;
    int n,m;
    vector<Element> ele_vec;
    vector<Opt> opt_vec;
    cerr << file_s << "\n";
    tmpf.open(file_s.c_str());
    int _k, _v; float _w;
    tmpf >> n;
    for (int i = 0; i < n; i++) {
        tmpf >> _k >> _v >> _w;
        //cerr << _k <<endl;
        element[i] = Element(_k, _v, _w);
        ele_vec.push_back(Element(_k, _v, _w));
    }
    int opt_type;
    tmpf >> m;
    cerr<<n<<' '<<m<<endl;
    for (int i = 0; i < m; i++) {
        tmpf >> opt_type;
        if (opt_type == 0) {
            tmpf >> _k >> _v>>_w;
            opt_vec.push_back(Opt(0, _k, _v, _w));
        }
        else if (opt_type == 1) {
            tmpf >> _k >> _v >> _w;
            opt_vec.push_back(Opt(1, _k, _v, _w));
        }
        else if (opt_type == 2) {
            tmpf >> _k;
            opt_vec.push_back(Opt(2, _k, 0, 0));
        }
    }

    //aliasTest(n, ele_vec, opt_vec);
//    WIRS_sgt_test(n, element, opt_vec);
//
    if(op==0) WIRS_chunk_test(n, element, opt_vec);
////    WIRS_sgt_bucket_test(n, element, opt_vec,U);
//
    if(op==1)WIRS_chunk_bucket_test(n, element, opt_vec,2.0*csize*(double)U);
//
    if(op==2)WIRS_bstTest(n, element, opt_vec);
//    WIRS_chunk_test(n, element, opt_vec);
//
//    WIRS_chunk_bucket_test(n, element, opt_vec,U*100.0);
//
//    WIRS_bstTest(n, element, opt_vec);
}


int getNum(char *num){
    int x = 0;
    int len = strlen(num);
    for(int i=0;i<len;i++)
        x = x*10 + num[i]-'0';
    return x;
}

int main(int argc,char *argv[]) {
//int main(){
    string filename;
    int n,m,ins_del_ratio,query_ratio,coverage,sample_times,U,op;
    if(argc == 11){
        filename = argv[1];
        n = getNum(argv[2]);
        m = getNum(argv[3]);
        ins_del_ratio = getNum(argv[4]);
        query_ratio = getNum(argv[5]);
        coverage = getNum(argv[6]);
        sample_times = getNum(argv[7]);
        U = getNum(argv[8])+100;
        op = getNum(argv[9]);
        csize = getNum(argv[10]);
        element = new Element[n+100];
        WIRS_produce(n, m, ins_del_ratio, query_ratio,coverage,sample_times,filename);

        string file_s = "/home/DynamicSetSampling/data/WIRS_n" +to_string(n) + "_m" + to_string(m) + "_r1_" +
                        to_string(ins_del_ratio) + "_r2_" + to_string(query_ratio) + "_coverage"+ to_string(coverage)+"_times"+
                        to_string(sample_times);

        WIRS_test_instance(file_s,U,op);
    }
//    filename = "/home/DynamicSetSampling/data/USARoad.data";
//    n = 23947347;
//    m = 10000;
////n = 20;
////m = 1;
//    ins_del_ratio = 100;
//    query_ratio = 0;
//    coverage = 50;
//    sample_times = 100000;
//    U = 10000000;
//    WIRS_produce(n, m, ins_del_ratio, query_ratio,coverage,sample_times,filename);
//    string file_s = "/home/DynamicSetSampling/data/WIRS_n" +to_string(n) + "_m" + to_string(m) + "_r1_" +
//                        to_string(ins_del_ratio) + "_r2_" + to_string(query_ratio) + "_coverage"+ to_string(coverage)+"_times"+
//                        to_string(sample_times);
//
//    WIRS_test_instance(file_s,U);
    return 0;
}
