#include <iostream>
#include <vector>
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
#include "QuickBucket.hpp"
#include "tableBucket.hpp"
using namespace std;

Element USAelement[40000000];
int key_max = 1000000000;
int value_max = 1000;
int SAMPLE_TIMES = 1000000;
// U max weight value, dis_type weight distribution type
void produce_data(int tot_size, int opt_size, int ins_del_ratio, int query_ratio, float U, int dis_type, string file_s) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> key_random(0, tot_size * 2);
    uniform_int_distribution<int> value_random(0, value_max);
    uniform_real_distribution<float> weight_random(1, U);
    if (dis_type == 1) {
        //weight_random;
        // exponential distribution
    }
    uniform_int_distribution<int> opt_random(1, 100);


    vector<Element> elements;
    linearMethod linear_elements;
    cerr << file_s << "\n";
    ofstream tmpf;
    tmpf.open(file_s.c_str());
    tmpf << tot_size << "\n";
    int key, value;
    float weight;
    for (int i = 0; i < tot_size; i++) {
        int key = key_random(gen);
        float weight = weight_random(gen);
        int value = value_random(gen);

        while (linear_elements.find_key(key)) {
            key = key_random(gen);
        }
        //cerr << key << " " << weight << " " << value << "\n";
        //key = i;
        //value = i;
        //weight = i;
        if (i % 10000000 == 0) {
            cerr << i / 10000000 << "\n";
        }
        //if(i> 16777215)
        //cerr <<"this:" << i << " " << key << "\n";
        linear_elements.insert_element(Element(key, value, weight));
        tmpf << key << " " << value << " " << weight << "\n";
    }
    tmpf << opt_size << "\n";
    cout << "sum:" << linear_elements.calc_proportion_sum() << "\n";


    for (int i = 0; i < opt_size; i++) {
        float tmp_rdm = opt_random(gen);
        if (tmp_rdm <= query_ratio) {
            tmpf << "0\n";
            // random sampling query opt
        }
        else {
            tmp_rdm = opt_random(gen);
            if (tmp_rdm <= ins_del_ratio) {
                // insert opt
                int key = key_random(gen);
                int value = value_random(gen);
                float weight = weight_random(gen);
                //cerr << key << "\n";
                while (linear_elements.find_key(key)) {
                    key = key_random(gen);
                }
                linear_elements.insert_element(Element(key, value, weight));
                tmpf << "1 " << key << " " << value << " " << weight << "\n";
            }
            else {
                //delete opt
                int key = linear_elements.random_choose();
                linear_elements.delete_element(key);
                tmpf << "2 " << key << "\n";
            }
        }
        //if()
    }
    cout << "sum2:" << linear_elements.calc_proportion_sum() << "\n";
}

void WIRS_produce_data(int tot_size, int opt_size, int ins_del_ratio, int query_ratio, float U, int dis_type, string file_s) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> key_random(0, key_max);
    uniform_int_distribution<int> value_random(0, value_max);
    uniform_real_distribution<float> weight_random(1, U);
    if (dis_type == 1) {
        weight_random;
        // exponential distribution
    }
    uniform_int_distribution<int> opt_random(1, 100);


    vector<Element> elements;
    linearMethod linear_elements;

    ofstream tmpf;
    tmpf.open(file_s.c_str());
    tmpf << tot_size << "\n";
    for (int i = 0; i < tot_size; i++) {
        int key = key_random(gen);
        float weight = weight_random(gen);
        int value = value_random(gen);
        //value = weight;
        while (linear_elements.find_key(key)) {
            key = key_random(gen);
        }
        linear_elements.insert_element(Element(key, value, weight));
        tmpf << key << " " << value << " " << weight << "\n";
    }
    tmpf << opt_size << "\n";
    cout << "sum:" << linear_elements.calc_proportion_sum() << "\n";


    for (int i = 0; i < opt_size; i++) {
        float tmp_rdm = opt_random(gen);
        if (tmp_rdm <= query_ratio) {
            int l = key_random(gen);
            int r = key_max;
            int t = SAMPLE_TIMES;
            //l = 1, r = key_max;
            tmpf << "0" << " " << l << " " << r << " " << t << "\n";
            // random sampling query opt
        }
        else {
            tmp_rdm = opt_random(gen);
            if (tmp_rdm <= ins_del_ratio) {
                // insert opt
                int key = key_random(gen);
                int value = value_random(gen);
                float weight = weight_random(gen);
                while (linear_elements.find_key(key)) {
                    key = key_random(gen);
                }
                linear_elements.insert_element(Element(key, value, weight));
                tmpf << "1 " << key << " " << value << " " << weight << "\n";
            }
            else {
                //delete opt
                int key = linear_elements.random_choose();
                linear_elements.delete_element(key);
                tmpf << "2 " << key << "\n";
            }
        }
    }
    cout << "sum2:" << linear_elements.calc_proportion_sum() << "\n";
}

void produce(int n, int m, int ins_del_ratio, int query_ratio, int U, int dis_type) {

    string file_s = "/home/fyzhang/DynamicSampling/data/" + to_string(n) + "_m" + to_string(m) + "_r1_" +
        to_string(ins_del_ratio) + "_r2_" + to_string(query_ratio) + "_U" + to_string(U) + "_dis" + to_string(dis_type);
    produce_data(n, m, ins_del_ratio, query_ratio, U, dis_type, file_s);
}

void WIRS_produce(int n, int m, int ins_del_ratio, int query_ratio, int U, int dis_type) {

    string file_s = "/home/fyzhang/DynamicSampling/data/WIRS_n" + to_string(n) + "_m" + to_string(m) + "_r1_" +
        to_string(ins_del_ratio) + "_r2_" + to_string(query_ratio) + "_U" + to_string(U) + "_dis" + to_string(dis_type);
    WIRS_produce_data(n, m, ins_del_ratio, query_ratio, U, dis_type, file_s);
}

void aliasTest(int n, vector<Element>& ele_vec, vector<Opt>& opt_vec) {
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();
    aliasMethod tmp_method(n, ele_vec);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "Alias build_time:" << elapsed_seconds.count() << "\n";
    int tot_query = 0;
    float tot_sum = 0;
    float rvalue = 0;
    start = chrono::system_clock::now();
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            rvalue = tmp_method.random_sample_value();
            //cerr << rvalue << "\n";
            tot_sum += rvalue;
            tot_query++;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            tmp_method.insert_element(tmp_ele);
        }
        if (i.opt_type == 2) {
            tmp_method.delete_element(i.key);
        }
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "Alias Operation_time:" << elapsed_seconds.count() << "\n";

    cout << "alias:" << tot_sum / tot_query << "\n";

    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cerr << "alias_memory :" << gms << "mb\n";
}

void WIRS_sgt_test(int n, Element* ele_vec, vector<Opt>& opt_vec, float U) {
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();

    sort(ele_vec, ele_vec + n, cmp_element_key);
    ScapegoatTree<sizeBlockMethod> scapegoatTree;
    scapegoatTree.init(ele_vec, n);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "sgt tree build_time:" << elapsed_seconds.count() << "\n";

    int tot_query = 0;
    double tot_sum = 0;
    start = chrono::system_clock::now();
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            //cerr << i.key << " " << i.value << " " << i.weight << "\n";
            scapegoatTree.ask(i.key, i.value, i.weight);
            //tot_sum += tmp_method.random_sample_value();
            //tot_query++;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            scapegoatTree.insert(tmp_ele);
            //tmp_method.insert_element(tmp_ele);
        }
        if (i.opt_type == 2) {
            scapegoatTree.erase(i.key);
            //tmp_method.delete_element(i.key);
        }
        //if (i.opt_type == 1)
                //tmp_method.print_blocks();
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "sgt Operation_time:" << elapsed_seconds.count() << "\n";

    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cerr << "sgt_memory " << gms << "mb\n";
    //tot_query = 0;
    //tot_sum = 0;
    //for (int i = 0; i < 1000000; i++) {
    //      tot_sum += tmp_method.random_sample_value();
    //      tot_query++;
    //}

    //~scapegoatTree;
    //cout << "sizeBlock:" << tot_sum / tot_query << "\n";
}

void WIRS_Value_chunk_test(int n, Element* ele_vec, vector<Opt>& opt_vec, float U) {
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();

    sort(ele_vec, ele_vec + n, cmp_element_key);
    ChunkScapegoatTree<valueBucketMethod> scapegoatTree;
    scapegoatTree.setChunkSize(30);
    scapegoatTree.init(ele_vec, n, U);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "sgt value chunk build_time:" << elapsed_seconds.count() << "\n";

    BF tmp_bf;
    tmp_bf.init(ele_vec, n);
    int tot_query = 0;
    double tot_sum = 0;
    start = chrono::system_clock::now();
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            vector<int> ans;
            double sum = 0;
            cerr << i.key << " " << i.value << " " << i.weight << "\n";
            ans = scapegoatTree.ask(i.key, i.value, i.weight);
            //ans.clear();
            for (auto i : ans) {
                sum += i;
            }
            ans.clear();
            cerr << sum / i.weight << " " << tmp_bf.ask_sum(i.key, i.value) << "\n";
            //tot_sum += tmp_method.random_sample_value();
            //tot_query++;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            scapegoatTree.insert(tmp_ele);
            //tmp_method.insert_element(tmp_ele);
        }
        if (i.opt_type == 2) {
            scapegoatTree.erase(i.key);
            //tmp_method.delete_element(i.key);
        }
        //if (i.opt_type == 1)
                //tmp_method.print_blocks();
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "value chunk Operation_time:" << elapsed_seconds.count() << "\n";
    //cout <<"internal sample time " << scapegoatTree.elapsed_seconds.count() << "\n";
    //tot_query = 0;
    //tot_sum = 0;
    //for (int i = 0; i < 1000000; i++) {
    //      tot_sum += tmp_method.random_sample_value();
    //      tot_query++;
    //}
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cerr << "value chunk_memory " << gms << "mb\n";
    //cout << "sizeBlock:" << tot_sum / tot_query << "\n";
}

void WIRS_chunk_test(int n, Element* ele_vec, vector<Opt>& opt_vec, float U) {
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cerr << "chunk_memory " << gms << "mb\n";

    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();

    sort(ele_vec, ele_vec + n, cmp_element_key);
    ChunkScapegoatTree<sizeBlockMethod> scapegoatTree;
    scapegoatTree.setChunkSize(100);
    scapegoatTree.init(ele_vec, n);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "sgt chunk build_time:" << elapsed_seconds.count() << "\n";
    BF tmp_bf;
    tmp_bf.init(ele_vec, n);
    int tot_query = 0;
    double tot_sum = 0;
    start = chrono::system_clock::now();
    //vector<int> tmp_ans;
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            vector<int> ans;
            double sum = 0;
            cerr << i.key << " " << i.value << " " << i.weight << "\n";
            ans = scapegoatTree.ask(i.key, i.value, i.weight);
            //ans.clear();
            for (auto i : ans) {
                sum += i;
            }
            ans.clear();
            cerr << sum / i.weight << " " << tmp_bf.ask_sum(i.key, i.value) << "\n";
            //tot_sum += tmp_method.random_sample_value();
            //tot_query++;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            scapegoatTree.insert(tmp_ele);
            //tmp_method.insert_element(tmp_ele);
        }
        if (i.opt_type == 2) {
            scapegoatTree.erase(i.key);
            //tmp_method.delete_element(i.key);
        }
        //if (i.opt_type == 1)
                //tmp_method.print_blocks();
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "chunk Operation_time:" << elapsed_seconds.count() << "\n";
    //cout <<"internal sample time " << scapegoatTree.elapsed_seconds.count() << "\n";
    //tot_query = 0;
    //tot_sum = 0;
    //for (int i = 0; i < 1000000; i++) {
    //      tot_sum += tmp_method.random_sample_value();
    //      tot_query++;
    //}
    getrusage(RUSAGE_SELF, &rUsage);
    ms = rUsage.ru_maxrss;
    gms = ms / 1024;
    cerr << "chunk_memory " << gms << "mb\n";
    //cout << "sizeBlock:" << tot_sum / tot_query << "\n";
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
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            bst.ask(i.key, i.value, i.weight);
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
        //if (i.opt_type == 1)
                //tmp_method.print_blocks();
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "bst Operation_time:" << elapsed_seconds.count() << "\n";
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cerr << "bst_memory " << gms << "mb\n";
    //tot_query = 0;
    //tot_sum = 0;
    //for (int i = 0; i < 1000000; i++) {
    //      tot_sum += tmp_method.random_sample_value();
    //      tot_query++;
    //}

    //cout << "sizeBlock:" << tot_sum / tot_query << "\n";
}

vector<double> opt_time;
void sizeBlockTest(int n, vector<Element>& ele_vec, vector<Opt>& opt_vec, int m) {
    chrono::time_point<chrono::system_clock> start, end, one_start, one_end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();
    sizeBlockMethod tmp_method(n, ele_vec);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "sizeBlock build_time:" << elapsed_seconds.count() << "\n";

    int tot_query = 0;
    double tot_sum = 0;
    start = chrono::system_clock::now();

    ofstream tmpf;
    int64_t max_opt = 0, min_opt = 999999999;
    tmpf.open("/home/fyzhang/DynamicSampling/data/query_time_per_opt");
    int64_t tot_cost = 0;
    for (auto i : opt_vec) {
        one_start = chrono::high_resolution_clock::now();

        if (i.opt_type == 0) {
            tot_sum += tmp_method.random_sample_value();
            tot_query++;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            tmp_method.insert_element(tmp_ele);
        }
        if (i.opt_type == 2) {
            tmp_method.delete_element(i.key);
        }
        //if(i.opt_type == 1)
        //tmp_method.print_blocks();
        one_end = chrono::high_resolution_clock::now();
        elapsed_seconds = one_end - one_start;

        auto duration = chrono::duration_cast<chrono::nanoseconds>(one_end - one_start);
        tot_cost += duration.count();
        max_opt = max(max_opt, duration.count());
        min_opt = min(min_opt, duration.count());
        opt_time.push_back(duration.count());
        //tmpf << elapsed_seconds.count() << "\n";
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    for (int i = 1; i <= 10000; i++) {
        tot_sum += tmp_method.random_sample_value();
        tot_query++;
    }
    for (auto i : opt_time) {
        tmpf << i << "\n";
    }
    tmpf.close();

    cout << "sizeBlock Operation_time:" << elapsed_seconds.count() << " " << tot_cost << "  " << tot_cost / 1000000000.0 << "\n";
    printf("max: %.10lf min:%.10lf", max_opt * m, min_opt * m);
    cout << "max:" << max_opt << " min:" << min_opt << " tst_cnt" << tmp_method.max_test_cnt << "\n";

    cout << "sizeBlock:" << tot_sum / tot_query << "\n";
}

void valueBucketTest(int n, int U, vector<Element>& ele_vec, vector<Opt>& opt_vec) {

    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;
    start = chrono::system_clock::now();
    valueBucketMethod tmp_method(n, U, ele_vec);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "valueBucket build_time:" << elapsed_seconds.count() << "\n";
    int tot_query = 0;
    float return_value = 0;
    float tot_sum = 0;

    start = chrono::system_clock::now();
    for (auto i : opt_vec) {
        if (i.opt_type == 0) {
            return_value = tmp_method.random_sample_value();
            tot_sum += return_value;
            //cerr << return_value << "\n";
            tot_query++;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            tmp_method.insert_element(tmp_ele);
        }
        if (i.opt_type == 2) {
            tmp_method.delete_element(i.key);
        }
    }

    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "valueBucket Operation_time:" << elapsed_seconds.count() << "\n";

    //tot_query = 0;
    //tot_sum = 0;
    //for (int i = 0; i < 1000000; i++) {
    //      tot_sum += tmp_method.random_sample_value();
    //      tot_query++;
    //}

    cout << tot_sum / tot_query << "\n";
    //cout << "Bucket Generate Random Number:" << tmp_method.random_generate_seconds.count() << "\n";

    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cerr << "value_bucket_memory " << gms << "mb\n";
}



void bstTest(int n, vector<Element>& ele_vec, vector<Opt>& opt_vec) {
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;
    start = chrono::system_clock::now();
    sort(ele_vec.begin(), ele_vec.end(), cmp_element_key);
    bstMethod tmp_method(n, ele_vec);
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "BinarySearchTree build_time:" << elapsed_seconds.count() << "\n";

    int tot_query = 0;
    double tot_sum = 0;
    start = chrono::system_clock::now();
    for (auto i : opt_vec) {
        //cerr << i.opt_type << "\n";
        if (i.opt_type == 0) {
            tot_sum += tmp_method.random_sample_value();
            //cerr << tot_sum;
            tot_query++;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            tmp_method.insert_element(tmp_ele);
        }
        if (i.opt_type == 2) {
            tmp_method.delete_element(i.key);
        }
    }
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;
    cout << "BinarySearchTree Operation_time:" << elapsed_seconds.count() << "\n";

    //tot_query = 0;
    //tot_sum = 0;
    for (int i = 0; i < 1000000; i++) {
          tot_sum += tmp_method.random_sample_value();
          tot_query++;
    }

    cout << "BinarySearchTree:" << tot_sum / tot_query << "\n";
    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    float gms = ms / 1024;
    cerr << "BinarySearchTree_WSS_memory " << gms << "mb\n";
}

void BucketTest(int n, vector<Element>& ele_vec, vector<Opt>& opt_vec, int m) {
    chrono::time_point<chrono::system_clock> start, end, one_start, one_end;
    chrono::duration<double> elapsed_seconds;

    start = chrono::system_clock::now();
    BucketMethod tmp_method(n, ele_vec);
    end = chrono::system_clock::now();
    //tmp_method.print_bucket();

    elapsed_seconds = end - start;
    cout << "Bucket build_time:" << elapsed_seconds.count() << "\n";

    int tot_query = 0;
    double tot_sum = 0;
    start = chrono::system_clock::now();

    double max_opt = 0, min_opt = 999;
    int k = 0;
    for (auto i : opt_vec) {
        //one_start = chrono::system_clock::now();
        k++;
        if (i.opt_type == 0) {
            tot_sum += tmp_method.random_sample_value();
            tot_query++;
        }
        if (i.opt_type == 1) {
            Element tmp_ele = Element(i.key, i.value, i.weight);
            tmp_method.insert_element(tmp_ele);
        }
        if (i.opt_type == 2) {
            tmp_method.delete_element(i.key);
        }
        //if(i.opt_type == 1)
        //tmp_method.print_blocks();
        //one_end = chrono::system_clock::now();
        //elapsed_seconds = one_end - one_start;
        //max_opt = max(max_opt, elapsed_seconds.count());
        //min_opt = min(min_opt, elapsed_seconds.count());
    }
    //tmp_method.print_bucket();
    end = chrono::system_clock::now();
    elapsed_seconds = end - start;


    cout << "Bucket Operation_time:" << elapsed_seconds.count() << "\n";

    cout << "Bucket:" << tot_sum / tot_query <<" tot_query:"<< tot_query << "\n";
}

void test(int n, int m, int ins_del_ratio, int query_ratio, int U, int dis_type) {
    string file_s = "/home/fyzhang/DynamicSampling/data/" + to_string(n) + "_m" + to_string(m) + "_r1_" +
        to_string(ins_del_ratio) + "_r2_" + to_string(query_ratio) + "_U" + to_string(U) + "_dis" + to_string(dis_type);
    cerr << file_s << "\n";
    ifstream tmpf;
    vector<Element> ele_vec;
    vector<Opt> opt_vec;
    tmpf.open(file_s.c_str());
    int _k, _v; float _w;
    tmpf >> n;
    for (int i = 0; i < n; i++) {
        tmpf >> _k >> _v >> _w;
        //if (i == 0)
                //cerr << _k << " " << _v << " " << _w << "\n";
        ele_vec.push_back(Element(_k, _v, _w));
    }
    int opt_type;
    tmpf >> m;

    for (int i = 0; i < m; i++) {
        tmpf >> opt_type;
        if (opt_type == 0) {
            opt_vec.push_back(Opt(0, 0, 0, 0));
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
    sizeBlockTest(n, ele_vec, opt_vec, m);

    //valueBucketTest(n, U, ele_vec, opt_vec);

    bstTest(n, ele_vec, opt_vec);

}

void WSS_test() {
    int n = 10000000;
    int m = 1000000;
    int ins_del_ratio = 0; // percent
    int query_ratio = 100;
    int U = 800000;
    int dis_type = 0;
    produce(n, m, ins_del_ratio, query_ratio, U, dis_type);
    int a[4] = { 0,0,0,0 };
    test(n, m, ins_del_ratio, query_ratio, U, dis_type);
}

void WIRS_test_instance(int n, int m, int ins_del_ratio, int query_ratio, int U, int dis_type) {
    string file_s = "/home/fyzhang/DynamicSampling/data/WIRS_n" + to_string(n) + "_m" + to_string(m) + "_r1_" +
        to_string(ins_del_ratio) + "_r2_" + to_string(query_ratio) + "_U" + to_string(U) + "_dis" + to_string(dis_type);
    ifstream tmpf;
    vector<Element> ele_vec;
    vector<Opt> opt_vec;
    cerr << file_s << "\n";
    tmpf.open(file_s.c_str());
    int _k, _v; float _w;
    tmpf >> n;
    for (int i = 0; i < n; i++) {
        tmpf >> _k >> _v >> _w;
        //cerr << _k <<endl;
        USAelement[i] = Element(_k, _v, _w);
        ele_vec.push_back(Element(_k, _v, _w));
    }
    int opt_type;
    tmpf >> m;

    for (int i = 0; i < m; i++) {
        tmpf >> opt_type;
        if (opt_type == 0) {
            tmpf >> _k >> _v >> _w;
            //_k = 1; _v = key_max; _w = 1e5;
            // _k left interval _v right interval _w sample times (float -> int)
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
    //WIRS_sgt_test(n, element, opt_vec, U);

    //WIRS_chunk_test(n, element, opt_vec, U);
    WIRS_Value_chunk_test(n, USAelement, opt_vec, 1000000000);

    //WIRS_bstTest(n, element, opt_vec);
}

void WIRS_test() {
    int n = 100000;
    int m = 100;
    int ins_del_ratio = 100; // percent
    int query_ratio = 0;
    int U = 1000;
    int dis_type = 0;
    WIRS_produce(n, m, ins_del_ratio, query_ratio, U, dis_type);
    int a[4] = { 0,0,0,0 };
    //vector<Element> tmp;
    //tmp.push_back(Element(1, 1, 1));
    //tmp.push_back(Element(2, 2, 2));
    //tmp.push_back(Element(3, 3, 3));
    //tmp.push_back(Element(4, 4, 4));
    //sizeBlockMethod tmpM(4, tmp);
    //tmpM.delete_element(1);
    //tmpM.delete_element(2);
    //tmpM.delete_element(3);
    //tmpM.delete_element(4);
    //tmpM.insert_element(Element(1, 1, 1));

    //for (int i = 1; i <= 100000; i++) {
    //      a[tmpM.random_sample_value() - 1]++;
    //}
    //cout << a[0] << " " << a[1] << " " << a[2] << " " << a[3];
    WIRS_test_instance(n, m, ins_del_ratio, query_ratio, U, dis_type);
}

int main(int argc, char* args[]) {
    //WIRS_test();
    //WSS_test();
    //exit(0);

    //vector<Element> a;
    //for (int i = 0; i < 1000; i++) {
    //    a.push_back(Element(i, i, i));
    //}
    //ConstantBucket tmp(1000, a);
    //struct rusage rUsage;
    //getrusage(RUSAGE_SELF, &rUsage);
    //long ms = rUsage.ru_maxrss;
    //float gms = ms / 1024;
    //cerr << "table :" << gms << "mb\n";
    //return 0;

    int cnt = 0;
    int algorithm_choose = 5;
    int _file = 1;
    //int action = 0;
    int elements_num;
    int query_num = 100000;
    float max_element = 0;
    string action = "query";
    string src_file = "";
    while (cnt < argc) {
        //cout << args[cnt] << " ";
        if (strcmp(args[cnt], "-algo") == 0) {
            //action = atoi(args[++cnt]);
            algorithm_choose = atoi(args[++cnt]);
        }
        else if (strcmp(args[cnt], "-file") == 0) {
            _file = atoi(args[++cnt]);
            //cout << "file:" << _file << "\n";
        }
        else if (strcmp(args[cnt], "-query") == 0) {
            action = "query";
            query_num = atoi(args[++cnt]);

        }
        else if (strcmp(args[cnt], "-insert") == 0) {
            action = "insert";
            query_num = atoi(args[++cnt]);

        }
        else if (strcmp(args[cnt], "-delete") == 0) {
            action = "delete";
            query_num = atoi(args[++cnt]);

        }
        cnt++;
    }
    vector<Element> ele_vec;
    vector<Opt> opt_vec;
    ifstream elements_file(src_file.c_str());
    robin_hood::unordered_map<int, int> key_choose;

    int _k, _v; float _w;
    cerr << src_file << " " << action << " " << algorithm_choose << " " << query_num << "\n";
    if (action == "insert" || action == "delete") {
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<int> key_random(0, elements_num - 1);
        for (int i = 0; i < query_num; i++) {
            int select_key = key_random(gen);
            while (key_choose.find(select_key) != key_choose.end()) {
                select_key = key_random(gen);
            }
            key_choose[select_key] = 1;
            if (action == "delete") {
                opt_vec.push_back(Opt(2, select_key, 0, 0));
            }
        }
    }
    double tot_weight = 0, tot_sum = 0;
    for (int i = 0; i < elements_num; i++) {
        elements_file >> _k >> _v >> _w;
        tot_weight += _w;
        tot_sum += _v * _w;
        //if (i == 0)
        max_element = max(max_element, _w);
        if (action == "insert" && key_choose.find(i) != key_choose.end()) {
            opt_vec.push_back(Opt(1, _k, _v, _w));
            continue;
        }
        ele_vec.push_back(Element(_k, _v, _w));
        if (algorithm_choose < 0) USAelement[i] = Element(_k, _v, _w);
        //int opt_type;
    }
    cerr<<"prop_res:"<<tot_sum/tot_weight<<"\n";
    auto rng = std::default_random_engine{};
    std::shuffle(std::begin(ele_vec), std::end(ele_vec), rng);
    std::shuffle(std::begin(opt_vec), std::end(opt_vec), rng);
    if (action == "query") {
        for (int i = 0; i < query_num; i++) {
            int l = rand() % ele_vec.size();
            int r = rand() % ele_vec.size();
            if (l > r) swap(l, r);
            opt_vec.push_back(Opt(0, l, r, 1000));
        }
    }

    opt_vec.push_back(Opt(0, 1, 1000, 1000));
    //cerr << src_file <<" " << action << " " << algorithm_choose << " " << query_num << "\n";
    cerr << opt_vec.size() << " " << ele_vec.size() << "\n";
    //opt_vec.push_back(Opt(0, 0, 0, 0));
    //opt_vec.p
    if (algorithm_choose == -1) {
        WIRS_chunk_test(ele_vec.size(), USAelement, opt_vec, 1e9);
    }
    if (algorithm_choose == -2) {
        WIRS_Value_chunk_test(ele_vec.size(), USAelement, opt_vec, 1e9);
    }
    if (algorithm_choose == 1) {
        aliasTest(ele_vec.size(), ele_vec, opt_vec);
    }

    if (algorithm_choose == 2) {
        sizeBlockTest(ele_vec.size(), ele_vec, opt_vec, opt_vec.size());
    }

    if (algorithm_choose == 3) {
        cerr << "max_U" << max_element << "\n";
        valueBucketTest(ele_vec.size(), max_element + 5, ele_vec, opt_vec);
    }

    if (algorithm_choose == 4) {
        cerr << elements_num << " " << opt_vec.size() << "\n";

        bstTest(ele_vec.size(), ele_vec, opt_vec);
    }

    if (algorithm_choose == 5) {
        cerr << elements_num << " " << opt_vec.size() <<" " <<ele_vec[5000].key<<" "<<ele_vec[5000].weight << "\n";
        BucketTest(ele_vec.size(), ele_vec, opt_vec,opt_vec.size());
    }


}