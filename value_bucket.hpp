#pragma once

#include "utility.h"
#include <queue>
#include <random>
#include <iostream>
//#include <unordered_map>
#include <algorithm>
#include <set>
#include <chrono>
#include "XoshiroCpp.hpp"
#include "robin_hood.h"
using namespace std;

class valueBucketMethod {
public:

    XoshiroCpp::Xoroshiro128Plus rng;
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> random_generate_seconds;
    //create the prob_bias and alias table
    string method_name() {
        return "value_bucket";
    }
    /*alias() {}*/
    struct Bucket {
        Bucket() {
            ele_size = 0;
            bucket_weight = 0;
            //max_weight = 0;
        }

        Element& query_element(const int& tmp_key) {
            return elements[position_map[tmp_key]];
        }

        int query_position(const int& tmp_key) {
            return position_map[tmp_key];
        }

        double query_weight(const int& tmp_key) {
            return elements[position_map[tmp_key]].weight;
        }

        bool find_key(const int& tmp_key) {
            return (position_map.find(tmp_key) != position_map.end());
        }

        void insert_element(const Element& ins_ele) {
            //max_weight = max(max_weight, ins_ele.weight);
            ele_size++;
            bucket_weight += ins_ele.weight;
            elements.push_back(ins_ele);
            position_map[ins_ele.key] = elements.size() - 1;
        }

        void delete_element(const int& del_key) {

            int del_pos = position_map[del_key];
            int tail_key = elements[ele_size - 1].key;
            bucket_weight -= elements[del_pos].weight;
            elements[del_pos] = elements[ele_size - 1];
            position_map[tail_key] = del_pos;
            position_map.erase(del_key);
            elements.pop_back();
            ele_size--;
        }

        set<Element>::iterator max_element;

        int ele_size;
        double bucket_weight;
        //float max_weight;
        robin_hood::unordered_map<int, int> position_map;
        vector<Element> elements;
    };


    //set<Element> weight_bst;
    int bucket_size;
    double tot_weight;
    vector<Bucket> buckets;

    int find_bucket(const double & _weight) {
        return log2(_weight);
    }

    valueBucketMethod(int num, vector<Element>& all_ele){}
    valueBucketMethod(int num, float U, vector<Element>& all_ele) : gen(rd()) {
        constexpr std::uint64_t seed = 777;
        rng = XoshiroCpp::Xoroshiro128Plus(seed);
        tot_weight = 0;
        int pow_2 = 2;
        
        //bucket_size = 0;
        for (int i = 0; i < 31; i++) {
            bucket_upper[i] = pow_2;
            if (U <= pow_2) {
                bucket_size = i;
                //cerr << U << " " << bucket_size << "\n";
                break;
            }
            pow_2 *= 2;

        }
        buckets.resize(32);
        //sort(all_ele.begin(), all_ele.end(), cmp_element_weight);
        for (int i = 0; i < num; i++) {
            int bucket_idx = find_bucket(all_ele[i].weight);
            tot_weight += all_ele[i].weight;
            //cerr << bucket_idx<<" "<<i<<" "<<all_ele[i].weight << "\n";
            buckets[bucket_idx].insert_element(all_ele[i]);
        }
        //ele_size = num;
    }
    //randomly generate
    int random_sample_value() {
        uniform_real_distribution<double> dr(0, tot_weight);
        double random_bucket_weight = dr(rng);
        int cur_block = bucket_size;
        //start = chrono::system_clock::now();
        while (random_bucket_weight > buckets[cur_block].bucket_weight) {
            random_bucket_weight -= buckets[cur_block].bucket_weight;
            cur_block--;
        }
        //end = chrono::system_clock::now();
        //random_generate_seconds += end - start;
        uniform_int_distribution<int> pos_rand(0, buckets[cur_block].ele_size - 1);
        uniform_int_distribution<int> dr2(1, int(bucket_upper[cur_block]));
        int tmp_pos;
        int rej_weight;
        while (1) {
            //start = chrono::system_clock::now();
            tmp_pos = pos_rand(rng);
            rej_weight = dr2(rng);
            //end = chrono::system_clock::now();
            //random_generate_seconds += end - start;
            //cerr << tmp_pos <<" "<<cur_block<<" "<<buckets[cur_block].ele_size << "\n";
            //cerr << rej_weight << " " << buckets[cur_block].elements[tmp_pos].weight <<" "<<tmp_pos << "\n";
            if (rej_weight > buckets[cur_block].elements[tmp_pos].weight)
                continue;
            else {
                return buckets[cur_block].elements[tmp_pos].value;
            }
        }
    }


    void delete_element(int del_key) {
        // can be optimized to O(1)
        int cur_block = 0;
        while (!buckets[cur_block].find_key(del_key))
            cur_block++;
        float del_weight = buckets[cur_block].query_weight(del_key);
        buckets[cur_block].delete_element(del_key);
        tot_weight -= del_weight;
        return;
    }

    void insert_element(Element& ins_ele) {
        int cur_block = find_bucket(ins_ele.weight);
        buckets[cur_block].insert_element(ins_ele);
        tot_weight += ins_ele.weight;
        return;
    }
private:
    //int ele_size, reconstruct_flag;
    //int block_capcity[30];
    float bucket_upper[32];
    random_device rd;
    mt19937 gen;
    uniform_int_distribution<int> pos_rand;
    //uniform_real_distribution<float> dr;
};
