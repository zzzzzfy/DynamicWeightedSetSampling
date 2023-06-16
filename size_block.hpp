#pragma once

#include "utility.h"
#include <queue>
#include <random>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <set>
#include "XoshiroCpp.hpp"
#include "robin_hood.h"
#include "basic_alias.hpp"
using namespace std;

class sizeBlockMethod {
public:
    int max_test_cnt = 0;
    XoshiroCpp::Xoroshiro128Plus rng;
    string method_name() {
        return "size_block";
    }
    struct Block {
        Block() {
            ele_size = 0;
            block_weight = 0;
            //max_weight = 0;
        }

        Element& query_element(const int& tmp_key) {
            return elements[position_map[tmp_key]];
        }

        int query_position(const int& tmp_key) {
            return position_map[tmp_key];
        }

        float query_weight(const int& tmp_key) {
            return elements[position_map[tmp_key]].weight;
        }

        bool find_key(const int& tmp_key) {
            return (position_map.find(tmp_key) != position_map.end());
        }

        void insert_element(const Element& ins_ele) {
            //max_weight = max(max_weight, ins_ele.weight);
            ele_size++;
            block_weight += ins_ele.weight;
            elements.push_back(ins_ele);
            position_map[ins_ele.key] = elements.size() - 1;
        }

        void delete_element(const int& del_key) {

            int del_pos = position_map[del_key];
            int tail_key = elements[ele_size - 1].key;
            block_weight -= elements[del_pos].weight;
            elements[del_pos] = elements[ele_size - 1];
            position_map[tail_key] = del_pos;
            position_map.erase(del_key);
            elements.pop_back();
            ele_size--;
        }

        set<Element>::iterator max_element;

        int ele_size;
        double block_weight;
        //float max_weight;
        robin_hood::unordered_map<int, int> position_map;
        vector<Element> elements;
    };


    set<Element> weight_bst;
    int block_size;
    double tot_weight;
    vector<Block> blocks;
    aliasMethod* block_alias = NULL;
    bool use_alias = false;
    sizeBlockMethod(int num, float U, vector<Element>& all_ele) {}
    sizeBlockMethod(int num, vector<Element>& all_ele) : gen(rd()) {
        cerr << "build block!\n";
        constexpr std::uint64_t seed = 777;
        rng = XoshiroCpp::Xoroshiro128Plus(seed);
        for (int i = 0; i < 31; i++) {
            block_capcity[i] = 1 << i;
        }
        tot_weight = 0;
        blocks.resize(32);
        block_size = 0;
        sort(all_ele.begin(), all_ele.end(), cmp_element_weight);
        for (int i = 0; i < num; i++) {
            tot_weight += all_ele[i].weight;
            if (blocks[block_size].ele_size == block_capcity[block_size]) {
                block_size++;
            }

            blocks[block_size].insert_element(all_ele[i]);
            weight_bst.insert(all_ele[i]);

            if (blocks[block_size].ele_size == 1) {
                // max_element
                blocks[block_size].max_element = weight_bst.find(all_ele[i]);
            }
        }
        if (block_size > alias_threshold) {
            use_alias = true;
        }

        //cerr << blocks[13].max_element->key << "\n";
        ele_size = num;
        check_block_alias();
        //print_blocks();
    }

    void check_block_alias() {
        if (use_alias) {
            if (block_alias != NULL) { delete block_alias; }
            vector<Element> block_elements;
            for (int i = 0; i <= block_size; i++) {
                block_elements.push_back(Element(i, i, blocks[i].block_weight));
                cerr << blocks[i].block_weight << " ";
            }
            block_alias = new aliasMethod(block_elements.size(), block_elements);
        }
    }
    //randomly generate

    int random_sample_value() {
        int test_cnt = 0;
        int cur_block = 0;
        if (use_alias) {
            cur_block = block_alias->random_sample_value();
        }
        else {
            uniform_real_distribution<float> dr(0, tot_weight);
            double random_block_weight = dr(rng);

            while (random_block_weight > blocks[cur_block].block_weight) {
                random_block_weight -= blocks[cur_block].block_weight;
                cur_block++;
            }
            //cerr << cur_block << "\n";
        }
        uniform_int_distribution<int> pos_rand(0, blocks[cur_block].ele_size - 1);
        uniform_real_distribution<float> dr2(0, blocks[cur_block].max_element->weight);
        int tmp_pos;
        float rej_weight;
        while (1) {
            tmp_pos = pos_rand(rng);
            rej_weight = dr2(rng);
            test_cnt++;
            if (rej_weight <= blocks[cur_block].elements[tmp_pos].weight) {
                max_test_cnt = max(max_test_cnt, test_cnt);
                return blocks[cur_block].elements[tmp_pos].value;
            }
        }
    }

    void print_blocks() {
        for (int i = 0; i <= block_size; i++) {
            cout << "block_" << i << " size:" << blocks[i].ele_size << " weight:" << blocks[i].block_weight << " " << blocks[i].max_element->key << " " << blocks[i].max_element->weight << " elements:";
            for (int j = 0; j < blocks[i].ele_size; j++) {
                cout << blocks[i].elements[j].key << ",";
            }
            cout << "\n";
        }
    }


    void delete_element(int del_key) {
        int cur_block = block_size;
        while (!blocks[cur_block].find_key(del_key))
            cur_block--;
        float del_weight = blocks[cur_block].query_weight(del_key);
        tot_weight -= del_weight;



        if (del_key == blocks[cur_block].max_element->key) {
            blocks[cur_block].delete_element(del_key);

            for (int i = cur_block; i < block_size; i++) {
                // move next block max element to this block
                blocks[i].insert_element(*blocks[i + 1].max_element);
            }

            for (int i = cur_block + 1; i <= block_size; i++) {
                // delete max element from this block
                blocks[i].delete_element(blocks[i].max_element->key);
                //blocks[cur_block].max_element++;
            }


            if (blocks[block_size].ele_size == 0) block_size--;
            for (int i = cur_block; i <= block_size; i++) {
                blocks[i].max_element++;
            }
            weight_bst.erase(Element(del_key, 0, del_weight));
        }
        else {
            blocks[cur_block].delete_element(del_key);
            for (int i = cur_block; i < block_size; i++) {
                // move next block max element to this block
                blocks[i].insert_element(*blocks[i + 1].max_element);
            }

            for (int i = cur_block + 1; i <= block_size; i++) {
                // delete max element from this block
                blocks[i].delete_element(blocks[i].max_element->key);
                //blocks[cur_block].max_element++;
            }

            if (blocks[block_size].ele_size == 0) block_size--;
            for (int i = cur_block + 1; i <= block_size; i++) {
                blocks[i].max_element++;
            }
            weight_bst.erase(Element(del_key, 0, del_weight));

        }
        //        for(int i=0;i<=block_size;i++)
        //            if(blocks[i].max_element->weight<1e-10)
        //                puts("No");
        if (block_size == -1) block_size = 0;
        check_block_alias();
        //cerr << "!!!del:" << del_key << " " << cur_block << "\n";
        //print_blocks();
    }

    void insert_element(const Element& ins_ele) {
        int cur_block = block_size;
        //        cerr << cur_block << "\n";
                // a < b means a is larger than b
        while (cur_block > 0 && ins_ele < *blocks[cur_block].max_element)
            cur_block--;
        //        if(cur_block == 3)
        //            puts("OK");
                //cur_block may be -1
        auto ins_ele_it = weight_bst.insert(ins_ele);

        tot_weight += ins_ele.weight;

        blocks[cur_block].insert_element(ins_ele);
        if (blocks[cur_block].ele_size == 1) {
            blocks[cur_block].max_element = ins_ele_it.first;
        }
        for (int i = cur_block; i < block_size; i++) {
            if (blocks[i].ele_size > block_capcity[i]) {
                auto tmp_it = blocks[i + 1].max_element;
                tmp_it--;
                blocks[i].delete_element(tmp_it->key);
                blocks[i + 1].insert_element(*tmp_it);

            }
        }

        for (int i = cur_block; i <= block_size; i++) {
            //cerr << "influenced_block:" << ins_ele.key<<" "<<blocks[c] << "\n";
            if (ins_ele < *blocks[i].max_element) {
                blocks[i].max_element--;
            }
        }

        if (blocks[block_size].ele_size > block_capcity[block_size]) {
            block_size++;
            // insert smallest element
            blocks[block_size].insert_element(blocks[block_size - 1].query_element(weight_bst.rbegin()->key));
            blocks[block_size - 1].delete_element(weight_bst.rbegin()->key);
            blocks[block_size].max_element = weight_bst.end();
            blocks[block_size].max_element--;
        }
        //        for(int i=0;i<=block_size;i++)
        ////            if(blocks[i].max_element->weight<1e-10)
        ////            if(blocks[i].max_element->key == 40593276)
        //            if(blocks[i].ele_size<10)
        //            for(auto x:blocks[i].elements)
        //                if(x.weight>blocks[i].max_element->weight)
        //                puts("No");
        check_block_alias();
        //cerr << "!!!ins:" << ins_ele.key << " " << cur_block << "\n";
        //cerr << tot_weight << "\n";
        //print_blocks();
    }

private:
    int ele_size;
    unsigned int block_capcity[40];

    random_device rd;
    mt19937 gen;
    uniform_int_distribution<int> pos_rand;
    uniform_real_distribution<double> dr;
    //int* alias;
    //double* prob_bias;

};
