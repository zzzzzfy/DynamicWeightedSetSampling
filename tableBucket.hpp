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

int powT[1000];
vector<vector<int> >sample_map;
//vector<int>sample_map[1000000];
class ConstantBucket {
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

        void clear() {
            position_map.clear();
            elements.clear();
        }

        int ele_size;
        int belong;
        double bucket_weight;
        //double max_weight;
        robin_hood::unordered_map<int, int> position_map;
        vector<Element> elements;
    };
    struct SecondGroup {
        SecondGroup() {}
        ~SecondGroup() { clear(); }
        int log2_ceil(const double& _weight) {
            return ceil(log2(_weight));
        }
        SecondGroup(int _t, int _width, int _upper_value) {
            // loglog{k}
            constexpr std::uint64_t seed = 777;
            rng = XoshiroCpp::Xoroshiro128Plus(seed);
            group_width = _width;
            upper_value = _upper_value;
            status = 0;
            base_pow = pow(2, _t * _width);
            ceil_tot_weight = 0;
            tot_weight = 0;
            second_group_pos_weight.resize(_width + 1, 0);
            ceil_pos_weight.resize(_width + 1, 0);
        }
        void insert_element(int _pos, double _weight) {
            tot_weight += _weight;

            _weight /= base_pow;
            ceil_tot_weight -= ceil_pos_weight[_pos];

            second_group_pos_weight[_pos] += _weight;
            int _tmp_ceil_weight = ceil(second_group_pos_weight[_pos]);

            ceil_pos_weight[_pos] = _tmp_ceil_weight;
            ceil_tot_weight += ceil_pos_weight[_pos];
        }
        void set_up() {
            status = 0;
            int tmp_base = 1;
            for (int i = 0; i < group_width; i++) {
                //                if (ceil_pos_weight[i] > 0) {
                int _tmp_ceil_weight = ceil_pos_weight[i];
                cerr << _tmp_ceil_weight << ' ';
                //                    if(_tmp_ceil_weight<1)
                //                        puts("OK");
                //                    _tmp_ceil_weight--;
                status = status + (_tmp_ceil_weight)*powT[i];
                //                    status = status *tmp_base + _tmp_ceil_weight-1;
                //                }
            }
            cerr << endl;
            cerr << base_pow << " " << status << " status\n";
            cerr << endl;
        }

        void change(int _pos, double _weight) {

            double oldWeight = second_group_pos_weight[_pos];

            tot_weight += _weight - oldWeight * base_pow;

            _weight /= base_pow;
            ceil_tot_weight -= ceil_pos_weight[_pos];

            int tmp = ceil_pos_weight[_pos];
            status = status - tmp * powT[_pos];

            second_group_pos_weight[_pos] = _weight;
            int _tmp_ceil_weight = ceil(second_group_pos_weight[_pos]);

            ceil_pos_weight[_pos] = _tmp_ceil_weight;
            ceil_tot_weight += ceil_pos_weight[_pos];


            tmp = ceil_pos_weight[_pos];
            status = status + tmp * powT[_pos];

        }
        int random_id() {
            uniform_int_distribution<int> pos_rand(1, ceil_tot_weight);
            while (true) {
                int rnd_weight = pos_rand(rng);

                int choose_pos = sample_map[status][rnd_weight];
                uniform_real_distribution<double> dr2(0, ceil_pos_weight[choose_pos]);
                double rej_weight = dr2(rng);
                //cerr<<base_pow<<" "<<choose_pos<<" "<<rnd_weight<<" "<<status<<" "<<rej_weight<<" "<<
                if (rej_weight < second_group_pos_weight[choose_pos]) {
                    return choose_pos;
                }
            }
        }
        void print_group() {
            cerr << base_pow << " " << group_width << " " << upper_value << " "
                << tot_weight << " " << ceil_tot_weight << " " << status << "\n";
        }

        void clear() {
            second_group_pos_weight.clear();
            ceil_pos_weight.clear();
        }
        XoshiroCpp::Xoroshiro128Plus rng;
        double base_pow;
        int group_width;// [log{k}], [log{log{k}}]
        int upper_value;
        int ceil_tot_weight;
        double tot_weight; // before normalize
        vector<double> second_group_pos_weight;
        vector<int> ceil_pos_weight;
        //vector<Bucket> group_buckets;
        int status;
    };
    struct FirstGroup {
        int find_second_idx(const double& _weight) {
            return log2(_weight);
        }
        int log2_ceil(const double& _weight) {
            return ceil(log2(_weight));
        }
        FirstGroup() {}
        ~FirstGroup() { clear(); }
        FirstGroup(int _t, int _width) {
            constexpr std::uint64_t seed = 777;
            rng = XoshiroCpp::Xoroshiro128Plus(seed);
            base_t = _t;
            group_width = _width;
            second_width = log2_ceil(group_width);

            tot_weight = 0;
            //bucket_upper.resize(_width);
            base_pow = pow(2, _t * _width);
            group_pos_weight.resize(_width, 0);
            double tmp_upper = 2;
            //for (int i = 0; i < _width; i++) {
                //bucket_upper[i] = tmp_upper;
                //tmp_upper *= 2;
            //}
        }
        void insert_element(int _pos, double _weight) {
            group_pos_weight[_pos] += _weight;
            tot_weight += _weight;
        }
        int find_first_group_id(int _id) {
            return _id / second_width;
        }
        void set_second_group() {
            check_t = 777;
            for (int i = 0; i < group_width; i++) {
                if (group_pos_weight[i] > 0) {
                    double norm_weight = group_pos_weight[i] / base_pow; // 1< <k^2
                    int _idx = find_second_idx(norm_weight);
                    if (group_buckets.find(_idx) == group_buckets.end()) {
                        group_buckets[_idx] = Bucket();
                    }
                    group_buckets[_idx].insert_element(Element(i, i, norm_weight));
                }
            }

            for (auto i : group_buckets) {
                int group_id = find_first_group_id(i.first);
                int pos_id = i.first - group_id * second_width;
                cerr << i.first << " " << group_id << " " << pos_id << "\n";
                if (second_group.find(group_id) == second_group.end()) {
                    second_group[group_id] = SecondGroup(group_id, second_width,
                        group_width * group_width);
                }
                second_group[group_id].insert_element(pos_id, i.second.bucket_weight);
            }

            for (auto i : second_group) {
                second_group[i.first].set_up();
                //cerr <<"second:" << i.first << "\n";
                //second_group[i.first].print_group();
            }


        }

        void change(int _pos, double _weight) {
            if (group_pos_weight[_pos] > 0) {
                double norm_weight = group_pos_weight[_pos] / base_pow; // 1< <k^2
                tot_weight -= group_pos_weight[_pos];
                int _idx = find_second_idx(norm_weight);
                group_buckets[_idx].delete_element(_pos);

                int group_id = find_first_group_id(_idx);
                int pos_id = _idx - group_id * second_width;

                second_group[group_id].change(pos_id, group_buckets[_idx].bucket_weight);
            }

            group_pos_weight[_pos] = _weight;

            if (_weight > 0) {
                double norm_weight = group_pos_weight[_pos] / base_pow; // 1< <k^2
                int _idx = find_second_idx(norm_weight);
                group_buckets[_idx].insert_element(Element(_pos, _pos, norm_weight));

                int group_id = find_first_group_id(_idx);
                int pos_id = _idx - group_id * second_width;

                tot_weight += group_pos_weight[_pos];

                if (second_group.find(group_id) == second_group.end()) {
                    second_group[group_id] = SecondGroup(group_id, second_width,
                        group_width * group_width);
                    second_group[group_id].insert_element(pos_id, group_buckets[_idx].bucket_weight);
                    second_group[group_id].set_up();
                }
                else {
                    second_group[group_id].change(pos_id, group_buckets[_idx].bucket_weight);
                }
            }

        }

        int find_r() {
            //cerr << log2((tot_weight / base_pow)) << " " << second_width;
            return log2((tot_weight / base_pow)) / second_width;
        }

        int rand_from_group(int r) {
            //cerr << "second_group_random_from_group:" << r <<"\n";
            int choose_bucket = second_group[r].random_id();
            choose_bucket = r * second_width + choose_bucket;
            uniform_int_distribution<int> pos_rand(0, group_buckets[choose_bucket].ele_size - 1);
            uniform_real_distribution<double> dr2(0, pow(2, choose_bucket + 1));
            int tmp_pos;
            double rej_weight;
            while (1) {
                //start = chrono::system_clock::now();
                tmp_pos = pos_rand(rng);
                rej_weight = dr2(rng);
                if (rej_weight > group_buckets[choose_bucket].elements[tmp_pos].weight)
                    continue;
                else {
                    return group_buckets[choose_bucket].elements[tmp_pos].value;
                }
            }
        }

        int rand_pos() {
            double norm_tot_weight = tot_weight / base_pow;
            uniform_real_distribution<double> dr(0, norm_tot_weight);
            double random_group_weight = dr(rng);

            //            double random_group_weight = 0.1;
            int r = find_r();
            //cerr <<"first_group_rand_pos:" << r << " " << random_group_weight 
                //<< " " << tot_weight << " " << base_pow << "\n";
            for (auto i : second_group) {
                //cerr <<
                    //"group_id" << i.first <<" "<<i.second.tot_weight << " group_info:";
                //i.second.print_group();
            }
            double first_weight = 0;
            if (second_group.find(r) != second_group.end()) {
                first_weight = second_group[r].tot_weight;
                //cerr << "entry second group:" << r <<" "<<tot_weight<<" "<<first_weight << "\n";
                //cerr<<
                if (random_group_weight > norm_tot_weight - first_weight) {
                    //sample from G_{r}
                    return rand_from_group(r);
                }
            }
            //cerr << first_weight << "!";
            if (second_group.find(r - 1) != second_group.end()) {
                first_weight += second_group[r - 1].tot_weight;
                //cerr << first_weight << "!";
                if (random_group_weight > norm_tot_weight - first_weight) {
                    return rand_from_group(r - 1);
                }
            }
            if (second_group.find(r - 2) != second_group.end()) {
                first_weight += second_group[r - 2].tot_weight;
                //cerr << first_weight << "!";
                if (random_group_weight > norm_tot_weight - first_weight) {
                    return rand_from_group(r - 2);
                }
            }
            if (second_group.find(r - 3) != second_group.end()) {
                first_weight += second_group[r - 3].tot_weight;
                //cerr << first_weight << "!";
                if (random_group_weight > norm_tot_weight - first_weight) {
                    return rand_from_group(r - 3);
                }
            }
            first_weight = 0;
            for (auto i : second_group) {
                first_weight += i.second.tot_weight;
                if (random_group_weight > norm_tot_weight - first_weight) {
                    // sample from i.second corresponding group
                    return rand_from_group(i.first);
                }
            }
        }

        void clear() {
            group_pos_weight.clear();
            for (auto i : group_buckets) {
                i.second.clear();
            }
            group_buckets.clear();
            for (auto i : second_group) {
                i.second.clear();
            }
            second_group.clear();
        }
        XoshiroCpp::Xoroshiro128Plus rng;

        int base_t; // mean t*[log{k}] .. (t+1)*[log{k}]-1
        double base_pow; // mean 2^{t*[log{k}]}
        int group_width, second_width;// [log{k}], [log{log{k}}]
        double tot_weight;
        //double tot_norm_weight;// use to locate second group
        vector<double> group_pos_weight;
        //vector<Bucket> group_buckets;
        //vector<double> bucket_upper;
        int check_t;
        unordered_map <int, Bucket> group_buckets;
        unordered_map <int, SecondGroup> second_group;
    };

public:
    ConstantBucket() {}

    int pow_int(int x, int y) {
        int res = 1;
        for (int i = 1; i <= y; i++) {
            res *= x;
        }
        return res;
    }
    int log2_floor(const double& _weight) {
        return log2(_weight);
    }
    int log2_ceil(const double& _weight) {
        return ceil(log2(_weight));
    }

    //    void

    void dfs_create_table(int pos, int status, int _sum, vector<int>& status_list, int _width, int upper_value) {
        if (pos == _width) {
            //            _sum+=pos;
            sample_map[status].resize(_sum + 1);
            //            if(status==28)
            //                puts("OK");
            int cur_k = 0;
            int cur_sum = status_list[cur_k];
            for (int i = 1; i <= _sum; i++) {
                while (cur_sum < i) {
                    cur_k++;
                    cur_sum += status_list[cur_k];
                }
                sample_map[status][i] = cur_k;
            }
            return;
        }
        for (int i = 0; i <= upper_value; i++) {
            status_list[pos] = i;
            dfs_create_table(pos + 1, status + i * powT[pos],
                _sum + i, status_list, _width, upper_value);
        }
    }
    void create_table(int _width, int upper_value) {
        powT[0] = 1;
        for (int i = 1; i <= _width; i++)
            powT[i] = powT[i - 1] * (upper_value + 1);
        sample_map.resize(powT[_width]);
        vector<int> status_list;
        status_list.resize(_width + 1, 0);
        dfs_create_table(0, 0, 0, status_list, _width, upper_value);
    }
    int find_first_group_id(int _x) {
        return _x / base_width;
    }

    void rebuild(int N) {
        clear();
        build(N);
    }

    void build(int N) {
        nowSize = N;
        base_width = log2_ceil(nowSize);
        tot_weight = 0;


        second_upper_value = base_width * base_width;
        second_width = log2_ceil(base_width);
        //cerr << "buckets:" << _k << " First_layer" << base_width <<
        //" Second_layer" << second_width << " " << second_upper_value << "\n";
        create_table(second_width, second_upper_value);

        for (int i = 0; i < nowNum; i++) {
            int _idx = log2_floor(elements[i].weight);
            tot_weight += elements[i].weight;
            if (second_bucket.find(_idx) == second_bucket.end()) {
                second_bucket[_idx] = Bucket();
            }
            belongTo[elements[i].key] = _idx;
            second_bucket[_idx].insert_element(elements[i]);
        }


        for (auto i : second_bucket) {
            //            cerr << i.first << " " << i.second.bucket_weight << "\n";
            int group_id = find_first_group_id(i.first);
            int pos_id = i.first - group_id * base_width;
            cerr << i.first << " " << group_id << " " << pos_id << "\n";
            if (first_group.find(group_id) == first_group.end()) {
                first_group[group_id] = FirstGroup(group_id, base_width);
            }
            //cerr << "group:" << group_id << " " << pos_id << "\n";
            first_group[group_id].insert_element(pos_id, i.second.bucket_weight);
        }


        for (auto i : first_group) {
            first_group[i.first].set_second_group();
        }

        //for(int i = )
//        for (auto i : first_group) {
//            //cerr << "Group_" << i.first <<" "<<i.second.check_t << "\n";
//            for (auto j: i.second.second_group) {
//                //cerr <<
//                //"group_id" << j.first << " ";
//                j.second.print_group();
//            }
//        }
    }

    void ins(Element element) {
        elements.push_back(element);
        post[element.key] = nowNum;
        nowNum++;
        tot_weight += element.weight;
        if (nowNum > nowSize)
        {
            rebuild(nowSize * 2);
            return;
        }
        int _idx = log2_floor(element.weight);
        if (second_bucket.find(_idx) == second_bucket.end()) {
            second_bucket[_idx] = Bucket();
            second_bucket[_idx].insert_element(element);
            belongTo[element.key] = _idx;

            int group_id = find_first_group_id(_idx);
            int pos_id = _idx - group_id * base_width;
            if (first_group.find(group_id) == first_group.end()) {
                first_group[group_id] = FirstGroup(group_id, base_width);
                first_group[group_id].insert_element(pos_id, second_bucket[_idx].bucket_weight);
                first_group[group_id].set_second_group();
            }
            else {
                first_group[group_id].change(pos_id, second_bucket[_idx].bucket_weight);
            }
        }
        else {
            second_bucket[_idx].insert_element(element);
            belongTo[element.key] = _idx;

            int group_id = find_first_group_id(_idx);
            int pos_id = _idx - group_id * base_width;
            first_group[group_id].change(pos_id, second_bucket[_idx].bucket_weight);
        }

    }

    void del(int key) {
        int _pos = post[key];
        Element element = elements[_pos];
        nowNum--;
        int tmpKey = elements[nowNum].key;
        post[tmpKey] = _pos;
        swap(elements[nowNum], elements[_pos]);
        elements.pop_back();
        tot_weight -= element.weight;
        if (nowNum <= nowSize / 4) {
            rebuild(nowNum * 2);
            return;
        }

        int _idx = log2_floor(element.weight);

        second_bucket[_idx].delete_element(key);

        int group_id = find_first_group_id(_idx);
        int pos_id = _idx - group_id * base_width;
        first_group[group_id].change(pos_id, second_bucket[_idx].bucket_weight);
    }

    ConstantBucket(int _k, vector<Element>& _elements) {
        elements = _elements;
        for (int i = 0; i < _k; i++)
            post[elements[i].key] = i;
        constexpr std::uint64_t seed = 777;
        rng = XoshiroCpp::Xoroshiro128Plus(seed);
        nowNum = _k;
        build(_k);
    }
    ~ConstantBucket() {
        clear();
    }
    int find_r() {
        return log2(tot_weight) / base_width;
    }
    int random_from_group(int r) {
        //cerr << "first_level_group_r:" << r << "\n";
        int choose_bucket = first_group[r].rand_pos();
        choose_bucket = r * base_width + choose_bucket;

        //cerr <<r<<" choose:" << choose_bucket << " \n";
        uniform_int_distribution<int> pos_rand(0, second_bucket[choose_bucket].ele_size - 1);
        uniform_real_distribution<double> dr2(0, pow(2, choose_bucket + 1));
        int tmp_pos;
        double rej_weight;
        while (1) {
            //start = chrono::system_clock::now();
            tmp_pos = pos_rand(rng);
            rej_weight = dr2(rng);
            if (rej_weight > second_bucket[choose_bucket].elements[tmp_pos].weight)
                continue;
            else {
                return second_bucket[choose_bucket].elements[tmp_pos].value;
            }
        }
    }
    int random_id() {
        uniform_real_distribution<double> dr(0, tot_weight);
        double random_group_weight = dr(rng);
        //cerr << "first_random_id:" << random_group_weight << "\n";
        int r = find_r();
        double first_weight = 0;
        if (first_group.find(r) != first_group.end()) {
            first_weight = first_group[r].tot_weight;
            if (random_group_weight > tot_weight - first_weight) {
                //sample from G_{r}
                return random_from_group(r);
            }
        }
        if (first_group.find(r - 1) != first_group.end()) {
            first_weight += first_group[r - 1].tot_weight;
            if (random_group_weight > tot_weight - first_weight) {
                // sample from G_{r-1}
                return random_from_group(r - 1);
            }
        }
        if (first_group.find(r - 2) != first_group.end()) {
            first_weight += first_group[r - 2].tot_weight;
            if (random_group_weight > tot_weight - first_weight) {
                // sample from G_{r-2}
                return random_from_group(r - 2);
            }
        }
        if (first_group.find(r - 3) != first_group.end()) {
            first_weight += first_group[r - 1].tot_weight;
            if (random_group_weight > tot_weight - first_weight) {
                // sample from G_{r-3}
                return random_from_group(r - 3);
            }
        }
        first_weight = 0;
        for (auto i : first_group) {
            first_weight += i.second.tot_weight;
            if (random_group_weight > tot_weight - first_weight) {
                // sample from i.second corresponding group
                return random_from_group(i.first);
            }
        }
    }


    void clear() {
        for (auto i : second_bucket) i.second.clear();
        for (auto i : first_group) i.second.clear();
        second_bucket.clear();
        first_group.clear();
        for (auto i : sample_map)
            vector<int>().swap(i);
        vector<vector<int> >().swap(sample_map);

    }
    XoshiroCpp::Xoroshiro128Plus rng;

    int num_k, nowSize, nowNum;
    int base_width, second_width;
    int second_upper_value;
    double tot_weight;
    unordered_map <int, Bucket> second_bucket;
    unordered_map <int, FirstGroup> first_group;

    robin_hood::unordered_map<int, int> belongTo;
    robin_hood::unordered_map<int, int> post;
    vector<Element> elements;
    //robin_hood<>;
};