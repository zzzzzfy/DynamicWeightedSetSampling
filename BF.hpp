#ifndef DYNAMICSETSAMPLING_BF_HPP
#define DYNAMICSETSAMPLING_BF_HPP
#include<cstdio>
#include<iostream>
#include<vector>
#include<map>
#include "utility.h"

class BF {
private:
    std::vector<Element> elements;
public:
    void init(Element *_elements,int num){
        for (int i = 0; i < num; i++)
            elements.push_back(_elements[i]);
        //elements = std::vector(_elements, _elements + num);
    }
    void insert(Element element){
        elements.push_back(element);
    }
    void erase(int key){
        for(auto it = elements.begin();it!=elements.end();it++)
        {
            if((*it).key==key)
            {
                elements.erase(it);
                return;
            }
        }
    }
    std::map<int,double> ask(int l,int r){
        std::map<int,double> mp;
        double totWeight = 0;
        for(Element element:elements)
            if(l<=element.key&&element.key<=r)
                totWeight+=element.weight;
        for(Element element:elements)
            if(l<=element.key&&element.key<=r)
                mp[element.value]+=element.weight/totWeight;
        return mp;
    }

    double ask_sum(int l, int r) {
        double totWeight = 0;
        double sum = 0;
        for (Element element : elements)
            if (l <= element.key && element.key <= r)
                totWeight += element.weight;
        for (Element element : elements)
            if (l <= element.key && element.key <= r)
                sum += element.value * (element.weight / totWeight);

        return sum;
    }
};


#endif //DYNAMICSETSAMPLING_BF_HPP
